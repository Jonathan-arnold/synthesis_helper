"""Bootstrap an organism-specific shell 0 from EnzymeMap.

Downloads EnzymeMap's processed_reactions.csv.gz, filters to a chosen organism,
maps each reaction's SMILES onto the existing good_chems.txt inventory (via
InChI), runs the BFS expansion to completion, and writes every reachable
chemical to a file in the ubiquitous_metabolites.txt format so it can be
consumed as shell 0 by main.py.

Usage:
    python bootstrap_ecoli_shell0.py                  # default: coli
    python bootstrap_ecoli_shell0.py --organism coli  # explicit
    python bootstrap_ecoli_shell0.py --organism "Saccharomyces cerevisiae"
"""

from __future__ import annotations

import argparse
import sys
import urllib.request
from pathlib import Path

import pandas as pd

try:
    from rdkit import Chem, RDLogger
    from rdkit.Chem.inchi import MolToInchi
except ImportError:
    sys.stderr.write("rdkit is required. Install with: pip install rdkit\n")
    raise

from synthesis_helper.models import Chemical, Reaction
from synthesis_helper.parser import parse_chemicals, parse_metabolite_list
from synthesis_helper.synthesize import synthesize

RDLogger.DisableLog("rdApp.*")

DATA_DIR = Path(__file__).parent / "data"
ENZYMEMAP_URL = (
    "https://github.com/hesther/enzymemap/raw/main/data/processed_reactions.csv.gz"
)
ENZYMEMAP_PATH = DATA_DIR / "enzymemap_processed_reactions.csv.gz"


def ensure_enzymemap_downloaded(path: Path) -> None:
    if path.exists():
        print(f"EnzymeMap CSV already at {path} ({path.stat().st_size / 1e6:.1f} MB)")
        return
    print(f"Downloading EnzymeMap CSV to {path} ...")
    path.parent.mkdir(parents=True, exist_ok=True)
    urllib.request.urlretrieve(ENZYMEMAP_URL, path)
    print(f"  done ({path.stat().st_size / 1e6:.1f} MB)")


def smiles_to_inchi(smiles: str, cache: dict[str, str | None]) -> str | None:
    if smiles in cache:
        return cache[smiles]
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        cache[smiles] = None
        return None
    inchi = MolToInchi(mol)
    cache[smiles] = inchi if inchi.startswith("InChI=") else None
    return cache[smiles]


def split_reaction_smiles(rxn_smiles: str) -> tuple[list[str], list[str]]:
    """'A.B>>C.D' -> (['A','B'], ['C','D']). Agents (middle of A>B>C) discarded."""
    if ">>" in rxn_smiles:
        left, right = rxn_smiles.split(">>", 1)
    else:
        parts = rxn_smiles.split(">")
        if len(parts) != 3:
            return [], []
        left, _, right = parts
    return [s for s in left.split(".") if s], [s for s in right.split(".") if s]


def build_reactions(
    df: pd.DataFrame,
    chemicals_by_inchi: dict[str, Chemical],
) -> list[Reaction]:
    """Convert EnzymeMap rows into Reaction objects matched to known chemicals.

    Skips any reaction whose molecules can't all be resolved to existing
    Chemical entries — tracking those keeps the reaction graph consistent with
    the rest of the corpus.
    """
    smiles_cache: dict[str, str | None] = {}
    reactions: list[Reaction] = []
    unmatched_smiles: dict[str, int] = {}
    bad_rxn_smiles = 0
    skipped_unmatched = 0
    seen_keys: set[tuple[frozenset[int], frozenset[int]]] = set()

    for next_id, (_, row) in enumerate(df.iterrows(), start=1):
        rxn_smiles = row["unmapped"]
        if not isinstance(rxn_smiles, str):
            bad_rxn_smiles += 1
            continue
        reactant_smi, product_smi = split_reaction_smiles(rxn_smiles)
        if not reactant_smi or not product_smi:
            bad_rxn_smiles += 1
            continue

        substrates: set[Chemical] = set()
        products: set[Chemical] = set()
        ok = True
        for smi_list, bucket in ((reactant_smi, substrates), (product_smi, products)):
            for smi in smi_list:
                inchi = smiles_to_inchi(smi, smiles_cache)
                chem = chemicals_by_inchi.get(inchi) if inchi else None
                if chem is None:
                    ok = False
                    unmatched_smiles[smi] = unmatched_smiles.get(smi, 0) + 1
                    break
                bucket.add(chem)
            if not ok:
                break
        if not ok or not substrates or not products:
            skipped_unmatched += 1
            continue

        key = (
            frozenset(c.id for c in substrates),
            frozenset(c.id for c in products),
        )
        if key in seen_keys:
            continue
        seen_keys.add(key)

        ecnum = row.get("ec_num", "")
        reactions.append(
            Reaction(
                id=next_id,
                substrates=frozenset(substrates),
                products=frozenset(products),
                ecnum=str(ecnum) if pd.notna(ecnum) else "",
            )
        )

    print(f"  Built {len(reactions)} unique reactions.")
    print(f"  Skipped {skipped_unmatched} rows with unmatched molecules.")
    print(f"  Skipped {bad_rxn_smiles} rows with bad/missing reaction SMILES.")
    print(f"  Unique unmatched SMILES: {len(unmatched_smiles)}")
    if unmatched_smiles:
        top = sorted(unmatched_smiles.items(), key=lambda kv: -kv[1])[:10]
        print("  Top unmatched (count, smiles):")
        for smi, count in top:
            print(f"    {count:>5}  {smi}")

    return reactions


def write_metabolite_list(
    reachables: list[Chemical],
    filepath: Path,
    descriptor: str,
) -> None:
    """Write in ubiquitous_metabolites.txt format: name\\tinchi\\tdescriptor."""
    with open(filepath, "w") as f:
        f.write("name\tinchi\tdescriptor\n")
        for chem in reachables:
            inchi = chem.inchi.strip('"')
            if "," in inchi or " " in inchi:
                inchi = f'"{inchi}"'
            f.write(f"{chem.name}\t{inchi}\t{descriptor}\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--organism",
        default="coli",
        help="Case-insensitive substring matched against EnzymeMap's `organism` "
        "column (default: coli).",
    )
    ap.add_argument(
        "--output",
        default=str(DATA_DIR / "ecoli_reachables_shell0.txt"),
        help="Output metabolite-list file.",
    )
    ap.add_argument(
        "--descriptor",
        default="ecoli_reachable",
        help="Value for the `descriptor` column in the output file.",
    )
    args = ap.parse_args()

    ensure_enzymemap_downloaded(ENZYMEMAP_PATH)

    print("Loading chemicals and metabolites ...")
    chemicals = parse_chemicals(DATA_DIR / "good_chems.txt")
    native = parse_metabolite_list(DATA_DIR / "minimal_metabolites.txt", chemicals)
    universal = parse_metabolite_list(DATA_DIR / "ubiquitous_metabolites.txt", chemicals)

    chemicals_by_inchi: dict[str, Chemical] = {}
    for chem in chemicals.values():
        clean = chem.inchi.strip('"')
        if clean:
            chemicals_by_inchi[clean] = chem
    print(
        f"  {len(chemicals)} chemicals, {len(chemicals_by_inchi)} with InChI; "
        f"{len(native)} minimal, {len(universal)} ubiquitous."
    )

    print(f"Loading EnzymeMap from {ENZYMEMAP_PATH} ...")
    df = pd.read_csv(ENZYMEMAP_PATH, compression="gzip", low_memory=False)
    print(f"  {len(df)} rows, columns: {list(df.columns)}")

    org_series = df["organism"].fillna("").astype(str)
    mask = org_series.str.contains(args.organism, case=False, regex=False)
    df_org = df[mask].copy()
    print(f"  {len(df_org)} rows matching organism filter '{args.organism}'.")
    if df_org.empty:
        sys.stderr.write("No rows matched — aborting.\n")
        sys.exit(1)

    print("Mapping reactions to existing chemicals ...")
    reactions = build_reactions(df_org, chemicals_by_inchi)
    if not reactions:
        sys.stderr.write("No reactions could be mapped — aborting.\n")
        sys.exit(1)

    print(f"Running BFS expansion on {len(reactions)} reactions ...")
    hg = synthesize(reactions, native, universal, verbose=True)

    reachables = sorted(hg.chemical_to_shell.keys(), key=lambda c: c.id)
    out_path = Path(args.output)
    print(f"Writing {len(reachables)} reachables to {out_path} ...")
    write_metabolite_list(reachables, out_path, args.descriptor)
    print("Done.")


if __name__ == "__main__":
    main()
