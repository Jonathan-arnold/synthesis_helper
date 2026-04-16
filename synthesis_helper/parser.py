"""Parsers for MetaCyc input files."""

from __future__ import annotations

from pathlib import Path

from synthesis_helper.models import Chemical, Reaction


def _read_lines(filepath: str | Path) -> list[str]:
    """Read a file, handling both \\n and \\r line endings."""
    text = Path(filepath).read_text(errors="replace")
    # Normalize line endings
    text = text.replace("\r\n", "\n").replace("\r", "\n")
    return [line.strip() for line in text.split("\n")]


def parse_chemicals(filepath: str | Path) -> dict[int, Chemical]:
    """Parse good_chems.txt → {id: Chemical}.

    Expected format (tab-separated, with header):
        id  name  inchi  smiles
    """
    chemicals: dict[int, Chemical] = {}
    for line in _read_lines(filepath):
        if not line or line.startswith("id"):
            continue
        parts = line.split("\t")
        chem_id = int(parts[0])
        name = parts[1] if len(parts) > 1 else ""
        inchi = parts[2].strip('"') if len(parts) > 2 else ""
        smiles = parts[3] if len(parts) > 3 else ""
        chemicals[chem_id] = Chemical(id=chem_id, name=name, inchi=inchi, smiles=smiles)
    return chemicals


def parse_reactions(
    filepath: str | Path, chemicals: dict[int, Chemical]
) -> list[Reaction]:
    """Parse good_reactions.txt → list[Reaction].

    Expected format (tab-separated, with header):
        rxnid  ecnum  substrates(space-separated ids)  products(space-separated ids)
    """
    reactions: list[Reaction] = []
    for line in _read_lines(filepath):
        if not line or line.startswith("rxnid"):
            continue
        parts = line.split("\t")
        rxn_id = int(parts[0])
        ecnum = parts[1] if len(parts) > 1 else ""
        substrate_ids = (
            [int(x) for x in parts[2].split()] if len(parts) > 2 and parts[2].strip() else []
        )
        product_ids = (
            [int(x) for x in parts[3].split()] if len(parts) > 3 and parts[3].strip() else []
        )
        substrates = frozenset(
            chemicals[sid] for sid in substrate_ids if sid in chemicals
        )
        products = frozenset(
            chemicals[pid] for pid in product_ids if pid in chemicals
        )
        reactions.append(
            Reaction(id=rxn_id, substrates=substrates, products=products, ecnum=ecnum)
        )
    return reactions


def parse_metabolite_list(
    filepath: str | Path, chemicals: dict[int, Chemical]
) -> set[Chemical]:
    """Parse a metabolite list file (minimal_metabolites.txt or ubiquitous_metabolites.txt).

    Expected format (tab-separated, with header, possibly \\r line endings):
        name  inchi  descriptor

    Matches metabolites to chemicals by InChI.
    """
    # Build InChI → Chemical lookup (strip quotes for matching)
    inchi_to_chem: dict[str, Chemical] = {}
    name_to_chem: dict[str, Chemical] = {}
    for chem in chemicals.values():
        clean_inchi = chem.inchi.strip('"')
        if clean_inchi:
            inchi_to_chem[clean_inchi] = chem
        if chem.name:
            name_to_chem[chem.name.lower()] = chem

    metabolites: set[Chemical] = set()
    unmatched: list[str] = []

    for line in _read_lines(filepath):
        if not line or line.startswith("name"):
            continue
        parts = line.split("\t")
        name = parts[0] if len(parts) > 0 else ""
        inchi = parts[1].strip('"') if len(parts) > 1 else ""

        # Try matching by InChI first, then by name
        matched = False
        if inchi and inchi in inchi_to_chem:
            metabolites.add(inchi_to_chem[inchi])
            matched = True
        elif name.lower() in name_to_chem:
            metabolites.add(name_to_chem[name.lower()])
            matched = True

        if not matched:
            unmatched.append(name)

    if unmatched:
        print(f"  Warning: {len(unmatched)} metabolites not matched to chemicals: "
              f"{unmatched[:5]}{'...' if len(unmatched) > 5 else ''}")

    return metabolites
