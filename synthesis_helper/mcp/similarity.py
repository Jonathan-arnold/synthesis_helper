"""Morgan-fingerprint based structural similarity over the chemical corpus.

Isolated from state.py so the RDKit import cost (~200-400 ms cold, then
~3-5 s to fingerprint 9361 molecules) stays strictly behind the first call
to `find_by_structure`. Baseline MCP startup never pays for it.
"""

from __future__ import annotations

import time
from typing import Iterable

from synthesis_helper.models import Chemical


_RADIUS = 2
_NBITS = 2048


def _lazy_rdkit():
    """Import RDKit on first call; silence parser warnings globally."""
    from rdkit import Chem, RDLogger  # noqa: F401
    from rdkit.Chem import AllChem

    RDLogger.DisableLog("rdApp.*")
    return Chem, AllChem


def _mol_from_chem(chem: Chemical, Chem) -> "object | None":
    """Prefer SMILES; fall back to InChI. Returns None if neither parses.

    The synthesis_helper corpus uses literal "null" / "" for missing
    SMILES, and both placeholders must be treated as absent before
    handing anything to RDKit.
    """
    smiles = chem.smiles
    if smiles and smiles != "null":
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            return mol
    inchi = chem.inchi
    if inchi and inchi != "null":
        mol = Chem.MolFromInchi(inchi)
        if mol is not None:
            return mol
    return None


def fingerprint_from_smiles(smiles: str):
    """Return a Morgan ExplicitBitVect for *smiles*. Raises ValueError on parse."""
    Chem, AllChem = _lazy_rdkit()
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles!r}")
    return AllChem.GetMorganFingerprintAsBitVect(mol, _RADIUS, nBits=_NBITS)


def build_index(chemicals: Iterable[Chemical]) -> tuple[list[tuple[int, object]], float, int]:
    """Build the corpus-wide Morgan fingerprint index.

    Returns ``(index, elapsed_ms, skipped)`` where ``index`` is a list of
    ``(chem_id, fingerprint)`` tuples, and ``skipped`` counts chemicals for
    which both SMILES and InChI failed to parse in RDKit.
    """
    Chem, AllChem = _lazy_rdkit()
    t0 = time.perf_counter()
    index: list[tuple[int, object]] = []
    skipped = 0
    for chem in chemicals:
        mol = _mol_from_chem(chem, Chem)
        if mol is None:
            skipped += 1
            continue
        try:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, _RADIUS, nBits=_NBITS)
        except Exception:
            skipped += 1
            continue
        index.append((chem.id, fp))
    elapsed_ms = (time.perf_counter() - t0) * 1000.0
    return index, elapsed_ms, skipped


def tanimoto_search(
    query_fp,
    index: list[tuple[int, object]],
    threshold: float,
    limit: int,
) -> list[tuple[int, float]]:
    """Return top-*limit* ``(chem_id, tanimoto)`` from *index* that meet *threshold*.

    Sorted by Tanimoto descending, then by chem_id ascending for determinism.
    """
    from rdkit import DataStructs

    scored: list[tuple[int, float]] = []
    for chem_id, fp in index:
        score = DataStructs.TanimotoSimilarity(query_fp, fp)
        if score >= threshold:
            scored.append((chem_id, float(score)))
    scored.sort(key=lambda t: (-t[1], t[0]))
    return scored[:limit]
