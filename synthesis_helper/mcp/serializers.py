"""JSON-safe DTO conversion for Chemical, Reaction, Cascade, Pathway."""

from __future__ import annotations

from typing import Any

from synthesis_helper.composition import AnnotatedEnzyme
from synthesis_helper.models import Cascade, Chemical, HyperGraph, Pathway, Reaction


def chemical_to_dto(chem: Chemical, hg: HyperGraph | None = None) -> dict[str, Any]:
    return {
        "id": chem.id,
        "name": chem.name,
        "inchi": chem.inchi,
        "smiles": chem.smiles,
        "shell": hg.chemical_to_shell.get(chem) if hg is not None else None,
    }


def reaction_to_dto(rxn: Reaction, hg: HyperGraph | None = None) -> dict[str, Any]:
    return {
        "id": rxn.id,
        "ecnum": rxn.ecnum,
        "substrate_ids": sorted(s.id for s in rxn.substrates),
        "product_ids": sorted(p.id for p in rxn.products),
        "shell": hg.reaction_to_shell.get(rxn) if hg is not None else None,
    }


def cascade_to_dto(cascade: Cascade, hg: HyperGraph) -> dict[str, Any]:
    """Flatten a Cascade for wire transfer.

    Includes the subset of metabolites appearing anywhere in the cascade,
    plus which of those are shell-0 (leaf / native).
    """
    reactions = sorted(cascade.reactions, key=lambda r: (hg.reaction_to_shell.get(r, 0), r.id))
    all_chems: set[Chemical] = set()
    for r in reactions:
        all_chems.update(r.substrates)
        all_chems.update(r.products)

    leaves = [c for c in all_chems if hg.chemical_to_shell.get(c) == 0]
    depth = hg.chemical_to_shell.get(cascade.target, 0)

    return {
        "target": chemical_to_dto(cascade.target, hg),
        "depth": depth,
        "reaction_count": len(reactions),
        "metabolite_count": len(all_chems),
        "reactions": [reaction_to_dto(r, hg) for r in reactions],
        "metabolites": sorted(
            (chemical_to_dto(c, hg) for c in all_chems),
            key=lambda d: (d["shell"] if d["shell"] is not None else -1, d["id"]),
        ),
        "leaf_metabolites": sorted(
            (chemical_to_dto(c, hg) for c in leaves),
            key=lambda d: d["id"],
        ),
    }


def pathway_to_dto(pathway: Pathway, index: int, hg: HyperGraph) -> dict[str, Any]:
    return {
        "index": index,
        "target": chemical_to_dto(pathway.target, hg),
        "length": len(pathway.reactions),
        "reactions": [reaction_to_dto(r, hg) for r in pathway.reactions],
        "metabolites": sorted(
            (chemical_to_dto(c, hg) for c in pathway.metabolites),
            key=lambda d: (d["shell"] if d["shell"] is not None else -1, d["id"]),
        ),
    }


def annotated_enzyme_to_dto(ae: AnnotatedEnzyme) -> dict[str, Any]:
    return {
        "step": ae.step,
        "ecnum": ae.ecnum,
        "enzyme_name": ae.name,
        "is_orphan": ae.is_orphan,
        "is_p450": ae.is_p450,
        "is_heme": ae.is_heme,
        "ec_class": ae.ec_class,
        "substrate_ids": list(ae.substrate_ids),
        "product_ids": list(ae.product_ids),
    }


def similarity_hit_to_dto(
    chem: Chemical, tanimoto: float, hg: HyperGraph | None = None
) -> dict[str, Any]:
    """chemical_to_dto + the Tanimoto similarity score (rounded to 4 dp)."""
    out = chemical_to_dto(chem, hg)
    out["tanimoto"] = round(float(tanimoto), 4)
    return out


def pathway_comparison_row_to_dto(
    pathway_index: int,
    pathway: Pathway,
    metrics: dict[str, Any],
    hg: HyperGraph | None = None,
) -> dict[str, Any]:
    """One row in `compare_pathways`. *metrics* is produced by the tool itself
    and merged into the row (step_count, unique_ecs, cofactor_uses, …)."""
    row: dict[str, Any] = {
        "pathway_index": pathway_index,
        "target": chemical_to_dto(pathway.target, hg),
        "step_count": len(pathway.reactions),
    }
    row.update(metrics)
    return row
