"""Convert a Pathway into an enzyme Composition."""

from __future__ import annotations

from dataclasses import dataclass

from synthesis_helper.models import Pathway


@dataclass
class Enzyme:
    name: str
    ecnum: str = ""


def pathway_to_composition(
    pathway: Pathway,
    enzyme_map: dict[str, str],
) -> list[Enzyme]:
    """Map each reaction in a pathway to an enzyme.

    Args:
        pathway: The pathway to convert.
        enzyme_map: Mapping from EC number → enzyme name.

    Returns:
        Ordered list of Enzymes (one per reaction step).
    """
    composition: list[Enzyme] = []
    for rxn in pathway.reactions:
        name = enzyme_map.get(rxn.ecnum, f"Unknown_EC_{rxn.ecnum or rxn.id}")
        composition.append(Enzyme(name=name, ecnum=rxn.ecnum))
    return composition
