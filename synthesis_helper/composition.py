"""Convert a Pathway into an enzyme Composition."""

from __future__ import annotations

from dataclasses import dataclass, field

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


# ---------------------------------------------------------------------------
# Annotated composition — engineering-relevant metadata per reaction step.
# ---------------------------------------------------------------------------

_EC_CLASSES: dict[str, str] = {
    "1": "Oxidoreductase",
    "2": "Transferase",
    "3": "Hydrolase",
    "4": "Lyase",
    "5": "Isomerase",
    "6": "Ligase",
    "7": "Translocase",
}


def _ec_class(ecnum: str) -> str:
    if not ecnum:
        return "Unknown"
    head = ecnum.split(".", 1)[0]
    return _EC_CLASSES.get(head, "Unknown")


def _is_p450(ecnum: str) -> bool:
    # EC 1.14.x covers oxygenases dominated by cytochrome P450s. Not a
    # precise filter — 1.14.13.x is the classic monooxygenase subclass —
    # but good enough for "this step is probably hard to express".
    return ecnum.startswith("1.14.")


@dataclass
class AnnotatedEnzyme:
    """One composition step annotated for research / engineering decisions."""
    step: int
    ecnum: str
    name: str
    is_orphan: bool  # ecnum blank OR not in enzyme_map
    is_p450: bool    # ecnum starts with 1.14.
    ec_class: str    # Oxidoreductase / Transferase / … / Unknown
    substrate_ids: list[int] = field(default_factory=list)
    product_ids: list[int] = field(default_factory=list)


def annotate_pathway(
    pathway: Pathway,
    enzyme_map: dict[str, str],
) -> list[AnnotatedEnzyme]:
    """Walk the pathway's reactions in order and return annotated enzyme rows.

    Unlike :func:`pathway_to_composition` which returns a minimal Enzyme,
    this returns per-step metadata the LLM / downstream caller can use to
    decide clonability: orphan ECs (no characterised enzyme), cytochrome
    P450s (notoriously hard to express in E. coli), and EC class.
    """
    out: list[AnnotatedEnzyme] = []
    for i, rxn in enumerate(pathway.reactions, start=1):
        ec = rxn.ecnum or ""
        name = enzyme_map.get(ec) if ec else None
        out.append(
            AnnotatedEnzyme(
                step=i,
                ecnum=ec,
                name=name or (f"Unknown_EC_{ec}" if ec else f"rxn_{rxn.id}"),
                is_orphan=(not ec) or (ec not in enzyme_map),
                is_p450=_is_p450(ec),
                ec_class=_ec_class(ec),
                substrate_ids=sorted(s.id for s in rxn.substrates),
                product_ids=sorted(p.id for p in rxn.products),
            )
        )
    return out
