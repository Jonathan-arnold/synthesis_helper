"""Core data model: Chemical, Reaction, HyperGraph, Cascade, Pathway."""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class Chemical:
    id: int
    name: str
    inchi: str = ""
    smiles: str = ""

    def __hash__(self) -> int:
        return hash(self.id)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Chemical):
            return NotImplemented
        return self.id == other.id


@dataclass(frozen=True)
class Reaction:
    id: int
    substrates: frozenset[Chemical]
    products: frozenset[Chemical]
    ecnum: str = ""

    def __hash__(self) -> int:
        return hash(self.id)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Reaction):
            return NotImplemented
        return self.id == other.id


@dataclass
class HyperGraph:
    """Result of BFS waveform expansion."""

    chemical_to_shell: dict[Chemical, int] = field(default_factory=dict)
    reaction_to_shell: dict[Reaction, int] = field(default_factory=dict)
    chemical_to_cascade: dict[Chemical, Cascade] = field(default_factory=dict)


@dataclass
class Cascade:
    """Tree of all reactions that can produce a given chemical."""

    target: Chemical
    reactions: set[Reaction] = field(default_factory=set)


@dataclass
class Pathway:
    """A single route from native metabolites to a target chemical."""

    target: Chemical
    reactions: list[Reaction] = field(default_factory=list)
    metabolites: set[Chemical] = field(default_factory=set)
