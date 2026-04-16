"""Tests for the BFS waveform expansion."""

from synthesis_helper.models import Chemical, Reaction
from synthesis_helper.synthesize import synthesize
from synthesis_helper.traceback import traceback


def _make_chemical(id: int, name: str = "") -> Chemical:
    return Chemical(id=id, name=name or f"chem_{id}")


def test_simple_expansion():
    """Two-step linear pathway: A → B → C."""
    a = _make_chemical(1, "A")
    b = _make_chemical(2, "B")
    c = _make_chemical(3, "C")

    r1 = Reaction(id=1, substrates=frozenset([a]), products=frozenset([b]))
    r2 = Reaction(id=2, substrates=frozenset([b]), products=frozenset([c]))

    hg = synthesize([r1, r2], native_metabolites={a}, universal_metabolites=set(), verbose=False)

    assert hg.chemical_to_shell[a] == 0
    assert hg.chemical_to_shell[b] == 1
    assert hg.chemical_to_shell[c] == 2
    assert len(hg.reaction_to_shell) == 2


def test_unreachable_reaction():
    """Reaction with unavailable substrate is never enabled."""
    a = _make_chemical(1, "A")
    b = _make_chemical(2, "B")
    x = _make_chemical(99, "X")

    r1 = Reaction(id=1, substrates=frozenset([a, x]), products=frozenset([b]))

    hg = synthesize([r1], native_metabolites={a}, universal_metabolites=set(), verbose=False)

    assert b not in hg.chemical_to_shell
    assert len(hg.reaction_to_shell) == 0


def test_traceback_simple():
    """Traceback finds the reaction that produces B from A."""
    a = _make_chemical(1, "A")
    b = _make_chemical(2, "B")

    r1 = Reaction(id=1, substrates=frozenset([a]), products=frozenset([b]))

    hg = synthesize([r1], native_metabolites={a}, universal_metabolites=set(), verbose=False)
    cascade = traceback(hg, b)

    assert cascade.target == b
    assert r1 in cascade.reactions
