"""Scenarios (a), (b), (c), (e), (f), (g) — BFS waveform expansion correctness."""

from __future__ import annotations

from synthesis_helper.synthesize import synthesize


def test_linear_chain_shells(small_cell, small_cell_hg):
    """(a) Shell propagates correctly across a 3-step linear chain N1 -> M1 -> M2 -> M_deep."""
    c = small_cell.chems
    assert small_cell_hg.chemical_to_shell[c["N1"]] == 0
    assert small_cell_hg.chemical_to_shell[c["M1"]] == 1
    assert small_cell_hg.chemical_to_shell[c["M2"]] == 2
    assert small_cell_hg.chemical_to_shell[c["M_deep"]] == 3


def test_multi_substrate_shell(small_cell, small_cell_hg):
    """(b) R3: {N1, N2} -> {M4}. Both natives (shell 0) -> M4 lands at shell 1."""
    c = small_cell.chems
    assert small_cell_hg.chemical_to_shell[c["M4"]] == 1


def test_universal_does_not_inflate_shell(small_cell, small_cell_hg):
    """(c) R4: {M4 (shell 1), U1 (shell 0)} -> M5 at shell 2 (max+1, U1 is transparent)."""
    c = small_cell.chems
    assert small_cell_hg.chemical_to_shell[c["U1"]] == 0
    assert small_cell_hg.chemical_to_shell[c["M5"]] == 2


def test_unreachable_chemical_absent(small_cell, small_cell_hg):
    """(e) R11: {X_missing} -> {M7}. X_missing is never introduced, so R11 never fires."""
    c = small_cell.chems
    assert c["X_missing"] not in small_cell_hg.chemical_to_shell
    assert c["M7"] not in small_cell_hg.chemical_to_shell
    assert small_cell.rxns["R11"] not in small_cell_hg.reaction_to_shell


def test_tied_shells_bfs_picks_shorter(small_cell, small_cell_hg):
    """(f) M3b has two producers:
          R8: {N1} -> {M3b}       (enables at shell 1)
          R9: {M1} -> {M3b}       (enables at shell 2)
       BFS must assign M3b the shorter shell (1).
    """
    c = small_cell.chems
    assert small_cell_hg.chemical_to_shell[c["M3b"]] == 1


def test_fed_chemical_enables_reaction(small_cell):
    """(g) F1 is not in the base native set, so R10 is disabled and M6 unreachable.
       When F1 is added to natives, R10 fires at shell 1 and M6 appears at shell 1.
    """
    c = small_cell.chems

    hg_base = synthesize(
        small_cell.all_reactions,
        native_metabolites=small_cell.natives,
        universal_metabolites=small_cell.universals,
        verbose=False,
    )
    assert c["M6"] not in hg_base.chemical_to_shell
    assert small_cell.rxns["R10"] not in hg_base.reaction_to_shell

    hg_fed = synthesize(
        small_cell.all_reactions,
        native_metabolites=small_cell.natives | small_cell.fed,
        universal_metabolites=small_cell.universals,
        verbose=False,
    )
    assert hg_fed.chemical_to_shell[c["F1"]] == 0
    assert hg_fed.chemical_to_shell[c["M6"]] == 1


def test_reaction_shells(small_cell, small_cell_hg):
    """Every enabled reaction lands at the shell predicted by ground-truth table."""
    r = small_cell.rxns
    expected = {
        "R1": 1, "R3": 1, "R5": 1, "R7": 1, "R8": 1,
        "R2": 2, "R4": 2, "R9": 2,
        "R2b": 3, "R6": 3,
    }
    for name, shell in expected.items():
        actual = small_cell_hg.reaction_to_shell.get(r[name])
        assert actual == shell, f"{name}: expected shell {shell}, got {actual}"
    assert r["R10"] not in small_cell_hg.reaction_to_shell
    assert r["R11"] not in small_cell_hg.reaction_to_shell


def test_target_reaches_shell_one_via_shorter_route(small_cell, small_cell_hg):
    """T has two producers: R7 (shell 1, direct from N2) and R6 (shell 3, via M5+M3a).
       BFS picks the shorter — T should be shell 1.
       R6 still enables (at shell 3), which is what keeps it available for the cascade.
    """
    c = small_cell.chems
    r = small_cell.rxns
    assert small_cell_hg.chemical_to_shell[c["T"]] == 1
    assert small_cell_hg.reaction_to_shell[r["R6"]] == 3
    assert small_cell_hg.reaction_to_shell[r["R7"]] == 1
