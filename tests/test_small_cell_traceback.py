"""Scenario (d-traceback): Cascade for a multi-producer target must include ALL routes."""

from __future__ import annotations

import pytest

from synthesis_helper.traceback import traceback


def test_cascade_for_T_contains_both_producers(small_cell, small_cell_hg):
    """Cascade(T) must include R7 (direct N2 -> T) AND R6 (M5+M3a -> T) plus
       R6's upstream chain: R4 (makes M5), R3 (makes M4 for R4), R5 (makes M3a).
    """
    c = small_cell.chems
    r = small_cell.rxns
    cascade = traceback(small_cell_hg, c["T"])

    assert cascade.target == c["T"]
    expected = {r["R7"], r["R6"], r["R4"], r["R3"], r["R5"]}
    assert cascade.reactions == expected, (
        f"Expected reactions {sorted(x.id for x in expected)}, "
        f"got {sorted(x.id for x in cascade.reactions)}"
    )


def test_cascade_stops_at_shell_zero(small_cell, small_cell_hg):
    """Cascade(T) must NOT include R1, R2, R2b — they produce chemicals unrelated to T."""
    c = small_cell.chems
    r = small_cell.rxns
    cascade = traceback(small_cell_hg, c["T"])
    for unrelated in ("R1", "R2", "R2b", "R8", "R9"):
        assert r[unrelated] not in cascade.reactions, (
            f"{unrelated} should not be in Cascade(T)"
        )


def test_cascade_for_intermediate_subset(small_cell, small_cell_hg):
    """Cascade(M5) should include only the chain producing M5 and its substrates:
       R4 (direct), R3 (makes M4). It must NOT include R5, R6, R7 — those are
       T-specific or M3a-specific.
    """
    c = small_cell.chems
    r = small_cell.rxns
    cascade = traceback(small_cell_hg, c["M5"])
    assert cascade.reactions == {r["R4"], r["R3"]}


def test_traceback_on_unreachable_raises(small_cell, small_cell_hg):
    """(e-traceback) Target not in the hypergraph raises ValueError."""
    c = small_cell.chems
    with pytest.raises(ValueError):
        traceback(small_cell_hg, c["M7"])


def test_traceback_max_producers_keeps_shortest_route(small_cell, small_cell_hg):
    """cap=1 on T must keep only R7 (shell 1, direct N2 -> T) and drop the longer
       R6 branch (shell 3, via M5+M3a). The producers index sorts by
       (reaction_shell, reaction_id), so shortest routes survive the cap.
    """
    c = small_cell.chems
    r = small_cell.rxns
    cascade = traceback(small_cell_hg, c["T"], max_producers_per_chemical=1)

    assert cascade.reactions == {r["R7"]}, (
        f"cap=1 should retain only the shortest producer R7; got "
        f"{sorted(x.id for x in cascade.reactions)}"
    )


def test_traceback_max_producers_none_matches_default(small_cell, small_cell_hg):
    """cap=None is explicitly equivalent to the default (no cap). Guards against
       accidental regressions where None is mis-treated as 0 or as truthy.
    """
    c = small_cell.chems
    default_cascade = traceback(small_cell_hg, c["T"])
    explicit_none = traceback(small_cell_hg, c["T"], max_producers_per_chemical=None)
    assert default_cascade.reactions == explicit_none.reactions


def test_traceback_max_producers_cap_large_enough_is_full_cascade(
    small_cell, small_cell_hg
):
    """Cap >= max fan-in (T has 2 producers, every other chemical has <=1) must be
       indistinguishable from uncapped. Confirms the cap is an upper bound, not a
       hard truncation.
    """
    c = small_cell.chems
    full = traceback(small_cell_hg, c["T"])
    capped = traceback(small_cell_hg, c["T"], max_producers_per_chemical=10)
    assert full.reactions == capped.reactions
