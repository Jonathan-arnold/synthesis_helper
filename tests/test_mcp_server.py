"""Tests for the synthesis_helper.mcp.server tools + resources.

We drive the small_cell fixture through the real lazy-loading state
(writing TSVs into a tmp data dir) and then call the tool functions
directly. The FastMCP @tool decorator keeps them as plain callables,
so exercising them this way covers all business logic without async
plumbing.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from synthesis_helper.mcp import server as mcp_server
from synthesis_helper.mcp import state
from synthesis_helper.mcp.server import (
    chemical_resource,
    describe_pathway,
    enumerate_pathways_for,
    find_chemical,
    get_cascade,
    get_shell,
    list_reachables,
    open_cascade_externally,
    open_pathway_externally,
    reachables_resource,
    reaction_resource,
    resynthesize_with_fed,
    stats_resource,
    visualize_cascade,
    visualize_pathway,
)
from tests.fixtures.small_cell import build_small_cell, write_small_cell_tsvs


@pytest.fixture(scope="module")
def small_cell_data_dir(tmp_path_factory) -> Path:
    """Write small_cell TSVs to a fresh tmp dir once per module."""
    sc = build_small_cell()
    d = tmp_path_factory.mktemp("mcp_small_cell_data")
    write_small_cell_tsvs(sc, d)
    return d


@pytest.fixture(autouse=True)
def bootstrap_from_small_cell(small_cell_data_dir: Path):
    """Reset + point state at the small_cell tmp dir before each test."""
    state.reset()
    state.configure(data_dir=small_cell_data_dir)
    yield
    state.reset()


# --- tools -------------------------------------------------------------------


def test_find_chemical_by_exact_name():
    hits = find_chemical("T")
    assert any(h["name"] == "T" and h["id"] == 300 for h in hits)


def test_find_chemical_case_insensitive_substring():
    hits = find_chemical("m_dee")
    names = {h["name"] for h in hits}
    assert "M_deep" in names


def test_get_shell_native():
    assert get_shell("N1")["shell"] == 0


def test_get_shell_target_shortcut():
    """T is reachable via R7 (direct N2 → T) at shell 1 even though R6 is shell 3."""
    assert get_shell("T")["shell"] == 1


def test_get_shell_unreachable_returns_null_shell():
    """M7 is never reachable (only producer needs X_missing)."""
    assert get_shell("M7")["shell"] is None


def test_list_reachables_shell_filter():
    out = list_reachables(max_shell=1)
    assert all(c["shell"] <= 1 for c in out)
    names = {c["name"] for c in out}
    assert {"N1", "N2", "U1"}.issubset(names)
    assert "M_deep" not in names  # shell 3


def test_list_reachables_name_filter():
    out = list_reachables(name_contains="M3")
    names = {c["name"] for c in out}
    assert names <= {"M3a", "M3b"}
    assert names  # at least one


def test_get_cascade_for_target():
    cascade = get_cascade("T")
    assert cascade["target"]["name"] == "T"
    # Every leaf must be shell 0
    assert all(leaf["shell"] == 0 for leaf in cascade["leaf_metabolites"])
    rxn_ids = {r["id"] for r in cascade["reactions"]}
    # Both direct (R7 id=8) and indirect (R6 id=7) routes appear in cascade
    assert 7 in rxn_ids and 8 in rxn_ids


def test_enumerate_pathways_for_target():
    pathways = enumerate_pathways_for("T", max_pathways=5)
    assert pathways, "expected at least one pathway to T"
    # every pathway's final reaction must produce T
    for p in pathways:
        last = p["reactions"][-1]
        assert 300 in last["product_ids"]


def test_describe_pathway_markdown_contains_arrow():
    md = describe_pathway("T", pathway_index=0)
    assert "→" in md
    assert "Pathway 0" in md


# --- visualization tools -----------------------------------------------------

_PNG_MAGIC = b"\x89PNG\r\n\x1a\n"


def _skip_if_no_graphviz():
    pytest.importorskip("graphviz")
    if shutil.which("dot") is None:
        pytest.skip("graphviz `dot` binary not on PATH")


def test_visualize_pathway_returns_png():
    _skip_if_no_graphviz()
    img = visualize_pathway("T", pathway_index=0)
    assert img.data[:8] == _PNG_MAGIC
    assert len(img.data) > 500


def test_visualize_cascade_returns_png():
    _skip_if_no_graphviz()
    img = visualize_cascade("T")
    assert img.data[:8] == _PNG_MAGIC
    assert len(img.data) > 500


def test_visualize_cascade_raises_when_too_large():
    _skip_if_no_graphviz()
    with pytest.raises(ValueError, match="exceeds max_reactions"):
        visualize_cascade("T", max_reactions=0)


# --- externally-opened viewers ----------------------------------------------


@pytest.fixture
def mock_viewer(monkeypatch):
    """Stub subprocess.run inside server so tests never spawn Preview.app."""
    calls: list[list[str]] = []

    def fake_run(cmd, **kwargs):  # noqa: ANN001
        calls.append(list(cmd))

        class _Result:
            returncode = 0

        return _Result()

    monkeypatch.setattr(mcp_server.subprocess, "run", fake_run)
    return calls


def test_open_pathway_externally_writes_png_and_launches_viewer(mock_viewer):
    _skip_if_no_graphviz()
    out = open_pathway_externally("T", pathway_index=0)
    p = Path(out["path"])
    try:
        assert p.exists()
        assert p.read_bytes()[:8] == _PNG_MAGIC
        assert out["opened_in_viewer"] is True
        assert out["chemical"] == "T"
        assert out["steps"] >= 1
        assert out["bytes"] > 500
        # Viewer command was invoked with the saved path
        assert mock_viewer, "expected subprocess.run to be called"
        assert str(p) in mock_viewer[-1]
    finally:
        p.unlink(missing_ok=True)


def test_open_cascade_externally_writes_png_and_launches_viewer(mock_viewer):
    _skip_if_no_graphviz()
    out = open_cascade_externally("T")
    p = Path(out["path"])
    try:
        assert p.exists()
        assert p.read_bytes()[:8] == _PNG_MAGIC
        assert out["opened_in_viewer"] is True
        assert out["reactions"] >= 1
    finally:
        p.unlink(missing_ok=True)


def test_open_pathway_externally_reports_viewer_failure(monkeypatch):
    _skip_if_no_graphviz()

    def fake_run(cmd, **kwargs):  # noqa: ANN001
        raise FileNotFoundError("no such viewer binary")

    monkeypatch.setattr(mcp_server.subprocess, "run", fake_run)
    out = open_pathway_externally("T", pathway_index=0)
    p = Path(out["path"])
    try:
        # PNG should still be written even when the viewer can't launch
        assert p.exists()
        assert out["opened_in_viewer"] is False
        assert out["viewer_error"]
    finally:
        p.unlink(missing_ok=True)


def test_open_pathway_externally_unreachable_raises(mock_viewer):
    _skip_if_no_graphviz()
    with pytest.raises(ValueError, match="not reachable"):
        open_pathway_externally("M7")


def test_open_pathway_externally_index_out_of_range_raises(mock_viewer):
    _skip_if_no_graphviz()
    with pytest.raises(ValueError, match="out of range"):
        open_pathway_externally("T", pathway_index=9999)


# --- reaction label formatting (enzyme name title + EC subtitle) ------------


def test_rxn_label_with_enzyme_name_has_title_and_subtitle():
    from synthesis_helper.mcp.visualize import _rxn_label
    from synthesis_helper.models import Reaction

    rxn = Reaction(id=42, substrates=frozenset(), products=frozenset(), ecnum="1.1.1.1")
    label = _rxn_label(rxn, {"1.1.1.1": "Alcohol dehydrogenase"})
    assert label.startswith("<") and label.endswith(">")
    assert "<B>" in label and "Alcohol dehydrogenase" in label
    assert "EC 1.1.1.1" in label


def test_rxn_label_without_enzyme_name_falls_back_to_ec_as_title():
    from synthesis_helper.mcp.visualize import _rxn_label
    from synthesis_helper.models import Reaction

    rxn = Reaction(id=42, substrates=frozenset(), products=frozenset(), ecnum="1.1.1.1")
    label = _rxn_label(rxn, {})
    assert "EC 1.1.1.1" in label
    assert "rxn #42" in label  # subtitle


def test_rxn_label_no_ec_no_name_uses_rxn_id_only():
    from synthesis_helper.mcp.visualize import _rxn_label
    from synthesis_helper.models import Reaction

    rxn = Reaction(id=42, substrates=frozenset(), products=frozenset())
    label = _rxn_label(rxn, {})
    assert "rxn #42" in label
    assert "<BR/>" not in label  # no subtitle


def test_rxn_label_escapes_html_special_chars_in_name():
    from synthesis_helper.mcp.visualize import _rxn_label
    from synthesis_helper.models import Reaction

    rxn = Reaction(id=1, substrates=frozenset(), products=frozenset(), ecnum="x")
    label = _rxn_label(rxn, {"x": "Weird & <name>"})
    assert "&amp;" in label and "&lt;" in label and "&gt;" in label


def test_ec_names_loaded_from_data_dir(tmp_path):
    """state._bootstrap picks up optional ec_names.tsv."""
    from synthesis_helper.mcp.server import _resolve_or_raise  # noqa: F401
    sc = build_small_cell()
    d = tmp_path / "cell_with_ec"
    d.mkdir()
    write_small_cell_tsvs(sc, d)
    (d / "ec_names.tsv").write_text("ecnum\tname\n1.1.1.1\tAlcohol dehydrogenase\n")

    state.reset()
    state.configure(data_dir=d)
    assert state.get_ec_names() == {"1.1.1.1": "Alcohol dehydrogenase"}
    state.reset()


def test_resynthesize_with_fed_enables_new_chemicals():
    baseline_hg = state.get_hypergraph()
    assert not any(c.name == "M6" for c in baseline_hg.chemical_to_shell)

    out = resynthesize_with_fed(fed_chemical_refs=["F1"])
    assert out["delta_chemicals"] >= 2  # F1 itself + M6
    new_names = {c["name"] for c in out["new_reachables"]}
    assert {"F1", "M6"}.issubset(new_names)


def test_resynthesize_unresolved_refs_are_reported():
    out = resynthesize_with_fed(fed_chemical_refs=["F1", "NO_SUCH_THING"])
    assert "NO_SUCH_THING" in out["unresolved"]


# --- resources ---------------------------------------------------------------


def test_stats_resource_shape():
    s = stats_resource()
    assert s["total_chemicals"] == 15  # small_cell has 15 chems (incl X_missing)
    assert s["reachable_chemicals"] >= 8  # most are reachable
    assert s["shells"] >= 3


def test_reachables_resource_tsv_header():
    tsv = reachables_resource()
    first = tsv.splitlines()[0]
    assert first == "id\tname\tinchi\tshell"
    # N1 (native, shell 0) should appear
    assert any(line.startswith("101\tN1\t") for line in tsv.splitlines())


def test_chemical_resource_produced_and_consumed():
    res = chemical_resource("300")  # T
    assert res["name"] == "T"
    # T is produced by R6 (id=7) and R7 (id=8)
    assert {7, 8}.issubset(set(res["produced_by_reaction_ids"]))


def test_reaction_resource_expands_chemicals():
    res = reaction_resource("8")  # R7: N2 → T
    assert 102 in res["substrate_ids"]
    assert 300 in res["product_ids"]
    sub_names = {s["name"] for s in res["substrates"]}
    prod_names = {p["name"] for p in res["products"]}
    assert sub_names == {"N2"}
    assert prod_names == {"T"}


def test_unknown_chemical_resource_raises():
    with pytest.raises(ValueError):
        chemical_resource("99999999")


def test_unknown_reaction_resource_raises():
    with pytest.raises(ValueError):
        reaction_resource("99999999")
