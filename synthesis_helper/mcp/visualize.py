"""Graphviz PNG rendering for Cascades and Pathways.

Bipartite KEGG-style layout:
  - Chemicals: ellipse nodes, color-coded by shell (green=native, salmon=target).
  - Reactions: box nodes (light yellow).
  - Edges: substrate -> reaction -> product.
  - Engine: dot (hierarchical, top-down), with rank=same subgraphs per shell
    so chemicals cluster by depth.

Requires the `graphviz` Python package AND the `dot` system binary
(brew install graphviz / apt install graphviz). Import is done lazily inside
each entry point so the MCP server still starts if the dep is missing —
the tool fails at call time with a clear message instead.
"""

from __future__ import annotations

from typing import Iterable

from synthesis_helper.models import Cascade, Chemical, HyperGraph, Pathway, Reaction


_GRAPH_ATTRS = {
    "rankdir": "TB",
    "splines": "spline",
    "nodesep": "0.35",
    "ranksep": "0.55",
    "bgcolor": "white",
    "fontname": "Helvetica",
}
_NODE_DEFAULTS = {"fontname": "Helvetica", "fontsize": "10"}
_EDGE_DEFAULTS = {"arrowsize": "0.7", "color": "gray40"}


def _chem_label(chem: Chemical, shell: int | None) -> str:
    name = chem.name if chem.name else f"#{chem.id}"
    shell_str = "?" if shell is None else str(shell)
    return f"{name}\\n(#{chem.id}, s={shell_str})"


def _rxn_label(rxn: Reaction) -> str:
    return f"EC {rxn.ecnum}" if rxn.ecnum else f"rxn #{rxn.id}"


def _chem_style(chem: Chemical, shell: int | None, is_target: bool) -> dict[str, str]:
    if is_target:
        return {
            "shape": "ellipse",
            "style": "filled,bold",
            "fillcolor": "salmon",
            "color": "firebrick",
            "penwidth": "2",
        }
    if shell == 0:
        return {
            "shape": "ellipse",
            "style": "filled",
            "fillcolor": "lightgreen",
            "color": "darkgreen",
        }
    return {
        "shape": "ellipse",
        "style": "filled",
        "fillcolor": "white",
        "color": "gray40",
    }


def _rxn_style() -> dict[str, str]:
    return {
        "shape": "box",
        "style": "filled,rounded",
        "fillcolor": "lightyellow",
        "color": "gray30",
    }


def _lazy_import_graphviz():
    try:
        import graphviz  # type: ignore
    except ImportError as e:
        raise ImportError(
            "The `graphviz` Python package is required for visualization tools. "
            "Install with `pip install graphviz` and ensure the `dot` system "
            "binary is on PATH (brew install graphviz / apt install graphviz)."
        ) from e
    return graphviz


def _build_digraph(
    reactions: Iterable[Reaction],
    target: Chemical,
    hg: HyperGraph,
):
    graphviz = _lazy_import_graphviz()
    dot = graphviz.Digraph(engine="dot")
    dot.attr(**_GRAPH_ATTRS)
    dot.attr("node", **_NODE_DEFAULTS)
    dot.attr("edge", **_EDGE_DEFAULTS)

    chems_used: dict[Chemical, int | None] = {}
    rxns_used: list[Reaction] = []
    for r in reactions:
        rxns_used.append(r)
        for c in r.substrates:
            chems_used.setdefault(c, hg.chemical_to_shell.get(c))
        for c in r.products:
            chems_used.setdefault(c, hg.chemical_to_shell.get(c))
    chems_used.setdefault(target, hg.chemical_to_shell.get(target))

    # Chemical nodes grouped by shell (rank=same) for layered look.
    by_shell: dict[int | None, list[Chemical]] = {}
    for c, s in chems_used.items():
        by_shell.setdefault(s, []).append(c)
    for shell_val, members in sorted(
        by_shell.items(), key=lambda kv: (kv[0] is None, kv[0] or 0)
    ):
        with dot.subgraph() as sg:
            sg.attr(rank="same")
            for c in sorted(members, key=lambda x: x.id):
                sg.node(
                    f"c{c.id}",
                    _chem_label(c, shell_val),
                    **_chem_style(c, shell_val, is_target=(c == target)),
                )

    for r in rxns_used:
        dot.node(f"r{r.id}", _rxn_label(r), **_rxn_style())
        for s in sorted(r.substrates, key=lambda x: x.id):
            dot.edge(f"c{s.id}", f"r{r.id}")
        for p in sorted(r.products, key=lambda x: x.id):
            dot.edge(f"r{r.id}", f"c{p.id}")

    return dot


def render_pathway_png(pathway: Pathway, hg: HyperGraph) -> bytes:
    """Render a single Pathway to PNG bytes."""
    dot = _build_digraph(pathway.reactions, pathway.target, hg)
    return dot.pipe(format="png")


def render_cascade_png(
    cascade: Cascade,
    hg: HyperGraph,
    max_reactions: int = 80,
) -> bytes:
    """Render a Cascade (full producer tree) to PNG bytes.

    Raises ValueError when the cascade exceeds `max_reactions` — big cascades
    produce unreadable giant images; the caller should lower
    `max_producers_per_chemical` instead.
    """
    n = len(cascade.reactions)
    if n > max_reactions:
        raise ValueError(
            f"Cascade has {n} reactions, exceeds max_reactions={max_reactions}. "
            "Lower max_producers_per_chemical, or raise max_reactions "
            "if you really want a giant diagram."
        )
    dot = _build_digraph(cascade.reactions, cascade.target, hg)
    return dot.pipe(format="png")
