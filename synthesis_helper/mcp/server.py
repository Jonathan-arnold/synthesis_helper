"""FastMCP server exposing synthesize / traceback / pathways to Claude Code."""

from __future__ import annotations

import os
from typing import Any

from mcp.server.fastmcp import FastMCP
from mcp.server.fastmcp.utilities.types import Image

from synthesis_helper.models import Chemical
from synthesis_helper.pathways import enumerate_pathways
from synthesis_helper.traceback import traceback as build_cascade
from synthesis_helper.mcp import lookup, state
from synthesis_helper.mcp.serializers import (
    cascade_to_dto,
    chemical_to_dto,
    pathway_to_dto,
    reaction_to_dto,
)
from synthesis_helper.mcp.visualize import render_cascade_png, render_pathway_png


mcp = FastMCP("synthesis-helper")


_MAX_PATHWAYS_CAP = int(os.environ.get("SYNTHESIS_MAX_PATHWAYS", "1000"))


def _resolve_or_raise(chemical_ref: str | int) -> Chemical:
    chems = state.get_chemicals()
    hit = lookup.resolve_one(chemical_ref, chems)
    if hit is None:
        raise ValueError(
            f"No chemical matches reference {chemical_ref!r}. "
            "Try an id, an exact name, or an InChI string."
        )
    return hit


# -----------------------------------------------------------------------------
# Tools
# -----------------------------------------------------------------------------


@mcp.tool()
def find_chemical(query: str, limit: int = 10) -> list[dict[str, Any]]:
    """Search chemicals by id, name (substring), or InChI.

    Returns up to *limit* matches as dicts with id, name, inchi, smiles, shell.
    Shell is null when the chemical is not reachable from the baseline cell.
    """
    hg = state.get_hypergraph()
    chems = state.get_chemicals()
    hits = lookup.resolve(query, chems, limit=limit)
    return [chemical_to_dto(c, hg) for c in hits]


@mcp.tool()
def get_shell(chemical_ref: str | int) -> dict[str, Any]:
    """Return the shell (min reaction steps from native metabolites) of a chemical.

    shell=0 means native/universal; shell=null means not reachable.
    """
    hg = state.get_hypergraph()
    chem = _resolve_or_raise(chemical_ref)
    return {
        "id": chem.id,
        "name": chem.name,
        "shell": hg.chemical_to_shell.get(chem),
    }


@mcp.tool()
def list_reachables(
    max_shell: int | None = None,
    name_contains: str | None = None,
    limit: int = 100,
) -> list[dict[str, Any]]:
    """List reachable chemicals in the baseline hypergraph, filtered.

    Sorted by (shell, id). Default cap is 100 entries — bump *limit* for more.
    """
    hg = state.get_hypergraph()
    needle = name_contains.lower() if name_contains else None
    out = []
    for chem, shell in sorted(
        hg.chemical_to_shell.items(), key=lambda kv: (kv[1], kv[0].id)
    ):
        if max_shell is not None and shell > max_shell:
            continue
        if needle and needle not in chem.name.lower():
            continue
        out.append(chemical_to_dto(chem, hg))
        if len(out) >= limit:
            break
    return out


@mcp.tool()
def get_cascade(
    chemical_ref: str | int,
    max_producers_per_chemical: int = 5,
) -> dict[str, Any]:
    """Build the full cascade (all reactions reaching a target) from the baseline hypergraph.

    The cascade is the unflattened tree: every reaction that can produce the target,
    recursively, with shell-0 metabolites as leaves. Returned DTO includes
    reaction list, all metabolites touched, and the subset of leaf (shell-0)
    metabolites.

    *max_producers_per_chemical* caps branching when a non-native chemical has
    many producers (shortest-shell producers kept first).
    """
    hg = state.get_hypergraph()
    chem = _resolve_or_raise(chemical_ref)
    if chem not in hg.chemical_to_shell:
        raise ValueError(
            f"{chem.name!r} (id={chem.id}) is not reachable from the baseline cell."
        )
    cascade = build_cascade(hg, chem, max_producers_per_chemical=max_producers_per_chemical)
    return cascade_to_dto(cascade, hg)


@mcp.tool()
def enumerate_pathways_for(
    chemical_ref: str | int,
    max_pathways: int = 10,
    max_producers_per_chemical: int = 5,
) -> list[dict[str, Any]]:
    """Return up to *max_pathways* individual routes from native metabolites to target.

    Each pathway is a topologically ordered list of reactions whose combined
    output produces the target, starting entirely from shell-0 inputs.
    """
    capped = min(max_pathways, _MAX_PATHWAYS_CAP)
    hg = state.get_hypergraph()
    chem = _resolve_or_raise(chemical_ref)
    if chem not in hg.chemical_to_shell:
        raise ValueError(
            f"{chem.name!r} (id={chem.id}) is not reachable from the baseline cell."
        )
    cascade = build_cascade(hg, chem, max_producers_per_chemical=max_producers_per_chemical)
    pathways = enumerate_pathways(cascade, hg, max_pathways=capped)
    return [pathway_to_dto(p, i, hg) for i, p in enumerate(pathways)]


@mcp.tool()
def describe_pathway(
    chemical_ref: str | int,
    pathway_index: int,
    max_producers_per_chemical: int = 5,
) -> str:
    """Render a human-readable Markdown trace of a single pathway.

    Re-enumerates pathways (baseline hypergraph only) and picks index.
    """
    hg = state.get_hypergraph()
    chems = state.get_chemicals()
    chem = _resolve_or_raise(chemical_ref)
    cascade = build_cascade(hg, chem, max_producers_per_chemical=max_producers_per_chemical)
    pathways = enumerate_pathways(cascade, hg, max_pathways=pathway_index + 1)
    if pathway_index >= len(pathways):
        raise ValueError(
            f"Only {len(pathways)} pathway(s) found; index {pathway_index} out of range."
        )
    p = pathways[pathway_index]

    lines = [
        f"# Pathway {pathway_index} → {chem.name} (id={chem.id}, shell={hg.chemical_to_shell.get(chem)})",
        f"*{len(p.reactions)} reaction step(s)*",
        "",
    ]
    for step, rxn in enumerate(p.reactions, 1):
        subs = " + ".join(
            f"{chems[s.id].name}(#{s.id})" if s.id in chems else f"#{s.id}"
            for s in sorted(rxn.substrates, key=lambda c: c.id)
        )
        prods = " + ".join(
            f"{chems[p_.id].name}(#{p_.id})" if p_.id in chems else f"#{p_.id}"
            for p_ in sorted(rxn.products, key=lambda c: c.id)
        )
        ec = f"EC {rxn.ecnum}" if rxn.ecnum else f"rxn #{rxn.id}"
        lines.append(f"{step}. {subs} → {prods}  *({ec})*")
    return "\n".join(lines)


@mcp.tool()
def visualize_pathway(
    chemical_ref: str | int,
    pathway_index: int = 0,
    max_producers_per_chemical: int = 5,
) -> Image:
    """Render one pathway to a PNG for inline display in Claude Desktop.

    Bipartite layered diagram: chemicals (ellipses, green=native, salmon=target)
    and reactions (yellow boxes), laid out top-down with `dot`. Requires the
    `graphviz` Python package + system binary (brew/apt install graphviz).
    """
    hg = state.get_hypergraph()
    chem = _resolve_or_raise(chemical_ref)
    if chem not in hg.chemical_to_shell:
        raise ValueError(
            f"{chem.name!r} (id={chem.id}) is not reachable from the baseline cell."
        )
    cascade = build_cascade(hg, chem, max_producers_per_chemical=max_producers_per_chemical)
    pathways = enumerate_pathways(cascade, hg, max_pathways=pathway_index + 1)
    if pathway_index >= len(pathways):
        raise ValueError(
            f"Only {len(pathways)} pathway(s) found; index {pathway_index} out of range."
        )
    png = render_pathway_png(pathways[pathway_index], hg, state.get_ec_names())
    return Image(data=png, format="png")


@mcp.tool()
def visualize_cascade(
    chemical_ref: str | int,
    max_producers_per_chemical: int = 5,
    max_reactions: int = 80,
) -> Image:
    """Render the full cascade (all producer reactions) as a PNG.

    Same visual scheme as `visualize_pathway`. Raises if the cascade has more
    than *max_reactions* reactions — large cascades produce illegible images;
    lower *max_producers_per_chemical* first.
    """
    hg = state.get_hypergraph()
    chem = _resolve_or_raise(chemical_ref)
    if chem not in hg.chemical_to_shell:
        raise ValueError(
            f"{chem.name!r} (id={chem.id}) is not reachable from the baseline cell."
        )
    cascade = build_cascade(hg, chem, max_producers_per_chemical=max_producers_per_chemical)
    png = render_cascade_png(
        cascade, hg, max_reactions=max_reactions, ec_to_name=state.get_ec_names()
    )
    return Image(data=png, format="png")


@mcp.tool()
def resynthesize_with_fed(
    fed_chemical_refs: list[str | int],
    include_new_reachables: bool = True,
    new_reachables_limit: int = 50,
) -> dict[str, Any]:
    """Re-run synthesize() with additional fed chemicals added to shell 0.

    Results are cached per fed-set so repeated calls with the same inputs
    are fast. Returns summary stats plus (optionally) the list of chemicals
    that became reachable *only* with the fed additions, up to *new_reachables_limit*.
    """
    chems = state.get_chemicals()
    fed: list[Chemical] = []
    unresolved: list[Any] = []
    for ref in fed_chemical_refs:
        hit = lookup.resolve_one(ref, chems)
        if hit is None:
            unresolved.append(ref)
        else:
            fed.append(hit)

    fed_ids = tuple(sorted({c.id for c in fed}))
    baseline = state.get_hypergraph()
    fed_hg = state.get_or_build_fed_hypergraph(fed_ids)

    shells = max(fed_hg.reaction_to_shell.values(), default=0)
    new_reachables: list[dict[str, Any]] = []
    if include_new_reachables:
        diff = [
            c for c in fed_hg.chemical_to_shell
            if c not in baseline.chemical_to_shell
        ]
        diff.sort(key=lambda c: (fed_hg.chemical_to_shell[c], c.id))
        new_reachables = [chemical_to_dto(c, fed_hg) for c in diff[:new_reachables_limit]]

    return {
        "fed_chemicals": [chemical_to_dto(c, fed_hg) for c in fed],
        "unresolved": [str(u) for u in unresolved],
        "shells": shells,
        "reachable_chemicals": len(fed_hg.chemical_to_shell),
        "reachable_reactions": len(fed_hg.reaction_to_shell),
        "baseline_reachable_chemicals": len(baseline.chemical_to_shell),
        "delta_chemicals": len(fed_hg.chemical_to_shell) - len(baseline.chemical_to_shell),
        "delta_reactions": len(fed_hg.reaction_to_shell) - len(baseline.reaction_to_shell),
        "new_reachables": new_reachables,
        "new_reachables_truncated": (
            include_new_reachables
            and (len(fed_hg.chemical_to_shell) - len(baseline.chemical_to_shell))
            > new_reachables_limit
        ),
    }


# -----------------------------------------------------------------------------
# Resources
# -----------------------------------------------------------------------------


@mcp.resource("synthesis://stats")
def stats_resource() -> dict[str, Any]:
    """Summary stats of the loaded hypergraph."""
    return state.get_stats()


@mcp.resource("synthesis://reachables")
def reachables_resource() -> str:
    """TSV dump of all reachables (id, name, inchi, shell)."""
    hg = state.get_hypergraph()
    lines = ["id\tname\tinchi\tshell"]
    for chem, shell in sorted(
        hg.chemical_to_shell.items(), key=lambda kv: (kv[1], kv[0].id)
    ):
        lines.append(f"{chem.id}\t{chem.name}\t{chem.inchi}\t{shell}")
    return "\n".join(lines) + "\n"


@mcp.resource("synthesis://chemical/{chem_id}")
def chemical_resource(chem_id: str) -> dict[str, Any]:
    """Chemical detail + which reactions produce / consume it."""
    hg = state.get_hypergraph()
    chems = state.get_chemicals()
    try:
        cid = int(chem_id)
    except ValueError:
        raise ValueError(f"chemical_id must be integer, got {chem_id!r}")
    if cid not in chems:
        raise ValueError(f"Unknown chemical id {cid}")
    chem = chems[cid]

    produced_by = [r.id for r in hg.reaction_to_shell if chem in r.products]
    consumed_by = [r.id for r in hg.reaction_to_shell if chem in r.substrates]
    return {
        **chemical_to_dto(chem, hg),
        "produced_by_reaction_ids": sorted(produced_by),
        "consumed_by_reaction_ids": sorted(consumed_by),
    }


@mcp.resource("synthesis://reaction/{rxn_id}")
def reaction_resource(rxn_id: str) -> dict[str, Any]:
    """Reaction detail with expanded substrate / product names."""
    hg = state.get_hypergraph()
    try:
        rid = int(rxn_id)
    except ValueError:
        raise ValueError(f"reaction_id must be integer, got {rxn_id!r}")
    for r in state.get_reactions():
        if r.id == rid:
            return {
                **reaction_to_dto(r, hg),
                "substrates": [chemical_to_dto(s, hg) for s in sorted(r.substrates, key=lambda c: c.id)],
                "products": [chemical_to_dto(p, hg) for p in sorted(r.products, key=lambda c: c.id)],
            }
    raise ValueError(f"Unknown reaction id {rid}")


# -----------------------------------------------------------------------------
# Entry point
# -----------------------------------------------------------------------------


def main() -> None:
    """Run the server on stdio (what Claude Code connects to)."""
    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
