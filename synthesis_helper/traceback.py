"""Traceback: build a Cascade for a target chemical from the HyperGraph."""

from __future__ import annotations

from synthesis_helper.models import Cascade, Chemical, HyperGraph, Reaction


def traceback(hypergraph: HyperGraph, target: Chemical) -> Cascade:
    """Walk backward from target to native metabolites, building a Cascade.

    For each reaction producing the target, recurse on that reaction's
    substrates until shell-0 metabolites are reached.
    """
    if target not in hypergraph.chemical_to_shell:
        raise ValueError(f"Chemical {target.name!r} (id={target.id}) is not reachable.")

    cascade = Cascade(target=target)
    _collect_reactions(hypergraph, target, cascade, visited=set())
    return cascade


def _collect_reactions(
    hg: HyperGraph,
    chemical: Chemical,
    cascade: Cascade,
    visited: set[int],
) -> None:
    """Recursively collect all reactions that contribute to producing *chemical*."""
    if hg.chemical_to_shell.get(chemical) == 0:
        return  # base case: native/universal metabolite

    for rxn in hg.reaction_to_shell:
        if rxn.id in visited:
            continue
        if chemical in rxn.products:
            visited.add(rxn.id)
            cascade.reactions.add(rxn)
            for substrate in rxn.substrates:
                _collect_reactions(hg, substrate, cascade, visited)
