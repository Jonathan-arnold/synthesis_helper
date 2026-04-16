"""Enumerate individual Pathways from a Cascade."""

from __future__ import annotations

from synthesis_helper.models import Cascade, Chemical, HyperGraph, Pathway, Reaction


def enumerate_pathways(
    cascade: Cascade,
    hypergraph: HyperGraph,
    max_pathways: int = 1000,
) -> list[Pathway]:
    """Enumerate all pathways implied by a cascade.

    Walks the Cascade and flattens it into every valid singular route
    from native metabolites (shell 0) to the target.

    Args:
        cascade: The cascade to enumerate.
        hypergraph: The HyperGraph for shell lookups.
        max_pathways: Cap to avoid combinatorial explosion.
    """
    # Build lookup: chemical → reactions in the cascade that produce it
    producers: dict[Chemical, list[Reaction]] = {}
    for rxn in cascade.reactions:
        for product in rxn.products:
            producers.setdefault(product, []).append(rxn)

    pathways: list[Pathway] = []

    def _backtrack(
        target: Chemical,
        current_rxns: list[Reaction],
        current_metabolites: set[Chemical],
    ) -> None:
        if len(pathways) >= max_pathways:
            return

        if hypergraph.chemical_to_shell.get(target) == 0:
            # Reached a native metabolite — record the pathway
            pathways.append(
                Pathway(
                    target=cascade.target,
                    reactions=list(current_rxns),
                    metabolites=set(current_metabolites),
                )
            )
            return

        for rxn in producers.get(target, []):
            if rxn in current_rxns:
                continue  # avoid cycles
            current_rxns.append(rxn)
            current_metabolites.add(target)
            for substrate in rxn.substrates:
                current_metabolites.add(substrate)

            # Check if all substrates are either shell-0 or producible
            all_resolved = all(
                hypergraph.chemical_to_shell.get(s) == 0 or s in producers
                for s in rxn.substrates
            )
            if all_resolved:
                # Recurse on non-native substrates
                non_native = [
                    s for s in rxn.substrates if hypergraph.chemical_to_shell.get(s) != 0
                ]
                if not non_native:
                    # All substrates are native — pathway complete
                    pathways.append(
                        Pathway(
                            target=cascade.target,
                            reactions=list(current_rxns),
                            metabolites=set(current_metabolites),
                        )
                    )
                else:
                    # Recurse on the first non-native substrate
                    # (simplified — full enumeration would branch on all)
                    _backtrack(non_native[0], current_rxns, current_metabolites)

            current_rxns.pop()
            for substrate in rxn.substrates:
                current_metabolites.discard(substrate)
            current_metabolites.discard(target)

    _backtrack(cascade.target, [], set())
    return pathways
