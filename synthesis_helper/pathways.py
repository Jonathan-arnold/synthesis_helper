"""Enumerate individual Pathways from a Cascade."""

from __future__ import annotations

import itertools
from itertools import islice

from synthesis_helper.models import Cascade, Chemical, HyperGraph, Pathway, Reaction


def enumerate_pathways(
    cascade: Cascade,
    hypergraph: HyperGraph,
    max_pathways: int = 1000,
) -> list[Pathway]:
    """Enumerate pathways implied by a cascade.

    A Pathway is a set of reactions whose combined products include the target
    and whose substrates are all either produced by some reaction in the set
    or native (shell 0). Reactions are returned in topological order — every
    reaction comes after the reactions producing its non-native substrates.

    For reactions with multiple non-native substrates (e.g. ligase reactions,
    aldol condensations), pathways branch across all substrates: the final
    pathway includes a sub-pathway to produce each one. This is the cartesian
    product of per-substrate sub-pathway sets, capped at max_pathways.

    Args:
        cascade: the cascade produced by traceback().
        hypergraph: the HyperGraph (used for shell lookups).
        max_pathways: global cap on returned pathways. Also bounds per-
            substrate enumeration during recursion, to keep combinatorial
            blowup from exhausting memory before the outer cap fires.
    """
    producers: dict[Chemical, list[Reaction]] = {}
    for rxn in cascade.reactions:
        for product in rxn.products:
            producers.setdefault(product, []).append(rxn)

    in_progress: set[int] = set()

    def paths_for(chem: Chemical) -> "itertools.chain[list[Reaction]]":
        """Yield lists of reactions that produce chem, in topological order.

        Empty list => chem is shell 0 (no reactions needed). Cycles are
        broken via in_progress: if chem is already being produced upstream
        in the current branch, we yield nothing.
        """
        if hypergraph.chemical_to_shell.get(chem) == 0:
            yield []
            return
        if chem.id in in_progress:
            return

        in_progress.add(chem.id)
        try:
            for rxn in producers.get(chem, []):
                non_native = sorted(
                    (s for s in rxn.substrates
                     if hypergraph.chemical_to_shell.get(s) != 0),
                    key=lambda s: s.id,
                )

                if not non_native:
                    yield [rxn]
                    continue

                sub_path_lists: list[list[list[Reaction]]] = []
                for sub in non_native:
                    sp = list(islice(paths_for(sub), max_pathways))
                    if not sp:
                        break
                    sub_path_lists.append(sp)
                if len(sub_path_lists) != len(non_native):
                    continue

                for combo in itertools.product(*sub_path_lists):
                    merged: list[Reaction] = []
                    seen: set[int] = set()
                    for sub_rxns in combo:
                        for r in sub_rxns:
                            if r.id not in seen:
                                merged.append(r)
                                seen.add(r.id)
                    if rxn.id not in seen:
                        merged.append(rxn)
                    yield merged
        finally:
            in_progress.discard(chem.id)

    results: list[Pathway] = []
    for rxns in islice(paths_for(cascade.target), max_pathways):
        metabolites: set[Chemical] = set()
        for rxn in rxns:
            metabolites.update(rxn.substrates)
            metabolites.update(rxn.products)
        results.append(
            Pathway(
                target=cascade.target,
                reactions=rxns,
                metabolites=metabolites,
            )
        )

    return results
