"""BFS waveform expansion — the core synthesis algorithm."""

from __future__ import annotations

from synthesis_helper.models import Chemical, HyperGraph, Reaction


def synthesize(
    all_reactions: list[Reaction],
    native_metabolites: set[Chemical],
    universal_metabolites: set[Chemical],
    verbose: bool = True,
) -> HyperGraph:
    """Run BFS waveform expansion and return the populated HyperGraph.

    Algorithm:
    1. Initialize shell 0 with native ∪ universal metabolites.
    2. Each round, enable reactions whose substrates are all reachable.
    3. Add new products at the current shell number.
    4. Repeat until no new reactions or chemicals are added (fixed point).
    """
    hg = HyperGraph()

    # Shell 0: starting pool
    for chem in native_metabolites | universal_metabolites:
        hg.chemical_to_shell[chem] = 0

    enabled_reactions: set[int] = set()
    shell = 1

    while True:
        new_reactions: list[Reaction] = []
        new_chemicals: dict[Chemical, int] = {}

        for rxn in all_reactions:
            if rxn.id in enabled_reactions:
                continue
            if all(s in hg.chemical_to_shell for s in rxn.substrates):
                new_reactions.append(rxn)
                enabled_reactions.add(rxn.id)
                hg.reaction_to_shell[rxn] = shell
                for product in rxn.products:
                    if product not in hg.chemical_to_shell:
                        new_chemicals[product] = shell

        if not new_reactions:
            break

        for chem, s in new_chemicals.items():
            hg.chemical_to_shell[chem] = s

        if verbose:
            print(
                f"Shell {shell}: +{len(new_reactions)} reactions, "
                f"+{len(new_chemicals)} chemicals "
                f"(total: {len(hg.reaction_to_shell)} rxns, "
                f"{len(hg.chemical_to_shell)} chems)"
            )

        shell += 1

    if verbose:
        print(
            f"Expansion complete: {len(hg.chemical_to_shell)} reachable chemicals, "
            f"{len(hg.reaction_to_shell)} enabled reactions, "
            f"{shell - 1} shells."
        )

    return hg
