"""I/O utilities: dump reachables, load data files."""

from __future__ import annotations

from pathlib import Path

from synthesis_helper.models import HyperGraph


def dump_reachables(hypergraph: HyperGraph, filepath: str | Path) -> None:
    """Write the reachables list to a tab-separated file.

    Output format: id\tname\tinchi\tshell
    """
    sorted_chems = sorted(
        hypergraph.chemical_to_shell.items(), key=lambda x: (x[1], x[0].id)
    )
    with open(filepath, "w") as f:
        f.write("id\tname\tinchi\tshell\n")
        for chem, shell in sorted_chems:
            f.write(f"{chem.id}\t{chem.name}\t{chem.inchi}\t{shell}\n")
