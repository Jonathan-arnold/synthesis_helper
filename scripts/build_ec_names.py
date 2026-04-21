"""Rebuild data/ec_names.tsv from ExPASy's enzyme.dat.

Usage:
    curl -sS https://ftp.expasy.org/databases/enzyme/enzyme.dat -o /tmp/enzyme.dat
    python scripts/build_ec_names.py /tmp/enzyme.dat

Parses ID + DE lines, resolves "Transferred entry" chains to the new EC's
name, skips deleted entries. Upstream license: CC BY 4.0
(https://enzyme.expasy.org/, SIB Swiss Institute of Bioinformatics).
"""

from __future__ import annotations

import re
import sys
from pathlib import Path


def parse(enzyme_dat: Path) -> dict[str, str]:
    src = enzyme_dat.read_text()
    records = src.split("\n//\n")

    raw_name: dict[str, str] = {}
    transfer_targets: dict[str, list[str]] = {}
    deleted: set[str] = set()

    transfer_re = re.compile(r"transferred entry:\s*(.+)", re.I)
    deleted_re = re.compile(r"^deleted entry", re.I)

    for rec in records:
        ec: str | None = None
        de_parts: list[str] = []
        in_de = False
        for line in rec.splitlines():
            if line.startswith("ID   "):
                ec = line[5:].strip()
            elif line.startswith("DE   "):
                de_parts.append(line[5:].strip())
                in_de = True
            elif in_de and line.startswith("     "):
                de_parts.append(line.strip())
            elif line.startswith(("AN ", "CA ", "CF ", "CC ", "PR ", "DR ")):
                in_de = False
        if not ec or not de_parts:
            continue
        de = " ".join(de_parts).rstrip(".").strip()
        m = transfer_re.match(de)
        if m:
            targets = re.findall(r"\d+\.\d+\.\d+\.\d+", m.group(1))
            if targets:
                transfer_targets[ec] = targets
                continue
        if deleted_re.match(de):
            deleted.add(ec)
            continue
        raw_name[ec] = de

    def resolve(ec: str, seen: set[str]) -> str | None:
        if ec in raw_name:
            return raw_name[ec]
        if ec in deleted or ec in seen:
            return None
        seen.add(ec)
        for tgt in transfer_targets.get(ec, []):
            hit = resolve(tgt, seen)
            if hit:
                return hit
        return None

    resolved: dict[str, str] = dict(raw_name)
    for old_ec in transfer_targets:
        name = resolve(old_ec, set())
        if name:
            resolved[old_ec] = name
    return resolved


def _sort_key(ec: str) -> list[int]:
    return [int(p) if p.isdigit() else 0 for p in ec.split(".")]


def main() -> None:
    if len(sys.argv) != 2:
        print("usage: build_ec_names.py <path/to/enzyme.dat>", file=sys.stderr)
        sys.exit(2)
    names = parse(Path(sys.argv[1]))
    out = Path(__file__).resolve().parents[1] / "data" / "ec_names.tsv"
    with out.open("w") as fh:
        fh.write("ecnum\tname\n")
        for ec in sorted(names, key=_sort_key):
            fh.write(f"{ec}\t{names[ec].replace(chr(9), ' ')}\n")
    print(f"wrote {len(names)} entries -> {out}", file=sys.stderr)


if __name__ == "__main__":
    main()
