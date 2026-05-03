"""Resolve a user-supplied chemical reference (id, name, or InChI) to Chemical."""

from __future__ import annotations

import re

from synthesis_helper.models import Chemical


_PROTON = re.compile(r"/p[+-]\d+")


def _strip_proton(inchi: str) -> str:
    return _PROTON.sub("", inchi)


def resolve_one(
    ref: str | int,
    chemicals: dict[int, Chemical],
) -> Chemical | None:
    """Best-effort resolve. Returns None when no candidate found."""
    matches = resolve(ref, chemicals, limit=1)
    return matches[0] if matches else None


def resolve(
    ref: str | int,
    chemicals: dict[int, Chemical],
    limit: int = 10,
) -> list[Chemical]:
    """Resolve a reference to up to *limit* Chemicals.

    Resolution order:
      1. Integer id (or str digits) → exact match.
      2. Exact case-insensitive name match.
      3. InChI strict match (ionization-normalized).
      4. Substring name match (case-insensitive).
    """
    # 1. id
    if isinstance(ref, int):
        if ref in chemicals:
            return [chemicals[ref]]
        return []

    s = ref.strip()
    if not s:
        return []

    if s.isdigit():
        cid = int(s)
        if cid in chemicals:
            return [chemicals[cid]]

    # 2. exact name (case-insensitive)
    s_lower = s.lower()
    exact = [c for c in chemicals.values() if c.name.lower() == s_lower]
    if exact:
        return exact[:limit]

    # 3. InChI strict (if the ref looks like one)
    if s.startswith("InChI="):
        key = _strip_proton(s)
        inchi_hits = [c for c in chemicals.values() if _strip_proton(c.inchi) == key]
        if inchi_hits:
            return inchi_hits[:limit]

    # 4. substring name
    subs = [c for c in chemicals.values() if s_lower in c.name.lower()]
    # Shorter names first (more relevant when query is specific), then by id
    subs.sort(key=lambda c: (len(c.name), c.id))
    return subs[:limit]
