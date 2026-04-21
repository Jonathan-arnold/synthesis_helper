"""Lazy-loaded singleton holding parsed data + the baseline HyperGraph."""

from __future__ import annotations

import os
import threading
import time
from dataclasses import dataclass, field
from pathlib import Path

from synthesis_helper.models import Chemical, HyperGraph, Reaction
from synthesis_helper.parser import (
    parse_chemicals,
    parse_metabolite_list,
    parse_reactions,
)
from synthesis_helper.synthesize import synthesize


DEFAULT_DATA_DIR = Path(__file__).resolve().parents[2] / "data"


@dataclass
class _State:
    data_dir: Path = DEFAULT_DATA_DIR
    chemicals: dict[int, Chemical] | None = None
    reactions: list[Reaction] | None = None
    natives: set[Chemical] | None = None
    universals: set[Chemical] | None = None
    hypergraph: HyperGraph | None = None
    fed_hypergraphs: dict[tuple[int, ...], HyperGraph] = field(default_factory=dict)
    load_time_ms: float = 0.0


_state = _State()
_lock = threading.Lock()


def configure(data_dir: Path | str | None = None) -> None:
    """Override the data directory. Must be called before first access."""
    env = os.environ.get("SYNTHESIS_DATA_DIR")
    if data_dir is not None:
        _state.data_dir = Path(data_dir)
    elif env:
        _state.data_dir = Path(env)


def reset() -> None:
    """Clear all cached state. Intended for tests."""
    with _lock:
        _state.chemicals = None
        _state.reactions = None
        _state.natives = None
        _state.universals = None
        _state.hypergraph = None
        _state.fed_hypergraphs.clear()
        _state.load_time_ms = 0.0


def _bootstrap() -> None:
    """Parse files, seed shell 0, run BFS. Idempotent under lock."""
    if _state.hypergraph is not None:
        return

    configure()
    data = _state.data_dir
    t0 = time.perf_counter()

    chemicals = parse_chemicals(data / "good_chems.txt")
    reactions = parse_reactions(data / "good_reactions.txt", chemicals)
    natives = parse_metabolite_list(data / "minimal_metabolites.txt", chemicals)
    universals = parse_metabolite_list(
        data / "ubiquitous_metabolites.txt", chemicals
    )

    hg = synthesize(reactions, natives, universals, verbose=False)

    _state.chemicals = chemicals
    _state.reactions = reactions
    _state.natives = natives
    _state.universals = universals
    _state.hypergraph = hg
    _state.load_time_ms = (time.perf_counter() - t0) * 1000.0


def get_hypergraph() -> HyperGraph:
    with _lock:
        if _state.hypergraph is None:
            _bootstrap()
        assert _state.hypergraph is not None
        return _state.hypergraph


def get_chemicals() -> dict[int, Chemical]:
    with _lock:
        if _state.chemicals is None:
            _bootstrap()
        assert _state.chemicals is not None
        return _state.chemicals


def get_reactions() -> list[Reaction]:
    with _lock:
        if _state.reactions is None:
            _bootstrap()
        assert _state.reactions is not None
        return _state.reactions


def get_natives() -> set[Chemical]:
    with _lock:
        if _state.natives is None:
            _bootstrap()
        assert _state.natives is not None
        return _state.natives


def get_universals() -> set[Chemical]:
    with _lock:
        if _state.universals is None:
            _bootstrap()
        assert _state.universals is not None
        return _state.universals


def get_stats() -> dict[str, float | int]:
    hg = get_hypergraph()
    shells = max(hg.reaction_to_shell.values(), default=0)
    return {
        "shells": shells,
        "total_chemicals": len(_state.chemicals or {}),
        "reachable_chemicals": len(hg.chemical_to_shell),
        "total_reactions": len(_state.reactions or []),
        "reachable_reactions": len(hg.reaction_to_shell),
        "load_time_ms": round(_state.load_time_ms, 1),
        "data_dir": str(_state.data_dir),
    }


def get_or_build_fed_hypergraph(fed_ids: tuple[int, ...]) -> HyperGraph:
    """Return a cached (or freshly computed) HyperGraph with fed chemicals."""
    with _lock:
        if fed_ids in _state.fed_hypergraphs:
            return _state.fed_hypergraphs[fed_ids]

        if _state.hypergraph is None:
            _bootstrap()
        assert _state.chemicals is not None
        assert _state.reactions is not None
        assert _state.natives is not None
        assert _state.universals is not None

        extra = {_state.chemicals[cid] for cid in fed_ids if cid in _state.chemicals}
        hg = synthesize(
            _state.reactions,
            native_metabolites=_state.natives | extra,
            universal_metabolites=_state.universals,
            verbose=False,
        )
        _state.fed_hypergraphs[fed_ids] = hg
        return hg
