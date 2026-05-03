"""Shared pytest fixtures for the synthesis_helper test suite."""

from __future__ import annotations

import pytest

from synthesis_helper.models import HyperGraph
from synthesis_helper.synthesize import synthesize
from tests.fixtures.small_cell import SmallCell, build_small_cell


@pytest.fixture
def small_cell() -> SmallCell:
    return build_small_cell()


@pytest.fixture
def small_cell_hg(small_cell: SmallCell) -> HyperGraph:
    return synthesize(
        small_cell.all_reactions,
        native_metabolites=small_cell.natives,
        universal_metabolites=small_cell.universals,
        verbose=False,
    )
