"""pytest configuration: the skill's flat modules on sys.path, the `live` marker, and the
offline lookup arguments (the same vendored slices `--demo` reads)."""

from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest

SKILL_DIR = Path(__file__).resolve().parent.parent
if str(SKILL_DIR) not in sys.path:
    sys.path.insert(0, str(SKILL_DIR))


def pytest_configure(config: pytest.Config) -> None:
    config.addinivalue_line("markers", "live: hits the live Pan-UKBB / FinnGen / Ensembl / Open Targets services")


def pytest_collection_modifyitems(config: pytest.Config, items: list[pytest.Item]) -> None:
    if os.environ.get("RUN_LIVE_TESTS") == "1":
        return
    skip_live = pytest.mark.skip(reason="live test (set RUN_LIVE_TESTS=1)")
    for item in items:
        if "live" in item.keywords and "snp-phewas" in str(item.fspath):
            item.add_marker(skip_live)


@pytest.fixture
def offline_kwargs():
    from snp_phewas_demo import demo_kwargs

    return demo_kwargs()
