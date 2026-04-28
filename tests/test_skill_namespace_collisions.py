"""Regression tests for cross-skill import namespace collisions."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent


def _clear_modules(*names: str) -> None:
    """Remove cached modules so each scenario starts fresh."""
    for name in names:
        sys.modules.pop(name, None)


def _load_module(module_name: str, path: Path) -> None:
    spec = importlib.util.spec_from_file_location(module_name, path)
    if spec is None or spec.loader is None:
        raise AssertionError(f"Could not build import spec for {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)


def test_struct_predictor_and_fine_mapping_import_in_same_process():
    _clear_modules("core", "struct_predictor", "fine_mapping")

    _load_module(
        "struct_predictor",
        ROOT / "skills" / "struct-predictor" / "struct_predictor.py",
    )
    _load_module(
        "fine_mapping",
        ROOT / "skills" / "fine-mapping" / "fine_mapping.py",
    )


def test_data_extractor_and_gwas_lookup_import_in_same_process():
    _clear_modules("core", "data_extractor", "gwas_lookup")

    _load_module(
        "data_extractor",
        ROOT / "skills" / "data-extractor" / "data_extractor.py",
    )
    _load_module(
        "gwas_lookup",
        ROOT / "skills" / "gwas-lookup" / "gwas_lookup.py",
    )
