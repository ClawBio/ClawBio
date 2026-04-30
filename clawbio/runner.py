"""Package-level bridge to the root clawbio.py runner module."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import ModuleType


_ROOT_RUNNER_MODULE_NAME = "clawbio._root_runner"
_ROOT_RUNNER_PATH = Path(__file__).resolve().parent.parent / "clawbio.py"


def _load_root_runner() -> ModuleType:
    """Load the root runner module once and return the cached module."""
    cached = sys.modules.get(_ROOT_RUNNER_MODULE_NAME)
    if cached is not None:
        return cached

    spec = importlib.util.spec_from_file_location(_ROOT_RUNNER_MODULE_NAME, _ROOT_RUNNER_PATH)
    if spec is None or spec.loader is None:
        raise ImportError(f"Unable to load ClawBio runner from {_ROOT_RUNNER_PATH}")

    module = importlib.util.module_from_spec(spec)
    sys.modules[_ROOT_RUNNER_MODULE_NAME] = module
    spec.loader.exec_module(module)
    return module


_ROOT_RUNNER = _load_root_runner()

run_skill = _ROOT_RUNNER.run_skill
list_skills = _ROOT_RUNNER.list_skills
upload_profile = _ROOT_RUNNER.upload_profile

__all__ = ["run_skill", "list_skills", "upload_profile"]
