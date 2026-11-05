from __future__ import annotations

import sys
from pathlib import Path

_SKILL_DIR = Path(__file__).resolve().parent


def purge_foreign_bare_modules(*names: str) -> None:
    """Drop cached top-level modules that do not belong to this skill.

    The skill is executed as a bare script (``python nfcore_rnaseq_wrapper.py``),
    so its modules are imported under top-level names like ``errors`` and
    ``schemas``. A sibling skill that ships an identically-named ``errors.py`` /
    ``schemas.py`` can therefore shadow ours when both run in one interpreter
    (most visibly, the shared test session). Each module calls this guard before
    importing its dependencies so the cache only ever holds *our* copy.

    Centralised here (instead of copy-pasted into every module) so the isolation
    rule has a single definition. Callers ``sys.modules.pop("_isolated_imports")``
    before importing this module, so a foreign copy of this file cannot shadow it
    either.
    """
    for name in names:
        module = sys.modules.get(name)
        module_file = Path(getattr(module, "__file__", "") or "")
        if module is not None and _SKILL_DIR not in module_file.parents and module_file != _SKILL_DIR / f"{name}.py":
            sys.modules.pop(name, None)
