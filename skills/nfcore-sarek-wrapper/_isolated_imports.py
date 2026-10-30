"""Keep this skill's bare modules isolated in a shared interpreter."""
from pathlib import Path
from clawbio.common.nfcore.imports import purge_foreign_bare_modules as _purge

_SKILL_DIR = Path(__file__).resolve().parent


def purge_foreign_bare_modules(*names: str) -> None:
    _purge(_SKILL_DIR, *names)
