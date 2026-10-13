"""Isolation for legacy bare imports, with an explicit owning skill."""

import sys
from pathlib import Path


def purge_foreign_bare_modules(skill_dir: Path, *names: str) -> None:
    for name in names:
        module = sys.modules.get(name)
        module_file = Path(getattr(module, "__file__", "") or "")
        if (
            module is not None
            and skill_dir not in module_file.parents
            and module_file != skill_dir / f"{name}.py"
        ):
            sys.modules.pop(name, None)
