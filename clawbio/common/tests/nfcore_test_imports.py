"""Explicit file imports for tests of wrappers with colliding bare module names."""

import hashlib
import importlib.util
import sys
from pathlib import Path
from types import ModuleType


def load_tested_module(skill_dir: Path, name: str) -> ModuleType:
    """Reuse the owner's canonical object so existing monkeypatches still apply.

    The wrapper conftest owns sibling-import isolation. Resolve the tested file
    explicitly so static analysis and foreign cache entries cannot select another
    wrapper's function signature. Include the full path to distinguish checkouts.
    """
    path = (skill_dir / f"{name}.py").resolve()
    key = f"_nfcore_test_{hashlib.sha256(str(path).encode()).hexdigest()}_{name}"
    for candidate in (name, key):
        module = sys.modules.get(candidate)
        if module is not None and Path(getattr(module, "__file__", "") or "").resolve() == path:
            sys.modules[key] = module
            sys.modules[name] = module
            return module
    spec = importlib.util.spec_from_file_location(key, path)
    module = importlib.util.module_from_spec(spec)
    # dataclasses looks up the defining module while executing its decorator.
    sys.modules[key] = module
    try:
        spec.loader.exec_module(module)
    except BaseException:
        sys.modules.pop(key, None)
        raise
    sys.modules[name] = module
    return module
