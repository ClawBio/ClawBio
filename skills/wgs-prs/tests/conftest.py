"""
conftest.py: fix clawbio package shadowing for wgs-prs tests.

The repo root contains clawbio.py (the CLI) which shadows the clawbio/
package directory. This conftest removes the stale module entry and re-imports
the package from the correct path before any tests run.
"""
import importlib
import sys
from pathlib import Path

# Project root: skills/wgs-prs/tests -> skills/wgs-prs -> skills -> project root
_PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent.parent

def _fix_clawbio_package() -> None:
    """Ensure 'clawbio' resolves to the clawbio/ package, not clawbio.py."""
    clawbio_pkg = _PROJECT_ROOT / "clawbio"
    if not clawbio_pkg.is_dir():
        return

    # Remove any stale clawbio module loaded from clawbio.py
    for key in list(sys.modules.keys()):
        if key == "clawbio" or key.startswith("clawbio."):
            del sys.modules[key]

    # Insert the project root at position 0 and load the package explicitly
    if str(_PROJECT_ROOT) not in sys.path:
        sys.path.insert(0, str(_PROJECT_ROOT))

    # Force Python to find clawbio/ (directory) not clawbio.py (file)
    # by loading it via its __init__.py path
    spec = importlib.util.spec_from_file_location(
        "clawbio",
        str(clawbio_pkg / "__init__.py"),
        submodule_search_locations=[str(clawbio_pkg)],
    )
    if spec:
        mod = importlib.util.module_from_spec(spec)
        sys.modules["clawbio"] = mod
        spec.loader.exec_module(mod)

_fix_clawbio_package()
