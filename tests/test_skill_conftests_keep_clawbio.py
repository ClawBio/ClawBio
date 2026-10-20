"""A skill's conftest must not unload or replace clawbio modules already imported.

pytest loads skills/<x>/tests/conftest.py during collection, after tests/ has
imported clawbio. A conftest that wiped sys.modules (wgs-prs did) left modules
imported earlier holding a second copy of clawbio.cli, so
test_mcp_server's identity check failed whenever the two shared a process.
"""

import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
CONFTESTS = sorted(ROOT.glob("skills/*/tests/conftest.py"))

CHECK = """
import runpy, sys
import clawbio.cli
before = {n: id(m) for n, m in sys.modules.items() if n == "clawbio" or n.startswith("clawbio.")}
runpy.run_path(sys.argv[1])
after = {n: id(m) for n, m in sys.modules.items() if n == "clawbio" or n.startswith("clawbio.")}
changed = sorted(n for n in before if after.get(n) != before[n])
if changed:
    sys.exit("removed or replaced: " + ", ".join(changed))
"""


@pytest.mark.parametrize("conftest", CONFTESTS, ids=lambda p: p.parts[-3])
def test_skill_conftest_leaves_loaded_clawbio_modules_alone(conftest):
    result = subprocess.run(
        [sys.executable, "-c", CHECK, str(conftest)],
        cwd=ROOT, capture_output=True, text=True, timeout=120,
    )
    assert result.returncode == 0, result.stderr[-2000:]
