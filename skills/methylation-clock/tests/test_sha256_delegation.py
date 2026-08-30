"""Checksum helper must delegate to the shared clawbio.common layer.

Kept separate from test_methylation_clock.py so it runs even where pyaging
(a heavy optional dependency) cannot be imported: pyaging is stubbed out
before importing the module, since checksum behaviour does not depend on it.
"""

import sys
from pathlib import Path
from unittest.mock import MagicMock

SKILL_DIR = Path(__file__).resolve().parent.parent
PROJECT_ROOT = SKILL_DIR.parents[1]
for p in (SKILL_DIR, PROJECT_ROOT):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

sys.modules.setdefault("pyaging", MagicMock())

import methylation_clock  # noqa: E402
from clawbio.common.checksums import sha256_file  # noqa: E402


def test_sha256_delegates_to_common_layer():
    assert methylation_clock._sha256 is sha256_file
