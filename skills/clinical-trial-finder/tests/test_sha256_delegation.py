"""Checksum helper must delegate to the shared clawbio.common layer."""

import sys
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parent.parent
PROJECT_ROOT = SKILL_DIR.parents[1]
for p in (SKILL_DIR, PROJECT_ROOT):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

import writers  # noqa: E402
from clawbio.common.checksums import sha256_file  # noqa: E402


def test_sha256_delegates_to_common_layer():
    assert writers._sha256 is sha256_file
