#!/usr/bin/env python3
"""gi-annotation: DNA annotation (gene/transcript) via Genomic Intelligence /v1/tasks/annotation/predict.

The GI API accepts either delivery mode on every task; this skill always
submits via ``Prefer: respond-async`` and polls until terminal, because the
annotation pipeline is long-running (~20 s for ~20 kbp).
"""

from __future__ import annotations

import sys
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SKILL_DIR.parent.parent))

from clawbio.gi.gi_runner import run_skill


def main() -> int:
    return run_skill(
        task="annotation",
        demo_path=SKILL_DIR / "example_data" / "annotation_tp53.fa",
        async_mode=True,
    )


if __name__ == "__main__":
    sys.exit(main())
