"""Every option the DnaSP skill declares must survive `clawbio run dnasp`.

The dispatcher drops any extra flag missing from the skill's allowlist, and it does
so silently: `clawbio run dnasp ... --n-sim 10000` once ran with no simulation at
all, because --n-sim, --sim-given and --sim-seed were declared in SKILL.md and
accepted by dnasp.py but absent from the allowlist. This test ties the two together,
so a flag added to the skill's metadata cannot be forgotten here.
"""
from __future__ import annotations

import re
from pathlib import Path

import yaml

from clawbio.cli import SKILLS

REPO = Path(__file__).resolve().parents[1]
RUNNER_OWNED = {"--input", "--output", "--demo"}   # supplied by the runner itself


def _declared_flags() -> set[str]:
    text = (REPO / "skills" / "dnasp" / "SKILL.md").read_text(encoding="utf-8")
    front = re.match(r"^---\s*\n(.*?)\n---", text, re.DOTALL).group(1)
    inputs = yaml.safe_load(front)["metadata"]["inputs"]
    return {i["cli_flag"] for i in inputs if i.get("cli_flag")}


def test_every_declared_dnasp_flag_is_forwarded_by_the_dispatcher():
    entry = SKILLS["dnasp"]
    forwarded = set(entry.get("allowed_extra_flags", set())) | set(entry.get("extra_path_flags", set()))
    missing = sorted(_declared_flags() - RUNNER_OWNED - forwarded)
    assert not missing, f"declared in SKILL.md but dropped by the dispatcher: {missing}"


def test_the_coalescent_options_specifically_are_forwarded():
    allowed = SKILLS["dnasp"]["allowed_extra_flags"]
    assert {"--n-sim", "--sim-given", "--sim-seed"} <= allowed
