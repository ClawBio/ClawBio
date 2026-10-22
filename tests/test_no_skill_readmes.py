"""A skill folder holds SKILL.md, not a README.

Agents load SKILL.md and references/, so a README beside it drifts unnoticed;
Anthropic's skills guide also says not to include one. READMEs describing
fixtures or example data in subfolders are fine; only the skill root is checked.
"""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# Skills that still carry a root README. Fix one per PR: fold anything unique
# into SKILL.md or references/, delete the README, and remove the name here.
LEGACY_READMES = {
    "nfcore-rnaseq-wrapper",
    "nfcore-rnastructurome-wrapper",
    "nfcore-sarek-wrapper",
    "nfcore-scrnaseq-wrapper",
    "seq-wrangler",
}


def _skills_with_root_readme() -> set[str]:
    return {
        p.parent.name
        for p in ROOT.glob("skills/*/*")
        if p.is_file() and p.name.lower().startswith("readme")
    }


def test_no_new_skill_readmes():
    new = _skills_with_root_readme() - LEGACY_READMES
    assert not new, (
        f"README in skill folder: {sorted(new)}. Put agent-facing detail in "
        "SKILL.md or references/, human-facing material in the repo README or docs/."
    )


def test_legacy_list_only_shrinks():
    fixed = LEGACY_READMES - _skills_with_root_readme()
    assert not fixed, f"README removed for {sorted(fixed)}; drop them from LEGACY_READMES."
