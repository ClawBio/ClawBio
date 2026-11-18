"""The Privacy section must match what the outputs actually contain.

`nutrigx_report.md` prints a per-SNP genotype column and `result.json` stores a
genotype call for every SNP, so the section cannot promise genotype-free output
(issue #793).
"""

from pathlib import Path

SKILL_DIR = Path(__file__).parent.parent


def privacy_section() -> str:
    text = (SKILL_DIR / "SKILL.md").read_text(encoding="utf-8")
    return text.split("## Privacy", 1)[1].split("\n## ", 1)[0]


def test_privacy_section_is_not_empty():
    # Guards the split above: a renamed heading must fail loudly, not vacuously.
    assert "genotype" in privacy_section().lower() or "local" in privacy_section().lower()


def test_privacy_section_does_not_promise_genotype_free_outputs():
    privacy = privacy_section().lower()
    assert "no raw genotype data" not in privacy
    assert "result.json" in privacy
