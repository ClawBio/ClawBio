"""A skill's version is written once, in SKILL.md `metadata.version`.

Skill code reads it with `clawbio.common.report.skill_version(__file__)`. On
2026-10-22, 43 skills also typed the version into their Python, and 8 of those
disagreed with SKILL.md (gwas-lookup's SKILL.md said 0.1.0 while its result.json
said 0.2.0), so the catalog and the reports named different versions of the
same code.

How to fix a failure: replace the literal with `skill_version(__file__)`. When
you migrate a skill listed in KNOWN_GAPS, delete it from the list.
"""

from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SKILLS = ROOT / "skills"

# `version="1.2.3"`, `skill_version="1.2.3"`, `"version": "1.2.3"`, `SKILL_VERSION = "1.2.3"`
LITERAL = re.compile(
    r"""\b(?:skill_)?version\s*=\s*["']\d+\.\d+\.\d+"""
    r"""|["']version["']\s*:\s*["']\d+\.\d+\.\d+"""
    r"""|\bSKILL_VERSION\s*=\s*["']\d+\.\d+\.\d+"""
)

# Skills not yet migrated. Some hits may be upstream or tool versions rather than
# the skill's own; classify each one when migrating it.
KNOWN_GAPS = {
    "affinity-proteomics", "analyze-fasta", "bigquery-public", "bio-orchestrator",
    "claw-ancestry-pca", "clinical-variant-reporter", "clinpgx", "data-extractor",
    "eqtl-catalogue-region-fetch", "equity-scorer", "fine-mapping", "flow-bio",
    "galaxy-bridge", "genome-compare", "gwas-catalog-region-fetch", "gwas-lookup", "gwas-pipeline",
    "gwas-prs", "gwas-sumstats-harmonize", "hla-typing", "illumina-bridge",
    "just-prs-mcp", "labstep", "ld-1000g-region-compute", "locuscompare-region-render",
    "mendelian-randomisation", "nfcore-rnaseq-wrapper", "nfcore-sarek-wrapper",
    "nfcore-scrnaseq-wrapper", "nutrigx", "omics-target-evidence-mapper",
    "pathway-enricher", "pharmgx-reporter", "profile-report", "proteomics-de",
    "pubmed-summariser", "rnaseq-de", "scrna-embedding", "scrna-orchestrator",
    "spatial-transcriptomics", "turingdb-graph", "ukb-ppp-region-fetch",
    "variant-annotation",
}


def skills_with_version_literals() -> dict[str, list[str]]:
    found: dict[str, list[str]] = {}
    for path in sorted(SKILLS.rglob("*.py")):
        rel = path.relative_to(SKILLS)
        if "tests" in rel.parts:
            continue
        for number, line in enumerate(path.read_text(errors="ignore").splitlines(), 1):
            if not line.lstrip().startswith("#") and LITERAL.search(line):
                found.setdefault(rel.parts[0], []).append(f"skills/{rel}:{number}")
    return found


def test_no_new_skill_types_its_version_in_code():
    new = {s: hits for s, hits in skills_with_version_literals().items() if s not in KNOWN_GAPS}
    assert not new, f"read the version with skill_version(__file__) instead: {new}"


def test_known_gaps_only_lists_skills_that_still_have_literals():
    fixed = KNOWN_GAPS - skills_with_version_literals().keys()
    assert not fixed, f"migrated; delete from KNOWN_GAPS: {sorted(fixed)}"
