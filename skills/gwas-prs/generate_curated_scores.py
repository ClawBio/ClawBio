#!/usr/bin/env python3
"""Regenerate curated_scores.json from CURATED_SCORES in gwas_prs.py.

`curated_scores.json` is read by no code, which is exactly why it drifted from
the Python dict and ended up citing a different paper for the BMI panel (issue
#356). SKILL.md called it "generated" before this script existed, which was a
claim rather than a mechanism. Now it is a mechanism.

    python3 skills/gwas-prs/generate_curated_scores.py [--check]

`--check` exits 1 if the committed file differs from what would be generated,
so CI or a pre-commit hook can enforce it. `tests/test_score_provenance.py`
pins the two against each other field by field regardless.
"""

from __future__ import annotations

import argparse
import json
import sys
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parent
TARGET = SKILL_DIR / "curated_scores.json"

# Mirrored into the JSON. `reference_distribution` last so the file reads in
# the same order as the Python dict.
FIELDS = (
    "name",
    "curated_panel_id",
    "legacy_pgs_id",
    "pgs_catalog_id",
    "trait_id",
    "trait",
    "variants_count",
    "publication",
    "pmid",
    "reference_distribution",
)


def load_curated() -> dict:
    spec = spec_from_file_location("gwas_prs", SKILL_DIR / "gwas_prs.py")
    mod = module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.CURATED_SCORES


def render(curated: dict) -> str:
    out = {
        pgs_id: {field: meta[field] for field in FIELDS}
        for pgs_id, meta in curated.items()
    }
    return json.dumps(out, indent=2) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check",
        action="store_true",
        help="exit 1 if the committed file is stale, instead of rewriting it",
    )
    args = parser.parse_args()

    rendered = render(load_curated())

    if args.check:
        current = TARGET.read_text() if TARGET.exists() else ""
        if current != rendered:
            print(
                f"{TARGET.name} is stale. Run: python3 "
                f"{Path(__file__).name}",
                file=sys.stderr,
            )
            return 1
        print(f"{TARGET.name} is up to date.")
        return 0

    TARGET.write_text(rendered)
    print(f"Wrote {TARGET.relative_to(SKILL_DIR.parent.parent)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
