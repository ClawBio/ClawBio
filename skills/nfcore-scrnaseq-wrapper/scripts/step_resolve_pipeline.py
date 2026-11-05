#!/usr/bin/env python3
"""clawbio step: resolve_pipeline_source.

Locates the nf-core/scrnaseq pipeline source that a real Nextflow run would
use: a verified local sibling `../scrnaseq` checkout at the requested
version, or (the expected outcome in this sandbox, which has no such
checkout) the pinned remote `nf-core/scrnaseq` tag. Ported unchanged from
ClawBio 1.x's ``pipeline_source.py`` — see that file for the full resolution
logic. This step does not download or execute anything; it only decides which
source string a later `nextflow run` invocation should reference.

This is a genuine hard-fail boundary (unlike the `preflight` step): an
unverifiable or dirty local checkout is a real configuration error the DAG
should stop on, matching ClawBio 1.x's behavior.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

_SCRIPTS_DIR = Path(__file__).resolve().parent
if str(_SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS_DIR))

from errors import SkillError
from pipeline_source import resolve_pipeline_source
from textio import write_text_lf


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Resolve the nf-core/scrnaseq pipeline source (local sibling checkout or pinned remote tag)."
    )
    parser.add_argument("--pipeline-version", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    try:
        source = resolve_pipeline_source(requested_version=args.pipeline_version)
    except SkillError as exc:
        print(json.dumps(exc.to_dict(), indent=2), file=sys.stderr)
        return 1

    write_text_lf(args.out, json.dumps(source, indent=2, default=str))
    print(f"[resolve_pipeline_source] {source['source_kind']} @ {source['resolved_version']}")
    print(json.dumps(source, indent=2, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
