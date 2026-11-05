#!/usr/bin/env python3
"""DAG step: resolve_source.

Thin CLI wrapper around the ported ``pipeline_source.resolve_pipeline_source``
(unchanged real logic from ClawBio 1.x): decide whether nf-core/rnaseq will be
run from a sibling local checkout or fetched remotely at the pinned version,
and enforce the pinned-version contract on a local checkout. No Nextflow is
invoked -- this only inspects the filesystem (and, for a local checkout, runs
``git`` to record branch/commit/dirty state).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

_SCRIPT_DIR = Path(__file__).resolve().parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))

from _isolated_imports import purge_foreign_bare_modules  # noqa: E402

purge_foreign_bare_modules("errors", "schemas", "pipeline_source")

from args_common import none_or  # noqa: E402
from errors import ErrorCode, SkillError  # noqa: E402
from pipeline_source import resolve_pipeline_source  # noqa: E402
from schemas import DEFAULT_PIPELINE_VERSION  # noqa: E402


def _check_pipeline_source_version(*, pipeline_version: str, allow_override: bool, pipeline_source: dict) -> None:
    """Port of nfcore_rnaseq_wrapper.py's ``_check_pipeline_source_version``."""
    if pipeline_source.get("source_kind") != "local_checkout":
        return
    manifest_version = str(pipeline_source.get("manifest_version", "") or "").strip()
    if manifest_version == DEFAULT_PIPELINE_VERSION:
        return
    if not manifest_version:
        print(
            "WARNING: could not determine the local nf-core/rnaseq checkout version "
            f"(expected {DEFAULT_PIPELINE_VERSION}).",
            file=sys.stderr,
        )
        return
    if allow_override:
        print(
            f"WARNING: local nf-core/rnaseq checkout reports version {manifest_version!r}, "
            f"differing from the pinned {DEFAULT_PIPELINE_VERSION} contract.",
            file=sys.stderr,
        )
        return
    raise SkillError(
        stage="validation",
        error_code=ErrorCode.PIPELINE_SOURCE_INVALID,
        message=f"Local nf-core/rnaseq checkout is version {manifest_version!r}, pinned to {DEFAULT_PIPELINE_VERSION}.",
        fix="Check out nf-core/rnaseq 3.26.0, or pass allow_pipeline_version_override=true.",
        details={"manifest_version": manifest_version, "contract_version": DEFAULT_PIPELINE_VERSION},
    )


def main() -> int:
    parser = argparse.ArgumentParser(description="Resolve the nf-core/rnaseq pipeline source.")
    parser.add_argument("--pipeline-version", required=True)
    parser.add_argument("--pipeline-local", required=True, help="Local checkout path, or the sentinel 'none'.")
    parser.add_argument("--allow-pipeline-version-override", required=True)
    parser.add_argument("--out", required=True, help="Path to write pipeline_source.json")
    args = parser.parse_args()

    requested_version = args.pipeline_version.strip() if args.pipeline_version.strip() != "none" else DEFAULT_PIPELINE_VERSION
    local_dir = none_or(args.pipeline_local)
    allow_override = args.allow_pipeline_version_override.strip().lower() in {"1", "true", "yes", "on"}

    try:
        source = resolve_pipeline_source(
            requested_version=requested_version,
            local_pipeline_dir=Path(local_dir).expanduser().resolve() if local_dir else None,
        )
        _check_pipeline_source_version(
            pipeline_version=requested_version, allow_override=allow_override, pipeline_source=source
        )
    except SkillError as exc:
        out_path = Path(args.out)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(json.dumps(exc.to_dict(), indent=2) + "\n", encoding="utf-8")
        print(json.dumps(exc.to_dict(), indent=2), file=sys.stderr)
        return 1

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(source, indent=2) + "\n", encoding="utf-8")
    print(f"[resolve_source] {source['source_kind']} -> {source['resolved_version']}")
    print(json.dumps(source, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
