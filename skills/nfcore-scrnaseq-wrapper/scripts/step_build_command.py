#!/usr/bin/env python3
"""clawbio step: build_command.

Composes the exact `nextflow run ...` command line a real execution would
use: the resolved pipeline source (local checkout path, or `nf-core/scrnaseq
-r <version>`), `-profile`, `-params-file`, and `-work-dir` if given. Ported
unchanged from ClawBio 1.x's ``command_builder.py``. This step never invokes
`subprocess` against Nextflow — it only assembles the argv list and its
shell-quoted string form, which is this DAG's final hand-off artifact for a
human (or a follow-up, out-of-DAG step) to actually run.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

_SCRIPTS_DIR = Path(__file__).resolve().parent
if str(_SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS_DIR))

from command_builder import build_nextflow_command
from textio import write_text_lf
from wrapper_args import parse_optional_str


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Compose the nextflow run command line (does not execute Nextflow)."
    )
    parser.add_argument("--pipeline-source", required=True, help="pipeline_source.json from resolve_pipeline_source")
    parser.add_argument("--params", required=True, help="params.yaml from build_params")
    parser.add_argument("--profile", required=True)
    parser.add_argument("--work-dir", default="none")
    parser.add_argument("--out", required=True, help="Path for the plain-text command")
    parser.add_argument("--out-json", required=True, help="Path for {argv, command_str}")
    args = parser.parse_args()

    pipeline_source = json.loads(Path(args.pipeline_source).read_text(encoding="utf-8"))
    work_dir = parse_optional_str(args.work_dir)

    argv, command_str = build_nextflow_command(
        pipeline_source=pipeline_source,
        profile=args.profile,
        params_path=Path(args.params).resolve(),
        resume=False,
        work_dir=work_dir,
        extra_configs=[],
    )

    write_text_lf(args.out, command_str + "\n")
    Path(args.out_json).write_text(
        json.dumps({"argv": argv, "command_str": command_str}, indent=2), encoding="utf-8"
    )

    print("[build_command] composed (not executed):")
    print(f"  {command_str}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
