#!/usr/bin/env python3
"""DAG step: build_command.

Thin CLI wrapper around the ported, unmodified ``command_builder.build_nextflow_command``:
composes the exact ``nextflow run ...`` invocation a real run would execute
(profile string, -params-file, -work-dir, -resume). This step only WRITES that
command to a file -- it never invokes ``nextflow``. Actually launching it
requires Nextflow, a container/env runtime, and (for non-demo runs) reference
genome data, none of which this sandbox has; see SKILL.md Dependencies for
the documented, out-of-DAG execution step.
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

purge_foreign_bare_modules("errors", "schemas", "command_builder", "args_common")

from args_common import build_args, none_or  # noqa: E402
from command_builder import build_nextflow_command  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description="Compose (but do not run) the nf-core/rnaseq Nextflow command.")
    parser.add_argument("--output", required=True)
    parser.add_argument("--demo", required=True)
    parser.add_argument("--profile", required=True)
    parser.add_argument("--prokaryotic", required=True)
    parser.add_argument("--work-dir", required=True, help="Nextflow work dir, or the sentinel 'none' for the default <output>/upstream/work.")
    parser.add_argument("--pipeline-source", required=True)
    parser.add_argument("--params-yaml", required=True)
    parser.add_argument("--out-txt", required=True)
    parser.add_argument("--out-json", required=True)
    parsed = parser.parse_args()

    args = build_args(
        output=parsed.output,
        demo=parsed.demo,
        profile=parsed.profile,
        prokaryotic=parsed.prokaryotic,
        work_dir=parsed.work_dir,
    )
    pipeline_source = json.loads(Path(parsed.pipeline_source).read_text(encoding="utf-8"))
    output_dir = Path(parsed.output).expanduser().resolve()
    work_dir = Path(args.work_dir) if args.work_dir else (output_dir / "upstream" / "work")

    command, command_str = build_nextflow_command(
        pipeline_source=pipeline_source,
        profile=args.profile,
        params_path=Path(parsed.params_yaml),
        resume=False,
        work_dir=work_dir,
        extra_configs=[],
        demo=bool(args.demo),
        prokaryotic=bool(args.prokaryotic),
        rapid_quant=False,
        arm=False,
    )

    # Nextflow launches with cwd=<output> in the real wrapper (see schemas.py /
    # SKILL.md "Results are written under a relative upstream/results"). Record
    # that alongside the command so the replay recipe is complete without
    # guessing the launch directory.
    replay = {
        "cwd": str(output_dir),
        "command": command,
        "command_str": command_str,
        "note": (
            "NOT executed by this DAG. To actually run it: cd into 'cwd', ensure "
            "nextflow/java/a container runtime are installed and (for docker) the "
            "daemon is running, then run 'command_str' from that directory."
        ),
    }

    out_txt = Path(parsed.out_txt)
    out_txt.parent.mkdir(parents=True, exist_ok=True)
    out_txt.write_text(f"cd {output_dir}\n{command_str}\n", encoding="utf-8")

    out_json = Path(parsed.out_json)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(replay, indent=2) + "\n", encoding="utf-8")

    print(f"[build_command] wrote {out_txt} (not executed):")
    print(f"  cd {output_dir}")
    print(f"  {command_str}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
