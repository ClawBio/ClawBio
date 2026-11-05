#!/usr/bin/env python3
# nfcore_sarek_wrapper / scripts/build_command.py
"""Stage 4 of the pipeline.yaml DAG: compose the Nextflow invocation.

Thin CLI wrapper around the faithfully-ported ``command_builder`` module (same
directory, unchanged from ClawBio 1.x): ``compose_profile()`` +
``build_nextflow_command()``. Writes the exact ``nextflow run ...`` command
line a real run would execute, plus a preview ``commands.sh`` — but this
script never invokes ``subprocess.run`` on it. That is the line this port
draws around "external": actually running it needs Nextflow, a container
runtime, and reference genome data, none of which this sandbox has (see
SKILL.md). Everything up to and including "what exact command would run" is
real, sandbox-verified Python logic.
"""
from __future__ import annotations

import argparse
import json
import stat
import sys
from pathlib import Path

_SCRIPT_DIR = Path(__file__).resolve().parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))

from command_builder import build_nextflow_command, compose_profile
from pipeline_source import resolve_pipeline_source


def _truthy(raw: str) -> bool:
    return str(raw).strip().lower() in ("1", "true", "yes", "on")


def _none_if_sentinel(raw: str) -> str | None:
    v = (raw or "").strip()
    return None if not v or v.lower() == "none" else v


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Compose (do not run) the Sarek Nextflow command")
    parser.add_argument("--params-yaml", required=True)
    parser.add_argument("--profile", required=True)
    parser.add_argument("--demo", required=True)
    parser.add_argument("--arm", required=True)
    parser.add_argument("--gpu", required=True)
    parser.add_argument("--spark-profile", required=True)
    parser.add_argument("--mutect-profile", required=True)
    parser.add_argument("--pipeline-version", required=True)
    parser.add_argument("--pipeline-local", required=True)
    parser.add_argument("--out-command-json", required=True)
    parser.add_argument("--out-commands-sh", required=True)
    args = parser.parse_args(argv)

    demo = _truthy(args.demo)
    arm = _truthy(args.arm)
    gpu = _truthy(args.gpu)
    spark = _truthy(args.spark_profile)
    mutect = _truthy(args.mutect_profile)

    pipeline_local = _none_if_sentinel(args.pipeline_local)
    source = resolve_pipeline_source(
        requested_version=args.pipeline_version,
        local_pipeline_dir=Path(pipeline_local).expanduser().resolve() if pipeline_local else None,
    )

    composed_profile = compose_profile(
        args.profile, demo=demo, arm=arm, gpu=gpu, spark=spark, mutect_profile=mutect,
    )

    command_list, command_str = build_nextflow_command(
        pipeline_source=source,
        profile=args.profile,
        params_path=Path(args.params_yaml).resolve(),
        resume=False,
        work_dir=None,
        extra_configs=None,
        demo=demo,
        arm=arm,
        gpu=gpu,
        spark=spark,
        mutect_profile=mutect,
    )

    out_json = Path(args.out_command_json)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(
        json.dumps(
            {
                "command": command_list,
                "command_str": command_str,
                "composed_profile": composed_profile,
                "pipeline_source": source,
                "executed": False,
                "note": (
                    "This command was composed but NOT executed. Running it for real "
                    "additionally requires Nextflow >=25.10.2, a container runtime "
                    "(Docker/Singularity/...), and — unless --demo — local reference "
                    "genome data. See SKILL.md 'What a real run additionally needs'."
                ),
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    sh_path = Path(args.out_commands_sh)
    sh_path.parent.mkdir(parents=True, exist_ok=True)
    sh_path.write_text(
        "#!/usr/bin/env bash\n"
        "# Composed by nfcore_sarek_wrapper's build_command stage. NOT executed by\n"
        "# clawbio: this DAG validates and prepares a real nf-core/sarek 3.8.1 run\n"
        "# without launching Nextflow. Requires Nextflow, a container runtime, Java\n"
        "# >=17, and (unless --demo) local reference genome data before running.\n"
        "set -euo pipefail\n"
        f"cd \"$(dirname \"${{BASH_SOURCE[0]}}\")\"\n"
        f"{command_str}\n",
        encoding="utf-8",
    )
    sh_path.chmod(sh_path.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)

    print(f"[build_command] composed (not executed): {command_str}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
