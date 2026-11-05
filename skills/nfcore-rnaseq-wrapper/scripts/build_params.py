#!/usr/bin/env python3
"""DAG step: build_params.

Thin CLI wrapper around the ported, unmodified
``params_builder.build_effective_params`` + ``write_params_yaml``: translates
the audited wrapper CLI surface into ``<output>/reproducibility/params.yaml``,
the exact file a real run would hand to Nextflow via ``-params-file``.
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

purge_foreign_bare_modules("errors", "schemas", "params_builder", "args_common")

from args_common import build_args  # noqa: E402
from errors import SkillError  # noqa: E402
from params_builder import build_effective_params, write_params_yaml  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description="Build reproducibility/params.yaml for nf-core/rnaseq.")
    parser.add_argument("--output", required=True)
    parser.add_argument("--demo", required=True)
    parser.add_argument("--input", required=True)
    parser.add_argument("--aligner", required=True)
    parser.add_argument("--genome", required=True)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--gff", required=True)
    parser.add_argument("--profile", required=True)
    parser.add_argument("--prokaryotic", required=True)
    parser.add_argument("--normalized-samplesheet", required=True)
    parser.add_argument("--preflight-result", required=True, help="preflight_result.json (for gencode_autodetected)")
    parser.add_argument("--out", required=True, help="Path to echo the params dict as JSON (provenance)")
    parsed = parser.parse_args()

    args = build_args(
        output=parsed.output,
        demo=parsed.demo,
        input=parsed.input,
        aligner=parsed.aligner,
        genome=parsed.genome,
        fasta=parsed.fasta,
        gtf=parsed.gtf,
        gff=parsed.gff,
        profile=parsed.profile,
        prokaryotic=parsed.prokaryotic,
    )
    output_dir = Path(args.output).expanduser().resolve()
    preflight_result = json.loads(Path(parsed.preflight_result).read_text(encoding="utf-8"))
    gencode_autodetected = bool(preflight_result.get("gencode_autodetected", False))

    try:
        params = build_effective_params(
            args,
            normalized_samplesheet=Path(parsed.normalized_samplesheet),
            output_dir=output_dir,
            gencode_autodetected=gencode_autodetected,
        )
        params_path = write_params_yaml(params, output_dir=output_dir)
    except SkillError as exc:
        out_path = Path(parsed.out)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(json.dumps(exc.to_dict(), indent=2) + "\n", encoding="utf-8")
        print(json.dumps(exc.to_dict(), indent=2), file=sys.stderr)
        return 1

    out_path = Path(parsed.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(params, indent=2) + "\n", encoding="utf-8")
    print(f"[build_params] wrote {params_path}")
    print(json.dumps(params, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
