#!/usr/bin/env python3
"""clawbio step: build_params.

Composes the effective nf-core/scrnaseq `params.yaml` a real run would pass
via `-params-file`: aligner from `--preset`, `input`/`outdir`, protocol,
reference fields, `igenomes_ignore` auto-set whenever an explicit genome
reference is supplied, and `skip_cellbender`. Ported unchanged from ClawBio
1.x's ``params_builder.py``. Does not require Java, Nextflow, or a container
runtime — pure string/YAML composition, so this step genuinely runs to
completion in any sandbox.
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
from params_builder import build_effective_params, serialize_params_yaml
from textio import write_text_lf
from wrapper_args import build_args, parse_bool, parse_optional_str


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Compose the effective nf-core/scrnaseq params.yaml (no Nextflow invocation)."
    )
    parser.add_argument("--samplesheet", required=True, help="Normalized samplesheet CSV from validate_samplesheet")
    parser.add_argument("--preset", required=True)
    parser.add_argument("--protocol", default="none")
    parser.add_argument("--genome", default="none")
    parser.add_argument("--fasta", default="none")
    parser.add_argument("--gtf", default="none")
    parser.add_argument("--skip-cellbender", default="false")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--out", required=True, help="Path for params.yaml")
    parser.add_argument("--out-json", required=True, help="Path for the same params as JSON")
    args = parser.parse_args()

    wrapper_args = build_args(
        preset=args.preset,
        protocol=parse_optional_str(args.protocol),
        genome=parse_optional_str(args.genome),
        fasta=parse_optional_str(args.fasta),
        gtf=parse_optional_str(args.gtf),
        skip_cellbender=parse_bool(args.skip_cellbender),
        demo=False,
    )

    try:
        params = build_effective_params(
            wrapper_args,
            normalized_samplesheet=Path(args.samplesheet).resolve(),
            output_dir=Path(args.output_dir).expanduser().resolve(),
        )
    except SkillError as exc:
        print(json.dumps(exc.to_dict(), indent=2), file=sys.stderr)
        return 1

    write_text_lf(args.out, serialize_params_yaml(params))
    Path(args.out_json).write_text(json.dumps(params, indent=2, default=str), encoding="utf-8")

    print(f"[build_params] wrote params.yaml with {len(params)} key(s): {sorted(params)}")
    print(json.dumps(params, indent=2, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
