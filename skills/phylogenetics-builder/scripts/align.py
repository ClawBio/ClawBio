#!/usr/bin/env python3
"""Stage 1: multiple sequence alignment (MSA).

If `--aligned true`, the input is already an alignment — this stage just
validates equal sequence length and passes it through untouched (the 1.x
script's `effective_aligned` branch). Otherwise it runs the chosen aligner
(mafft default; muscle/clustalw/kalign/tcoffee/prank also supported).

Demo fallback (`--demo true`, mirrors 1.x's MANDATORY "never refuse"
contract): if the aligner binary is missing or the alignment fails, the
stage does not error out — it passes the input through untouched, matching
1.x's `pipeline_steps.append("msa:...-demo-fallback")` behaviour. Outside
demo mode a missing aligner or alignment failure is a hard error, exactly
like 1.x.
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

import common


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True, type=Path, help="Input FASTA (unaligned or pre-aligned)")
    p.add_argument("--aligned", required=True, help="'true' if --input is already an alignment")
    p.add_argument("--aligner", required=True, choices=common.VALID_ALIGNERS)
    p.add_argument("--demo", default="false", help="'true' to fall back gracefully instead of erroring")
    p.add_argument("--skill-dir", required=True, type=Path)
    p.add_argument("--out", required=True, type=Path, help="Output (aligned) FASTA path")
    p.add_argument("--out-meta", required=True, type=Path, help="Output stage-metadata JSON path")
    args = p.parse_args()

    args.out.parent.mkdir(parents=True, exist_ok=True)

    aligned_flag = common.as_bool(args.aligned)
    demo_flag = common.as_bool(args.demo)

    print(f"[align] input={args.input} aligned={aligned_flag} aligner={args.aligner}")

    if not args.input.exists():
        print(f"Error: input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    try:
        sequences = common.parse_fasta(args.input, require_aligned=aligned_flag)
    except Exception as exc:
        print(f"Error validating input: {exc}", file=sys.stderr)
        sys.exit(1)

    stage = ""
    aligner_used = None

    if aligned_flag:
        shutil.copy(args.input, args.out)
        stage = "msa:skipped(pre-aligned)"
    else:
        aligner_used = args.aligner
        try:
            common.run_msa(args.input, args.aligner, args.out)
            stage = f"msa:{args.aligner}"
        except Exception as exc:
            if demo_flag:
                shutil.copy(args.input, args.out)
                stage = f"msa:failed({exc.__class__.__name__})-demo-fallback"
                aligner_used = None
                print(f"  [demo fallback] {exc}", file=sys.stderr)
            else:
                print(f"Error running MSA ({args.aligner}): {exc}", file=sys.stderr)
                sys.exit(1)

    meta = {
        "stage": stage,
        "aligner": aligner_used,
        "aligned": aligned_flag,
        "num_taxa": len(sequences),
    }
    args.out_meta.write_text(json.dumps(meta, indent=2))
    print(f"  [{stage}] taxa={len(sequences)} -> {args.out}")


if __name__ == "__main__":
    main()
