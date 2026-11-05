#!/usr/bin/env python3
"""Stage 2: alignment trimming via trimAl -automated1.

Never hard-fails, exactly like 1.x: if `--trim false` (the 1.x `--no-trim`
flag) the alignment passes through untouched; if trimAl isn't on PATH the
stage warns and passes through untouched; if trimAl itself errors the stage
warns and passes through untouched. Only a missing/unreadable input file is
a hard error.
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

import common

TRIMAL_BIN = shutil.which("trimal")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True, type=Path, help="Aligned FASTA from the align step")
    p.add_argument("--trim", required=True, help="'false' to skip trimming (1.x --no-trim)")
    p.add_argument("--out", required=True, type=Path, help="Output (trimmed) FASTA path")
    p.add_argument("--out-meta", required=True, type=Path, help="Output stage-metadata JSON path")
    args = p.parse_args()

    args.out.parent.mkdir(parents=True, exist_ok=True)

    if not args.input.exists():
        print(f"Error: input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    trim_flag = common.as_bool(args.trim)
    print(f"[trim] input={args.input} trim={trim_flag}")

    trimmed = False
    stage = ""

    if not trim_flag:
        shutil.copy(args.input, args.out)
        stage = "trim:skipped(--no-trim)"
    else:
        if TRIMAL_BIN is None:
            shutil.copy(args.input, args.out)
            stage = "trim:skipped(trimal-not-found)"
        else:
            try:
                common.run_trimal(args.input, args.out)
                trimmed = True
                stage = "trim:trimal"
            except Exception as exc:
                print(
                    f"Warning: trimAl failed ({exc}). Continuing without trimming.",
                    file=sys.stderr,
                )
                shutil.copy(args.input, args.out)
                stage = f"trim:failed({exc.__class__.__name__})"

    meta = {"stage": stage, "trimmed": trimmed}
    args.out_meta.write_text(json.dumps(meta, indent=2))
    print(f"  [{stage}] -> {args.out}")


if __name__ == "__main__":
    main()
