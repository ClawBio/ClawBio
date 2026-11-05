#!/usr/bin/env python3
"""Stage 3: substitution model selection.

`--model auto` (the pipeline.yaml sentinel for 1.x's "no --model given")
runs IQ-TREE2 ModelFinder (`-m MFP`) and parses `Best-fit model according to
BIC:` out of the `.iqtree` log. Any other `--model` value skips ModelFinder
and uses that model directly — exactly 1.x's `if not model:` branch. If
ModelFinder isn't available (iqtree2/iqtree missing) or fails, this falls
back to `GTR+G`, matching 1.x — it never hard-fails.
"""
from __future__ import annotations

import argparse
import json
import sys
import tempfile
from pathlib import Path

import common


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True, type=Path, help="Trimmed (or aligned) FASTA")
    p.add_argument("--model", required=True, help="Substitution model, or 'auto' to run ModelFinder")
    p.add_argument("--out", required=True, type=Path, help="Output model file (plain text)")
    p.add_argument("--out-meta", required=True, type=Path, help="Output stage-metadata JSON path")
    args = p.parse_args()

    args.out.parent.mkdir(parents=True, exist_ok=True)

    if not args.input.exists():
        print(f"Error: input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    user_model = None if args.model.strip().lower() == common.SENTINEL_AUTO else args.model.strip()
    print(f"[select_model] input={args.input} model={args.model}")

    if user_model:
        model = user_model
        stage = "model:user-provided"
    else:
        iqtree_bin = common.iqtree_binary()
        if iqtree_bin:
            with tempfile.TemporaryDirectory() as tmpdir:
                mf_prefix = Path(tmpdir) / "modelfinder"
                try:
                    model = common.run_modelfinder(iqtree_bin, args.input, mf_prefix)
                    stage = f"modelfinder:{model}"
                except Exception as exc:
                    model = "GTR+G"
                    stage = f"modelfinder:failed({exc.__class__.__name__})-using-GTR+G"
                    print(f"  Warning: ModelFinder failed ({exc}). Using GTR+G.", file=sys.stderr)
        else:
            model = "GTR+G"
            stage = "model:default-GTR+G(no-iqtree)"

    args.out.write_text(model + "\n")
    meta = {"stage": stage, "model": model}
    args.out_meta.write_text(json.dumps(meta, indent=2))
    print(f"  [{stage}] -> {args.out}")


if __name__ == "__main__":
    main()
