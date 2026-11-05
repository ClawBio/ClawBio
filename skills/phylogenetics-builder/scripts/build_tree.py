#!/usr/bin/env python3
"""Stage 4: maximum-likelihood tree inference (IQ-TREE2 or RAxML-NG).

Reads the model chosen by select_model.py from `--model-file`. Outgroup
rooting (`-o` in IQ-TREE2, `--outgroup` in RAxML-NG) is applied here, at
build time, not as a separate step — 1.x treats it as a tree-builder flag,
not a post-processing pass (only midpoint rooting is post-processing; see
root.py).

Demo fallback (`--demo true`, mirrors 1.x's MANDATORY "never refuse"
contract): if the engine binary is missing or inference fails, falls back
to the bundled pre-computed 12-taxon primate tree (`examples/demo_tree.nwk`)
instead of erroring — 1.x's `get_demo_fallback()`. Outside demo mode a
missing engine or inference failure is a hard error, exactly like 1.x.
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
    p.add_argument("--model-file", required=True, type=Path, help="Model file from select_model")
    p.add_argument("--engine", required=True, choices=common.VALID_ENGINES)
    p.add_argument("--bootstrap", required=True, choices=common.VALID_BOOTSTRAP)
    p.add_argument("--outgroup", default=common.SENTINEL_NONE, help="Comma-separated outgroup taxa, or 'none'")
    p.add_argument("--threads", required=True, type=int)
    p.add_argument("--seed", required=True, type=int)
    p.add_argument("--demo", default="false", help="'true' to fall back gracefully instead of erroring")
    p.add_argument("--skill-dir", required=True, type=Path)
    p.add_argument("--out", required=True, type=Path, help="Output Newick tree path")
    p.add_argument("--out-meta", required=True, type=Path, help="Output stage-metadata JSON path")
    args = p.parse_args()

    args.out.parent.mkdir(parents=True, exist_ok=True)

    if not args.input.exists():
        print(f"Error: input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)
    if not args.model_file.exists():
        print(f"Error: model file not found: {args.model_file}", file=sys.stderr)
        sys.exit(1)

    model = args.model_file.read_text().strip()
    outgroup = common.none_if_sentinel(args.outgroup)
    demo_flag = common.as_bool(args.demo)

    print(
        f"[build_tree] input={args.input} engine={args.engine} model={model} "
        f"bootstrap={args.bootstrap} outgroup={outgroup} demo={demo_flag}"
    )

    newick = ""
    stage = ""
    engine_used = args.engine

    with tempfile.TemporaryDirectory() as tmpdir:
        tree_prefix = Path(tmpdir) / "tree"

        if args.engine == "iqtree2":
            iqtree_bin = common.iqtree_binary()
            if iqtree_bin:
                try:
                    newick = common.run_iqtree_main(
                        args.input, model, tree_prefix, args.bootstrap, outgroup,
                        args.threads, args.seed,
                    )
                    engine_used = "iqtree2"
                    stage = f"tree:iqtree2:{args.bootstrap}"
                except Exception as exc:
                    if demo_flag:
                        newick, model = common.get_demo_fallback(args.skill_dir)
                        engine_used = "precomputed"
                        stage = "tree:precomputed-fallback"
                        print(f"  [demo fallback] {exc}", file=sys.stderr)
                    else:
                        print(f"Error running IQ-TREE: {exc}", file=sys.stderr)
                        sys.exit(1)
            else:
                if demo_flag:
                    newick, model = common.get_demo_fallback(args.skill_dir)
                    engine_used = "precomputed"
                    stage = "tree:precomputed-fallback(no-iqtree)"
                else:
                    print("Error: IQ-TREE not found in PATH.", file=sys.stderr)
                    sys.exit(1)

        elif args.engine == "raxml-ng":
            import shutil
            raxml_bin = shutil.which("raxml-ng")
            if raxml_bin:
                try:
                    newick = common.run_raxml_main(
                        args.input, model, tree_prefix, args.bootstrap, outgroup,
                        args.threads, args.seed,
                    )
                    engine_used = "raxml-ng"
                    stage = f"tree:raxml-ng:{args.bootstrap}"
                except Exception as exc:
                    if demo_flag:
                        newick, model = common.get_demo_fallback(args.skill_dir)
                        engine_used = "precomputed"
                        stage = "tree:precomputed-fallback"
                        print(f"  [demo fallback] {exc}", file=sys.stderr)
                    else:
                        print(f"Error running RAxML-NG: {exc}", file=sys.stderr)
                        sys.exit(1)
            else:
                if demo_flag:
                    newick, model = common.get_demo_fallback(args.skill_dir)
                    engine_used = "precomputed"
                    stage = "tree:precomputed-fallback(no-raxml-ng)"
                else:
                    print("Error: raxml-ng not found in PATH.", file=sys.stderr)
                    sys.exit(1)

    args.out.write_text(newick + "\n")
    meta = {
        "stage": stage,
        "engine_used": engine_used,
        "model": model,
        "bootstrap": args.bootstrap,
        "outgroup": outgroup,
        "run_mode": "demo-fallback" if engine_used == "precomputed" else "live",
    }
    args.out_meta.write_text(json.dumps(meta, indent=2))
    print(f"  [{stage}] -> {args.out}")


if __name__ == "__main__":
    main()
