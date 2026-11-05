#!/usr/bin/env python3
# nfcore_sarek_wrapper / scripts/build_params.py
"""Stage 3 of the pipeline.yaml DAG: assemble the effective Sarek params.yaml.

Thin CLI wrapper around the faithfully-ported ``params_builder`` module (same
directory): ``build_effective_params()`` + ``write_params_yaml()``, unchanged
from ClawBio 1.x except for the ``write_text_lf`` import (see the header note
in params_builder.py). Builds an ``argparse.Namespace`` shaped the way
``params_builder.build_effective_params(args, ...)`` expects (it reads
``getattr(args, field, None)`` for a long list of Sarek passthrough fields —
see 1.x's ``_SAREK_PASSTHROUGH_PARAMS``); fields this port does not expose as
first-class pipeline.yaml inputs simply stay unset (None/False), matching
their upstream schema defaults. The ``--extra-params-json`` escape hatch
(shared with the preflight step) lets a caller supply any additional native
Sarek param this port has not wired a dedicated flag for — the same role
1.x's repeatable ``--extra-param key=value`` played.

Only reached after `preflight` has passed (prereqs gate this in
pipeline.yaml), so params.yaml written here corresponds to an
already-validated run configuration.
"""
from __future__ import annotations

import argparse
import json
import sys
from argparse import Namespace
from pathlib import Path

_SCRIPT_DIR = Path(__file__).resolve().parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))

from errors import SkillError
from params_builder import build_effective_params, serialize_params_yaml, write_params_yaml


def _truthy(raw: str) -> bool:
    return str(raw).strip().lower() in ("1", "true", "yes", "on")


def _none_if_sentinel(raw: str) -> str | None:
    v = (raw or "").strip()
    return None if not v or v.lower() == "none" else v


def _tool_list(raw: str) -> str | None:
    v = _none_if_sentinel(raw)
    if not v:
        return None
    return ",".join(t.strip().lower() for t in v.split(",") if t.strip())


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Build the effective Sarek params.yaml")
    parser.add_argument("--samplesheet-report", required=True)
    parser.add_argument("--step", required=True)
    parser.add_argument("--tools", required=True)
    parser.add_argument("--skip-tools", required=True)
    parser.add_argument("--aligner", required=True)
    parser.add_argument("--genome", required=True)
    parser.add_argument("--igenomes-ignore", required=True)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--wes", required=True)
    parser.add_argument("--demo", required=True)
    parser.add_argument("--extra-params-json", required=True)
    parser.add_argument("--out-dir", required=True, help="Directory to write reproducibility/params.yaml under")
    parser.add_argument("--out", required=True, help="Copy of the written params.yaml path, for pipeline.yaml emits")
    args = parser.parse_args(argv)

    samplesheet_report = json.loads(Path(args.samplesheet_report).read_text(encoding="utf-8"))
    normalized_path_raw = samplesheet_report.get("normalized_path")

    ns = Namespace()
    ns.step = args.step
    ns.demo = _truthy(args.demo)
    ns._noinput = False
    ns.build_only_index = False
    ns.aligner = _none_if_sentinel(args.aligner)
    ns.use_gatk_spark = None
    ns.genome = _none_if_sentinel(args.genome)
    ns.igenomes_base = None
    ns.igenomes_ignore = _truthy(args.igenomes_ignore)
    ns.fasta = _none_if_sentinel(args.fasta)
    ns.wes = _truthy(args.wes)
    ns.tools = _tool_list(args.tools)
    ns.skip_tools = _tool_list(args.skip_tools)
    ns.params_file = None
    ns.outdir_cache = None

    # Every other REFERENCE_PATH_FIELDS / ANNOTATION_PATH_FIELDS / boolean /
    # numeric / string field params_builder inspects via getattr(..., None)
    # simply stays absent from this Namespace, which is equivalent to None —
    # not exposing a dedicated CLI flag for every one of the 1.x wrapper's 150+
    # passthrough params (see SKILL.md) is exactly what --extra-param /
    # --extra-params-json exists to cover.
    extras: dict[str, object] = {}
    extras_path = _none_if_sentinel(args.extra_params_json)
    if extras_path:
        extras.update(json.loads(Path(extras_path).read_text(encoding="utf-8")))
    ns._extras = extras

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if ns.demo or normalized_path_raw is None:
        normalized_samplesheet = out_dir / "no-samplesheet-in-demo-mode.csv"
    else:
        normalized_samplesheet = Path(normalized_path_raw)

    try:
        params = build_effective_params(
            ns, normalized_samplesheet=normalized_samplesheet, output_dir=out_dir,
        )
        params_path = write_params_yaml(params, output_dir=out_dir)
    except SkillError as exc:
        print(json.dumps(exc.to_dict(), indent=2), file=sys.stderr)
        return 1

    # pipeline.yaml declares a single {out} emit path for this step; mirror the
    # canonical reproducibility/params.yaml there too so `clawbio done` can
    # hash a stable, predictable path regardless of --out-dir's location.
    out_copy = Path(args.out)
    out_copy.parent.mkdir(parents=True, exist_ok=True)
    out_copy.write_text(serialize_params_yaml(params), encoding="utf-8")

    print(f"[build_params] wrote {params_path} ({len(params)} keys)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
