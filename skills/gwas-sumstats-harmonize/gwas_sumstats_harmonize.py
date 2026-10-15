#!/usr/bin/env python3
"""
gwas_sumstats_harmonize.py — harmonize GWAS summary statistics to one canonical
format (SNP CHR BP EA NEA EAF BETA SE P N), optionally aligned to a reference.

The work is done by a Snakemake workflow in ./workflow (config-by-concern
layout: data/analysis/software.yaml + one rule and one standalone script per
stage). This wrapper writes the run's config to <output>/config/, then runs the
workflow with Snakemake when available, or runs the same stage scripts in
order with plain Python (identical commands, identical output).

Usage:
    python gwas_sumstats_harmonize.py --input sumstats.tsv.gz --output out/ [--reference ref.tsv]
    python gwas_sumstats_harmonize.py --input manifest.yaml --output out/
    python gwas_sumstats_harmonize.py --demo --output /tmp/harmonize_demo
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
from datetime import date
from pathlib import Path

import yaml

_SKILL_DIR = Path(__file__).resolve().parent
_PROJECT_ROOT = _SKILL_DIR.parents[1]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from clawbio.common.reproducibility import (  # noqa: E402
    write_checksums,
    write_commands_sh,
    write_environment_yml,
)

DISCLAIMER = (
    "ClawBio is a research and educational tool. It is not a medical device "
    "and does not provide clinical diagnoses. Consult a healthcare "
    "professional before making any medical decisions."
)

SKILL_NAME = "gwas-sumstats-harmonize"
SKILL_VERSION = "0.1.0"
WORKFLOW = _SKILL_DIR / "workflow"
SCRIPTS = WORKFLOW / "scripts"
DEMO_MANIFEST = _SKILL_DIR / "examples" / "demo_manifest.yaml"
STAGES = ["map_columns", "derive_effects", "qc_filter", "align_reference"]
_NAME_RE = re.compile(r"[^A-Za-z0-9_-]+")


class RunError(RuntimeError):
    """A stage failed or the inputs could not be set up."""


def _load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_lib = _load_module(SCRIPTS / "lib" / "harmonize_lib.py", "gwas_sumstats_harmonize_lib")
_config_loader = _load_module(WORKFLOW / "config" / "config_loader.py", "gwas_sumstats_harmonize_config")


# ── Inputs -> run config ──────────────────────────────────────────────────────


def _dataset_name(path: Path) -> str:
    name = _NAME_RE.sub("_", path.name.split(".")[0]).strip("_")
    return name or "sumstats"


def _read_mapping_file(path: str | None) -> dict:
    if not path:
        return {}
    data = yaml.safe_load(Path(path).read_text(encoding="utf-8")) or {}
    if not isinstance(data, dict):
        raise RunError(f"--columns file {path} must be a mapping of canonical name -> header name")
    return {str(k).upper(): str(v) for k, v in data.items()}


def collect_inputs(args) -> tuple[dict, str | None]:
    """Return ({name: {path, label, columns}}, reference_path)."""
    if args.demo:
        manifest = DEMO_MANIFEST
    elif args.input and Path(args.input).suffix.lower() in (".yaml", ".yml"):
        manifest = Path(args.input)
    else:
        manifest = None

    if manifest is not None:
        data = yaml.safe_load(manifest.read_text(encoding="utf-8")) or {}
        base = manifest.resolve().parent
        datasets = {}
        for name, entry in (data.get("datasets") or {}).items():
            if _NAME_RE.search(str(name)):
                raise RunError(f"dataset name '{name}' may contain only letters, digits, '_' or '-'")
            entry = dict(entry or {})
            path = Path(entry["path"])
            datasets[name] = {
                "path": str(path if path.is_absolute() else (base / path).resolve()),
                "label": entry.get("label", name),
                "columns": {str(k).upper(): str(v) for k, v in (entry.get("columns") or {}).items()},
            }
        if not datasets:
            raise RunError(f"manifest {manifest} lists no datasets")
        ref = data.get("reference")
        if ref and not Path(ref).is_absolute():
            ref = str((base / ref).resolve())
        return datasets, (args.reference or ref)

    if not args.input:
        raise RunError("give --input <sumstats file or manifest.yaml>, or --demo")
    path = Path(args.input).resolve()
    if not path.exists():
        raise RunError(f"input not found: {path}")
    name = _dataset_name(path)
    return {name: {"path": str(path), "label": name, "columns": _read_mapping_file(args.columns)}}, args.reference


def write_run_config(output_dir: Path, datasets: dict, reference: str | None, args) -> Path:
    config_dir = output_dir / "config"
    config_dir.mkdir(parents=True, exist_ok=True)
    data = {"datasets": datasets}
    analysis = {
        "output_dir": (output_dir / "work").as_posix(),
        "targets": list(datasets),
        "resources": {"reference": Path(reference).resolve().as_posix() if reference else None},
        "map_columns": {},
        "derive_effects": {},
        "qc_filter": {"palindromic": args.palindromic, "min_maf": args.min_maf, "keep_indels": not args.drop_indels},
        "align_reference": {"drop_unmatched": not args.keep_unmatched},
    }
    software = {"threads": args.cores}
    header = "# Written by gwas_sumstats_harmonize.py for this run; read only by workflow/config/config_loader.py\n"
    for fname, payload in (("data.yaml", data), ("analysis.yaml", analysis), ("software.yaml", software)):
        (config_dir / fname).write_text(header + yaml.safe_dump(payload, sort_keys=False), encoding="utf-8")
    return config_dir


# ── Engines ───────────────────────────────────────────────────────────────────


def _bool(v) -> str:
    return str(v).lower()


def stage_commands(cfg: dict, name: str) -> list[tuple[str, list[str]]]:
    """The exact per-stage commands, mirroring workflow/rules/*.smk."""
    out = cfg["analysis"]["output_dir"]
    qc = cfg["analysis"]["qc_filter"]
    al = cfg["analysis"]["align_reference"]
    ref = cfg["analysis"]["resources"]["reference"]
    py = sys.executable
    s = SCRIPTS.as_posix()
    cmds = [
        ("map_columns", [py, f"{s}/map_columns.py", "--input", cfg["data"]["datasets"][name]["path"],
                         "--columns-json", f"{out}/tmp/{name}.columns.json",
                         "--out", f"{out}/map_columns/{name}.tsv",
                         "--summary-json", f"{out}/map_columns/{name}.summary.json"]),
        ("derive_effects", [py, f"{s}/derive_effects.py", "--input", f"{out}/map_columns/{name}.tsv",
                            "--out", f"{out}/derive_effects/{name}.tsv",
                            "--summary-json", f"{out}/derive_effects/{name}.summary.json"]),
        ("qc_filter", [py, f"{s}/qc_filter.py", "--input", f"{out}/derive_effects/{name}.tsv",
                       "--out", f"{out}/qc_filter/{name}.tsv",
                       "--summary-json", f"{out}/qc_filter/{name}.summary.json",
                       "--palindromic", str(qc["palindromic"]), "--min-maf", str(qc["min_maf"]),
                       "--keep-indels", _bool(qc["keep_indels"])]),
        ("align_reference", [py, f"{s}/align_reference.py", "--input", f"{out}/qc_filter/{name}.tsv"]
         + (["--reference", ref] if ref else [])
         + ["--out", f"{out}/align_reference/{name}.tsv",
            "--summary-json", f"{out}/align_reference/{name}.summary.json",
            "--drop-unmatched", _bool(al["drop_unmatched"])]),
    ]
    return cmds


def _write_columns_json(cfg: dict) -> None:
    out = cfg["analysis"]["output_dir"]
    os.makedirs(f"{out}/tmp", exist_ok=True)
    for name, entry in cfg["data"]["datasets"].items():
        path = f"{out}/tmp/{name}.columns.json"
        text = json.dumps(entry["columns"], indent=2, sort_keys=True)
        if not os.path.exists(path) or Path(path).read_text() != text:
            Path(path).write_text(text)


def _log_tail(path: str, n: int = 8) -> str:
    try:
        return "\n".join(Path(path).read_text(encoding="utf-8", errors="replace").splitlines()[-n:])
    except OSError:
        return ""


def run_python_engine(cfg: dict) -> None:
    out = cfg["analysis"]["output_dir"]
    _write_columns_json(cfg)
    for name in cfg["analysis"]["targets"]:
        for stage, cmd in stage_commands(cfg, name):
            log = f"{out}/logs/{stage}/{name}.log"
            os.makedirs(os.path.dirname(log), exist_ok=True)
            with open(log, "w") as fh:
                proc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT)
            if proc.returncode != 0:
                raise RunError(f"stage {stage} failed for dataset '{name}':\n{_log_tail(log)}")
            done = Path(f"{out}/done/{stage}_{name}.done")
            done.parent.mkdir(parents=True, exist_ok=True)
            done.touch()


def run_snakemake_engine(cfg: dict, output_dir: Path, config_dir: Path, cores: int) -> None:
    cmd = [sys.executable, "-m", "snakemake", "-s", str(WORKFLOW / "Snakefile"), "-d", str(output_dir),
           "--cores", str(cores), "--config", f"config_dir={config_dir.as_posix()}"]
    if os.name == "nt":
        cmd.append("--drop-metadata")  # .snakemake/metadata names can overflow MAX_PATH on Windows
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        out = cfg["analysis"]["output_dir"]
        details = []
        for stage in STAGES:
            for name in cfg["analysis"]["targets"]:
                if not Path(f"{out}/done/{stage}_{name}.done").exists():
                    tail = _log_tail(f"{out}/logs/{stage}/{name}.log")
                    if tail:
                        details.append(f"stage {stage} failed for dataset '{name}':\n{tail}")
        raise RunError("\n".join(details) or (proc.stdout + proc.stderr)[-3000:])


def choose_engine(requested: str) -> str:
    if requested != "auto":
        return requested
    return "snakemake" if importlib.util.find_spec("snakemake") else "python"


# ── Results ───────────────────────────────────────────────────────────────────


def _summary(out: str, stage: str, name: str) -> dict:
    return json.loads(Path(f"{out}/{stage}/{name}.summary.json").read_text())


def collect_results(cfg: dict, output_dir: Path) -> dict:
    out = cfg["analysis"]["output_dir"]
    harmonized = output_dir / "harmonized"
    harmonized.mkdir(exist_ok=True)
    results = {}
    for name in cfg["analysis"]["targets"]:
        final = harmonized / f"{name}.tsv"
        shutil.copyfile(f"{out}/align_reference/{name}.tsv", final)
        m = _summary(out, "map_columns", name)
        d = _summary(out, "derive_effects", name)
        q = _summary(out, "qc_filter", name)
        a = _summary(out, "align_reference", name)
        pvals = [row["P"] for row in _lib.read_table(str(final))]
        alignment = dict(a["alignment"])
        dropped = dict(m.get("dropped", {}))
        dropped.update(q["dropped"])
        if alignment.get("mismatch"):
            dropped["allele_mismatch"] = alignment["mismatch"]
        if alignment.get("palindromic_unresolved"):
            dropped["palindromic_unresolved"] = alignment["palindromic_unresolved"]
        if a["params"]["drop_unmatched"] and alignment.get("not_in_reference"):
            dropped["not_in_reference"] = alignment["not_in_reference"]
        results[name] = {
            "label": cfg["data"]["datasets"][name]["label"],
            "input": cfg["data"]["datasets"][name]["path"],
            "output": final.relative_to(output_dir).as_posix(),
            "rows_in": m["rows_in"],
            "rows_out": a["rows_out"],
            "build": m.get("build") or cfg_build(cfg),
            "mapping": m["mapping"],
            "notes": m["notes"],
            "derived": d["derived"],
            "dropped": dropped,
            "alignment": alignment,
            "lambda_gc": round(_lib.lambda_gc(pvals), 4) if pvals else None,
        }
    return results


def cfg_build(cfg: dict) -> str | None:
    return cfg.get("declared_build")


def check_builds(results: dict, declared: str | None) -> None:
    """A declared build that contradicts a build-labelled column is an error, not a guess."""
    for name, v in results.items():
        if declared and v["build"] and v["build"] != declared:
            raise RunError(f"dataset '{name}': --build {declared} contradicts its column names, "
                           f"which say {v['build']} ({v['mapping'].get('BP')})")


def write_summary_table(output_dir: Path, results: dict) -> Path:
    path = output_dir / "tables" / "harmonize_summary.tsv"
    path.parent.mkdir(exist_ok=True)
    reasons = sorted({r for v in results.values() for r in v["dropped"]})
    cols = ["dataset", "rows_in", "rows_out", "lambda_gc"] + [f"dropped_{r}" for r in reasons] + [
        "swapped", "strand_flipped", "eaf_filled", "not_in_reference"]
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for name, v in results.items():
            al = v["alignment"]
            w.writerow([name, v["rows_in"], v["rows_out"], v["lambda_gc"]]
                       + [v["dropped"].get(r, 0) for r in reasons]
                       + [al.get("swapped", 0) + al.get("strand_flipped_swapped", 0),
                          al.get("strand_flipped", 0) + al.get("strand_flipped_swapped", 0),
                          al.get("eaf_filled", 0), al.get("not_in_reference", 0)])
    return path


def write_report(output_dir: Path, results: dict, reference: str | None, engine: str, cfg: dict) -> Path:
    qc = cfg["analysis"]["qc_filter"]
    drop_unmatched = cfg["analysis"]["align_reference"]["drop_unmatched"]
    unaligned = sum(v["alignment"].get("not_in_reference", 0) for v in results.values())
    if not reference:
        ea_note = "."
    elif drop_unmatched or not unaligned:
        ea_note = " (EA = reference ALT)."
    else:
        ea_note = (f". EA = reference ALT for matched variants; {unaligned} variant(s) absent from the "
                   "reference were kept as reported and are not aligned (--keep-unmatched).")
    builds = sorted({v["build"] for v in results.values() if v["build"]})
    lines = [
        "# GWAS Summary Statistics Harmonization Report",
        "",
        f"**Date**: {date.today().isoformat()}  ",
        f"**Datasets**: {len(results)} · **Engine**: {engine} · "
        f"**Reference**: {Path(reference).name if reference else 'none (alleles not aligned)'}  ",
        f"**QC**: palindromic={qc['palindromic']}, min_maf={qc['min_maf']}, keep_indels={qc['keep_indels']}  ",
        f"**Genome build**: {', '.join(builds) if builds else 'not declared and not in the column names; confirm it matches the reference'}",
        "",
        "Canonical output columns: `SNP CHR BP EA NEA EAF BETA SE P N`" + ea_note,
        "",
        "## Summary",
        "",
        "| Dataset | Rows in | Rows out | Dropped | Swapped | Strand-flipped | EAF filled | λGC |",
        "|---------|--------:|---------:|--------:|--------:|---------------:|-----------:|----:|",
    ]
    for name, v in results.items():
        al = v["alignment"]
        lines.append(
            f"| `{name}` | {v['rows_in']} | {v['rows_out']} | {v['rows_in'] - v['rows_out']} | "
            f"{al.get('swapped', 0) + al.get('strand_flipped_swapped', 0)} | "
            f"{al.get('strand_flipped', 0) + al.get('strand_flipped_swapped', 0)} | "
            f"{al.get('eaf_filled', 0)} | {v['lambda_gc']} |"
        )
    lines += ["", "## Per dataset", ""]
    for name, v in results.items():
        mapping = ", ".join(f"{k}←`{h}`" for k, h in v["mapping"].items())
        drops = ", ".join(f"{k}: {n}" for k, n in v["dropped"].items()) or "none"
        derived = ", ".join(f"{k}: {n}" for k, n in v["derived"].items()) or "none"
        lines += [f"### {name} ({v['label']})", "", f"- **Columns**: {mapping}",
                  f"- **Derived**: {derived}", f"- **Dropped**: {drops}"]
        if v["alignment"].get("not_in_reference"):
            lines.append(f"- **Not in reference**: {v['alignment']['not_in_reference']}")
        lines += [f"- **Note**: {n}" for n in v["notes"]]
        lines.append("")
    lines += [
        "## Interpretation notes",
        "",
        "- λGC is computed on the harmonized variants only. On a small or pre-filtered file it is",
        "  not an estimate of genome-wide inflation; use LD score regression for that.",
        "- Palindromic A/T and C/G variants are strand-resolved by comparing EAF with the reference",
        "  ALT frequency, and dropped (palindromic_unresolved) when either is missing or within 0.4-0.6.",
        "  Without a reference they are passed through as reported, not strand-checked.",
        "- BETA_SE_from_Z, when present, is on the standardised (per-SD) scale, not the original trait scale.",
        "- p-values below the float64 minimum (about 1e-308) are kept as exact strings.",
        "",
        "---",
        "",
        f"*{DISCLAIMER}*",
        "",
    ]
    path = output_dir / "report.md"
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


# ── CLI ───────────────────────────────────────────────────────────────────────


def parse_args(argv=None):
    p = argparse.ArgumentParser(description="Harmonize GWAS summary statistics to a canonical format.")
    p.add_argument("--input", help="Summary statistics file, or a manifest .yaml listing several datasets")
    p.add_argument("--output", required=True, help="Output directory")
    p.add_argument("--demo", action="store_true", help="Run on four synthetic cohorts in four formats")
    p.add_argument("--reference", help="Reference TSV (CHR BP REF ALT [AF]) to align alleles to")
    p.add_argument("--columns", help="YAML/JSON {canonical: header_name} mapping for a single --input file")
    p.add_argument("--palindromic", choices=["ambiguous", "all", "none"], default="ambiguous",
                   help="Drop palindromic SNVs: ambiguous (EAF missing or 0.4-0.6, default), all, none")
    p.add_argument("--min-maf", type=float, default=0.0, help="Drop variants with MAF below this (default 0)")
    p.add_argument("--drop-indels", action="store_true", help="Drop multi-base alleles")
    p.add_argument("--keep-unmatched", action="store_true",
                   help="Keep variants absent from the reference, unaligned (default: drop them)")
    p.add_argument("--build", choices=["GRCh37", "GRCh38"],
                   help="Genome build of the input positions; recorded, and checked against build-labelled columns")
    p.add_argument("--engine", choices=["auto", "snakemake", "python"], default="auto",
                   help="Run with Snakemake (if installed) or plain Python; output is identical")
    p.add_argument("--cores", type=int, default=1, help="Cores for the Snakemake engine")
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    output_dir = Path(args.output).resolve()
    try:
        datasets, reference = collect_inputs(args)
        output_dir.mkdir(parents=True, exist_ok=True)
        config_dir = write_run_config(output_dir, datasets, reference, args)
        cfg = _config_loader.load_config(str(config_dir))
        cfg["declared_build"] = args.build
        engine = choose_engine(args.engine)
        if engine == "snakemake":
            run_snakemake_engine(cfg, output_dir, config_dir, args.cores)
        else:
            run_python_engine(cfg)
        results = collect_results(cfg, output_dir)
        check_builds(results, args.build)
    except (RunError, _config_loader.ConfigError, OSError, yaml.YAMLError, KeyError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    table = write_summary_table(output_dir, results)
    report = write_report(output_dir, results, cfg["analysis"]["resources"]["reference"], engine, cfg)
    result_path = output_dir / "result.json"
    result_path.write_text(json.dumps({
        "skill": SKILL_NAME,
        "version": SKILL_VERSION,
        "engine": engine,
        "reference": cfg["analysis"]["resources"]["reference"],
        "build": args.build,
        "settings": {"qc_filter": cfg["analysis"]["qc_filter"], "align_reference": cfg["analysis"]["align_reference"]},
        "datasets": results,
    }, indent=2), encoding="utf-8")

    invocation = "python skills/gwas-sumstats-harmonize/gwas_sumstats_harmonize.py " + shlex.join(
        sys.argv[1:] if argv is None else list(argv))
    stage_lines = []
    for name in cfg["analysis"]["targets"]:
        for stage, cmd in stage_commands(cfg, name):
            stage_lines.append(f"# {stage} / {name}\n" + shlex.join(cmd))
    write_commands_sh(output_dir, invocation + "\n\n# Equivalent per-stage commands:\n" + "\n".join(stage_lines))
    write_environment_yml(output_dir, env_name="clawbio-gwas-sumstats-harmonize",
                          pip_deps=["pyyaml>=6.0", "snakemake>=8"])
    harmonized = [output_dir / v["output"] for v in results.values()]
    write_checksums([report, result_path, table] + harmonized, output_dir, anchor=output_dir)

    print(f"Harmonized {len(results)} dataset(s) with the {engine} engine -> {output_dir / 'harmonized'}")
    print(f"Report: {report}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
