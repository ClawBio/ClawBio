#!/usr/bin/env python3
"""Build a curated, offline research-demo plugin from explicit tracked files.

Nothing is installed or submitted. The archive contains three skills and their
source runtime, never the repository's personal genomes, profiles or eval inputs.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import zipfile

SKILLS = {
    "pharmgx-reporter": ("pharmgx", "Demonstrate pharmacogenomic reporting and missing-genotype abstention with synthetic data."),
    "clinical-variant-reporter": ("acmg", "Demonstrate ACMG-style variant classification and evidence self-audits with the bundled public GIAB-derived panel and cached evidence."),
    "equity-scorer": ("equity", "Measure population representation, heterozygosity and FST in a synthetic cohort; inspect the exploratory HEIM composite and its limitations."),
}

SYNTHETIC_PGX = """# Synthetic 23andMe-format research fixture, GRCh37.
# Constructed genotypes; no person's genome. Deliberately incomplete coverage.
# rsid chromosome position genotype
rs4244285\t10\t96541616\tGG
rs4986893\t10\t96540410\tGG
rs12248560\t10\t96522463\tCC
rs1799853\t10\t96702047\tCC
rs1057910\t10\t96741053\tAA
rs9923231\t16\t31107689\tCT
"""


def encoded(value: object) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode()


def skill_text(name: str, alias: str, description: str) -> bytes:
    return f"""---
name: {name}
description: {description} Use when asked to run or explain a ClawBio research demo. Private genomic input is unsupported by this plugin.
license: MIT
---

# ClawBio research demo: {alias}

1. Explain that this plugin executes bundled demonstration data only. For real
   private data, use a locally installed ClawBio workflow; never request an upload
   to a hosted chat or infer missing genotypes. Do not activate for diagnosis or
   treatment decisions. Do not substitute a narrative for a requested execution.
2. Locate this installed skill's directory. The plugin root is two parents above
   it. Locate `scripts/run_research_demo.py` under that root. Use an available
   Python 3.11+ interpreter. Check `requirements.txt` for dependencies. Explain
   missing execution access or dependencies; do not silently install packages.
3. Choose a fresh output directory in writable workspace, outside the plugin.
   Execute `python <plugin-root>/scripts/run_research_demo.py {alias} --output
   <fresh-output-dir>` using structured process arguments or properly quoted
   paths. Use only this launcher: the retained methodology reference documents
   upstream capabilities beyond this package's supported demonstration scope.
4. On failure, show the error and stop; do not invent results. In particular,
   `PACKAGE_INTEGRITY_FAILED` requires rebuilding from the intended source.
5. Read `report.md`, `result.json` and `plugin_run.json` from that run. Show the
   report and figures if available. Attribute values to the executed workflow.
   Explain missing/indeterminate findings, cached evidence dates, and incomplete
   coverage. Call the HEIM composite exploratory: it is not a validated measure
   of institutional health equity, clinical performance or population identity.
6. Read [upstream methodology](references/methodology.md) only for algorithmic
   detail. Do not treat words such as clinical-grade in upstream prose as proof
   of independent clinical validation. A checksum verifies bytes against this
   package manifest; it does not establish scientific correctness or publisher
   authenticity. Do not invent an overall trust score or safety certification.
7. Include the ClawBio disclaimer with any interpretation: ClawBio is a research
   and educational tool. It is not a medical device and does not provide clinical
   diagnoses. Consult a healthcare professional before making any medical decisions.
""".encode()


def build(root: Path, output: Path) -> Path:
    root = root.resolve()
    tracked = subprocess.check_output(
        ["git", "-C", str(root), "ls-files", "-z"],
    ).decode().split("\0")
    payload: dict[str, bytes] = {}
    origins = {}

    def copy(source: str, target: str | None = None) -> None:
        path = root / source
        if path.is_symlink() or not path.is_file():
            raise ValueError(f"Missing or symlinked source: {source}")
        data = path.read_bytes()
        key = target or source
        payload[key] = data
        origins[key] = {"path": source, "sha256": hashlib.sha256(data).hexdigest()}

    for name in tracked:
        parts = Path(name).parts
        if not parts:
            continue
        if parts[0] == "clawbio" and name.endswith(".py") and "tests" not in parts:
            # Shared source runtime; no package data or genomes.
            copy(name, f"runtime/{name}")
        elif len(parts) >= 3 and parts[0] == "skills" and parts[1] in SKILLS:
            if "tests" in parts or "evals" in parts:
                continue
            if name.endswith(".py"):
                copy(name, f"runtime/{name}")
    for name in SKILLS:
        original = f"skills/{name}/SKILL.md"
        copy(original, f"skills/{name}/references/methodology.md")
        # Fail closed on an upstream licence change before redistributing code.
        if "\nlicense: MIT\n" not in payload[f"skills/{name}/references/methodology.md"].decode():
            raise ValueError(f"Review changed licence before packaging {name}")
        alias, description = SKILLS[name]
        payload[f"skills/{name}/SKILL.md"] = skill_text(name, alias, description)
    for name in ("giab_acmg_panel.vcf", "demo_evidence_cache.json"):
        path = f"skills/clinical-variant-reporter/example_data/{name}"
        copy(path, f"runtime/{path}")
    for name in ("demo_populations.vcf", "demo_population_map.csv"):
        copy(f"examples/{name}", f"runtime/examples/{name}")
    payload["runtime/skills/pharmgx-reporter/demo_patient.txt"] = SYNTHETIC_PGX.encode()
    catalog = json.loads((root / "skills/catalog.json").read_text())
    catalog["skills"] = [entry for entry in catalog["skills"] if entry["name"] in SKILLS]
    if {entry["name"] for entry in catalog["skills"]} != set(SKILLS):
        raise ValueError("Curated skills missing from catalog")
    catalog["skill_count"] = len(catalog["skills"])
    payload["runtime/skills/catalog.json"] = encoded(catalog)
    copy("packaging/openai/run_research_demo.py", "scripts/run_research_demo.py")
    copy("packaging/openai/plugin.json", "plugin.json")
    copy("packaging/openai/requirements.txt", "requirements.txt")
    copy("LICENSE")
    provenance = {
        "schema_version": 1,
        "upstream_repository": "https://github.com/ClawBio/ClawBio",
        "upstream_base_commit": subprocess.check_output(
            ["git", "-C", str(root), "rev-parse", "HEAD"], text=True,
        ).strip(),
        "source_files": origins,
        "private_input_supported": False,
        "datasets": {
            "pharmgx": "Constructed synthetic genotypes; deliberately incomplete coverage",
            "acmg": "Upstream public GIAB-derived demonstration panel; cached evidence, not live database updates",
            "equity": "Upstream synthetic population VCF and population map",
        },
        "sha256": {name: hashlib.sha256(data).hexdigest() for name, data in payload.items()},
    }
    payload["package_provenance.json"] = encoded(provenance)
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(output, "w", zipfile.ZIP_DEFLATED) as archive:
        for name, data in sorted(payload.items()):
            item = zipfile.ZipInfo(name, date_time=(2026, 1, 1, 0, 0, 0))
            item.compress_type = zipfile.ZIP_DEFLATED
            item.external_attr = 0o100644 << 16
            archive.writestr(item, data)
    return output


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("dist/clawbio-research-0.1.0.zip"))
    args = parser.parse_args()
    print(build(Path(__file__).resolve().parents[1], args.output))
