#!/usr/bin/env python3
"""Execute one fixed, bundled research demo. Private inputs are unsupported."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
from importlib.metadata import version
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def verify_package() -> dict:
    provenance = json.loads((ROOT / "package_provenance.json").read_text())
    for name, expected in provenance["sha256"].items():
        path = ROOT / name
        if path.is_symlink() or not path.is_file() or not path.resolve().is_relative_to(ROOT):
            raise ValueError(f"PACKAGE_INTEGRITY_FAILED: {name}")
        if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            raise ValueError(f"PACKAGE_INTEGRITY_FAILED: {name}")
    return provenance


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workflow", choices=("pharmgx", "acmg", "equity"))
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    try:
        provenance = verify_package()
    except (ValueError, KeyError, OSError) as exc:
        print(str(exc), file=sys.stderr)
        return 2
    output = args.output.resolve()
    if output.exists() or output.is_relative_to(ROOT):
        print("OUTPUT_REFUSED: choose a fresh directory outside the plugin", file=sys.stderr)
        return 2
    runtime = ROOT / "runtime"
    workflows = {
        "pharmgx": ["skills/pharmgx-reporter/pharmgx_reporter.py", "--input",
                    str(runtime / "skills/pharmgx-reporter/demo_patient.txt"), "--no-enrich"],
        "acmg": ["skills/clinical-variant-reporter/clinical_variant_reporter.py", "--demo"],
        "equity": ["skills/equity-scorer/equity_scorer.py", "--input",
                   str(runtime / "examples/demo_populations.vcf"), "--pop-map",
                   str(runtime / "examples/demo_population_map.csv")],
    }
    command = workflows[args.workflow]
    command = [sys.executable, str(runtime / command[0]), *command[1:], "--output", str(output)]
    env = dict(os.environ)
    # Never enable optional outbound telemetry for these offline demos.
    env.pop("CLAWBIO_OTLP_ENDPOINT", None)
    try:
        result = subprocess.run(command, cwd=runtime, env=env, timeout=300)
    except subprocess.TimeoutExpired:
        print("DEMO_TIMEOUT: inspect partial output before retrying", file=sys.stderr)
        return 2
    if result.returncode:
        return result.returncode
    if not (output / "report.md").is_file() or not (output / "result.json").is_file():
        print("OUTPUT_CONTRACT_FAILED", file=sys.stderr)
        return 2
    # Preserve upstream findings while adding the package scope and canonical disclaimer.
    sys.path.insert(0, str(runtime))
    from clawbio.common.report import DISCLAIMER
    report = output / "report.md"
    with report.open("a") as stream:
        stream.write("\n\n## Plugin demonstration scope\n\n")
        stream.write(provenance["datasets"][args.workflow] + ".\n\n")
        stream.write("Bundled demonstration data only. Successful execution does not establish "
                     "independent clinical validation or an institutional equity certification.\n\n")
        stream.write(DISCLAIMER + "\n")
    receipt = {
        "schema_version": 1,
        "workflow": args.workflow,
        "mode": "bundled-demo",
        "dataset": provenance["datasets"][args.workflow],
        "private_input_supported": False,
        "package_verified": True,
        "package_provenance_sha256": hashlib.sha256((ROOT / "package_provenance.json").read_bytes()).hexdigest(),
        "upstream_base_commit": provenance["upstream_base_commit"],
        "completed_at": datetime.now(timezone.utc).isoformat(),
        "command": command,
        "scientific_validation": "Not established by package integrity or demo execution",
    }
    (output / "plugin_run.json").write_text(json.dumps(receipt, indent=2) + "\n")
    from clawbio.common.reproducibility import write_checksums, write_commands_sh, write_environment_yml
    # Retain upstream replay commands, and cover the final packaged report/receipt.
    existing_commands = output / "reproducibility/commands.sh"
    if existing_commands.exists():
        existing_commands.rename(output / "reproducibility/upstream_commands.sh")
    replay = shlex.join([sys.executable, str(ROOT / "scripts/run_research_demo.py"),
                        args.workflow, "--output"])
    write_commands_sh(output, 'set -eu\n' + replay + ' "${1:?Provide a fresh output directory}"')
    existing_environment = output / "reproducibility/environment.yml"
    if existing_environment.exists():
        existing_environment.rename(output / "reproducibility/upstream_environment.yml")
    distributions = ("numpy", "pandas", "matplotlib", "scikit-learn", "PyYAML", "opentelemetry-sdk")
    write_environment_yml(
        output, "clawbio-research-demo", [f"{name}=={version(name)}" for name in distributions],
        python_version=".".join(map(str, sys.version_info[:3])), channels=["conda-forge", "nodefaults"],
    )
    existing_checksums = output / "reproducibility/checksums.sha256"
    if existing_checksums.exists():
        existing_checksums.rename(output / "reproducibility/upstream_checksums.sha256")
    files = sorted(path for path in output.rglob("*") if path.is_file())
    write_checksums(files, output, anchor=output)
    print(f"Research demo complete: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
