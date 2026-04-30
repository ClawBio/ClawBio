"""
repro_bundle.py — Creates reproducibility artefacts for NutriGx Advisor
Outputs: commands.sh, environment.yml, checksums.txt
"""

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


CONDA_ENV = """name: nutrigx-advisor
channels:
  - conda-forge
  - defaults
dependencies:
  - python=3.11
  - numpy>=1.26
  - pandas>=2.2
  - matplotlib>=3.8
  - seaborn>=0.13
  - pip
  - pip:
    - clawbio==0.1.0
"""


def sha256_file(filepath: str) -> str:
    h = hashlib.sha256()
    try:
        with open(filepath, "rb") as f:
            for chunk in iter(lambda: f.read(8192), b""):
                h.update(chunk)
        return h.hexdigest()
    except FileNotFoundError:
        return "FILE_NOT_FOUND"


def create_reproducibility_bundle(input_file: str, output_dir: str, panel_path: str, args: dict):
    output_dir = Path(output_dir)
    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    # commands.sh
    cmd_args = " ".join(f"--{k.replace('_', '-')} {v}" for k, v in args.items() if v and k != "synthetic")
    commands = f"""#!/usr/bin/env bash
# NutriGx Advisor — Reproducibility Script
# Generated: {timestamp}
# ClawBio NutriGx Advisor v0.1.0

set -euo pipefail

# 1. Create conda environment
conda env create -f environment.yml
conda activate nutrigx-advisor

# 2. Run analysis
python nutrigx_advisor.py {cmd_args}

# 3. Verify checksums
sha256sum -c checksums.txt
"""
    (output_dir / "commands.sh").write_text(commands)

    # environment.yml
    (output_dir / "environment.yml").write_text(CONDA_ENV)

    # checksums.txt
    files_to_checksum = [
        input_file,
        panel_path,
        str(output_dir / "nutrigx_report.md"),
    ]
    checksum_lines = [f"# NutriGx Advisor checksums — {timestamp}"]
    for fp in files_to_checksum:
        chk = sha256_file(fp)
        checksum_lines.append(f"{chk}  {Path(fp).name}")

    (output_dir / "checksums.txt").write_text("\n".join(checksum_lines) + "\n")

    # provenance.json
    provenance = {
        "tool": "ClawBio NutriGx Advisor",
        "version": "0.1.0",
        "timestamp": timestamp,
        "input_file": Path(input_file).name,
        "args": args,
    }
    (output_dir / "provenance.json").write_text(json.dumps(provenance, indent=2))
