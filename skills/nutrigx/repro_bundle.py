"""
repro_bundle.py — Creates reproducibility artefacts for NutriGx Advisor

Delegates to the shared clawbio.common reproducibility layer.
Outputs (in <output_dir>/reproducibility/): commands.sh, environment.yml,
checksums.sha256, provenance.json
"""

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

_PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from clawbio.common.checksums import sha256_file  # noqa: E402
from clawbio.common.reproducibility import (  # noqa: E402
    write_checksums,
    write_commands_sh,
    write_environment_yml,
)

VERSION = "0.1.0"


def create_reproducibility_bundle(input_file: str, output_dir: str, panel_path: str, args: dict):
    output_dir = Path(output_dir)
    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    cmd_args = " ".join(
        f"--{k.replace('_', '-')} {v}" for k, v in args.items() if v and k != "synthetic"
    )
    write_commands_sh(
        output_dir,
        f"# NutriGx Advisor — Reproducibility Script\n"
        f"# Generated: {timestamp}\n"
        f"# ClawBio NutriGx Advisor v{VERSION}\n"
        f"set -euo pipefail\n"
        f"\n"
        f"# 1. Create conda environment\n"
        f"conda env create -f environment.yml\n"
        f"conda activate nutrigx\n"
        f"\n"
        f"# 2. Run analysis\n"
        f"python nutrigx.py {cmd_args}\n"
        f"\n"
        f"# 3. Verify checksums\n"
        f"sha256sum -c checksums.sha256",
    )

    write_environment_yml(
        output_dir,
        env_name="nutrigx",
        python_version="3.11",
        conda_deps=["numpy>=1.26", "pandas>=2.2", "matplotlib>=3.8", "seaborn>=0.13"],
        pip_deps=["clawbio==0.1.0"],
    )

    write_checksums(
        [input_file, panel_path, output_dir / "nutrigx_report.md"],
        output_dir,
    )

    provenance = {
        "tool": "ClawBio NutriGx Advisor",
        "version": VERSION,
        "timestamp": timestamp,
        "input_file": Path(input_file).name,
        "args": args,
    }
    (output_dir / "reproducibility" / "provenance.json").write_text(
        json.dumps(provenance, indent=2)
    )
