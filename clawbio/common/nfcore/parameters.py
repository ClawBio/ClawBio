"""Byte-stable nf-core params YAML; no scientific defaults."""

import math
from pathlib import Path

import yaml

from clawbio.common.textio import write_text_lf

# PyYAML's default emitter width, which wrapped long scalars onto continuation
# lines before params.yaml was written unwrapped. Kept to recognise old bundles.
LEGACY_LINE_WIDTH = 80


def serialize_params_yaml(
    params: dict[str, object], *, pipeline_name: str, line_width: float = math.inf
) -> str:
    """params.yaml text with one value per line.

    Unwrapped so the stdlib-only remap_paths.py, which edits params.yaml line by
    line, always sees each reference path whole.
    """
    header = (
        "# 'input' and 'outdir' are relative to the Nextflow launch directory (output_dir).\n"
        "# Use reproducibility/commands.sh for the recommended replay path.\n"
        f"# Manual replay: cd <output_dir> && nextflow run {pipeline_name} -params-file reproducibility/params.yaml\n"
    )
    return header + yaml.dump(
        params, allow_unicode=True, sort_keys=False, width=line_width
    )


def write_params_yaml(
    params: dict[str, object], *, output_dir: Path, pipeline_name: str
) -> Path:
    repro_dir = output_dir / "reproducibility"
    repro_dir.mkdir(parents=True, exist_ok=True)
    params_path = repro_dir / "params.yaml"
    write_text_lf(
        params_path, serialize_params_yaml(params, pipeline_name=pipeline_name)
    )
    return params_path
