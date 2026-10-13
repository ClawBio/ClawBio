"""Byte-stable nf-core params YAML; no scientific defaults."""
from pathlib import Path
import yaml
from clawbio.common.textio import write_text_lf


def serialize_params_yaml(params: dict[str, object], *, pipeline_name: str) -> str:
    header = (
        "# 'input' and 'outdir' are relative to the Nextflow launch directory (output_dir).\n"
        "# Use reproducibility/commands.sh for the recommended replay path.\n"
        f"# Manual replay: cd <output_dir> && nextflow run {pipeline_name} -params-file reproducibility/params.yaml\n"
    )
    return header + yaml.dump(params, allow_unicode=True, sort_keys=False)


def write_params_yaml(params: dict[str, object], *, output_dir: Path, pipeline_name: str) -> Path:
    repro_dir = output_dir / "reproducibility"
    repro_dir.mkdir(parents=True, exist_ok=True)
    params_path = repro_dir / "params.yaml"
    write_text_lf(params_path, serialize_params_yaml(params, pipeline_name=pipeline_name))
    return params_path
