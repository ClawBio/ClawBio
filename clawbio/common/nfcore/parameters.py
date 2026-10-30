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


def load_unique_yaml(text: str):
    """Load safe YAML with explicit duplicate mapping keys rejected."""
    import yaml

    class UniqueLoader(yaml.SafeLoader):
        def construct_mapping(self, node, deep=False):
            seen = set()
            for key_node, _ in node.value:
                if key_node.tag == "tag:yaml.org,2002:merge":
                    continue
                key = self.construct_object(key_node, deep=deep)
                if not isinstance(key, str) or key in seen:
                    raise yaml.constructor.ConstructorError(
                        "while reading a mapping",
                        node.start_mark,
                        f"Expected a unique string key, got {key!r}",
                        key_node.start_mark,
                    )
                seen.add(key)
            return super().construct_mapping(node, deep=deep)

    return yaml.load(text, Loader=UniqueLoader)

def reject_non_finite_numbers(params: dict) -> None:
    """Refuse NaN and infinity anywhere in the parameters.

    Every other rule belongs to nf-schema, which checks the pinned
    nextflow_schema.json when Nextflow starts. It cannot judge these: it serialises
    the parameters to JSON first, which has no NaN or infinity, so it fails with a
    serialisation error that names no parameter.
    """

    def non_finite(value) -> bool:
        if isinstance(value, float):
            return not math.isfinite(value)
        if isinstance(value, dict):
            return any(non_finite(item) for item in value.values())
        if isinstance(value, (list, tuple)):
            return any(non_finite(item) for item in value)
        return False

    for name, value in params.items():
        if non_finite(value):
            raise ValueError(f"{name!r} holds NaN or infinity; use a finite number.")


def load_unique_json(text: str):
    """Reject repeated JSON keys instead of silently taking the final value."""
    import json

    def mapping(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise ValueError(f"Duplicate JSON mapping key: {key!r}")
            result[key] = value
        return result

    return json.loads(text, object_pairs_hook=mapping)
