"""Byte-stable nf-core params YAML; no scientific defaults."""

import math
import re
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

def validate_parameter_scalars(
    params: dict, contract: dict, *, patterns: bool = True
) -> None:
    """Check pinned scalar constraints; leave relationships to preflight.

    nf-schema ignores null/false unset values, including reference sentinels.
    Quoted booleans stay strings so they cannot silently become truthy flags.
    """
    for name, value in params.items():
        spec = contract.get(name)
        if not spec or value is None or value is False or value == "":
            continue
        kind = spec.get("type")
        allowed = [kind] if isinstance(kind, str) else (kind or [])
        numeric = isinstance(value, (int, float)) and not isinstance(value, bool)
        matches = {
            "string": isinstance(value, str),
            "boolean": isinstance(value, bool),
            "integer": numeric
            and (not isinstance(value, float) or math.isfinite(value))
            and int(value) == value,
            "number": numeric
            and (not isinstance(value, float) or math.isfinite(value)),
            "object": isinstance(value, dict),
            "array": isinstance(value, list),
        }
        if allowed and not any(matches.get(option, False) for option in allowed):
            raise ValueError(
                f"{name} must have type {kind}; got {type(value).__name__}. Numeric values must be finite."
            )
        if patterns:
            if "enum" in spec and value not in spec["enum"]:
                raise ValueError(f"{name} must be one of {spec['enum']}.")
            if isinstance(value, str):
                if "pattern" in spec and not re.search(spec["pattern"], value):
                    raise ValueError(
                        f"{name} does not match the upstream pattern {spec['pattern']!r}."
                    )
                if "minLength" in spec and len(value) < spec["minLength"]:
                    raise ValueError(
                        f"{name} is shorter than its upstream minimum length."
                    )
                if "maxLength" in spec and len(value) > spec["maxLength"]:
                    raise ValueError(f"{name} exceeds its upstream maximum length.")
        if numeric:
            if "minimum" in spec and value < spec["minimum"]:
                raise ValueError(f"{name} must be >= {spec['minimum']}.")
            if "maximum" in spec and value > spec["maximum"]:
                raise ValueError(f"{name} must be <= {spec['maximum']}.")
            branches = spec.get("oneOf", [])
            if branches and not any(
                ("const" in branch and value == branch["const"])
                or ("minimum" in branch and value >= branch["minimum"])
                for branch in branches
            ):
                raise ValueError(
                    f"{name} does not satisfy its numeric alternatives: {branches}."
                )

def validate_against_contract(params: dict, contract: dict, *, patterns: bool = True) -> None:
    """Validate against a pinned contract, skipping what upstream skips.

    nf-schema does not validate the pipeline's validation.defaultIgnoreParams:
    nf-core/sarek 3.8.1 lists freebayes_filter and vep_cache_version there because
    its schema types them as strings while its own defaults are integers.
    """
    ignored = set(contract.get("ignored_parameters", ()))
    checked = {name: value for name, value in params.items() if name not in ignored}
    validate_parameter_scalars(checked, contract["parameters"], patterns=patterns)


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
