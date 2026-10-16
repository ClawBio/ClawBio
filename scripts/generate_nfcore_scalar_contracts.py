#!/usr/bin/env python3
"""Generate scalar type/bounds contracts from locally downloaded pinned sources.

Each contract also records the parameters the pinned nextflow.config lists in
validation.defaultIgnoreParams: nf-schema does not validate them, so neither do
the wrappers.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PIPELINES = ("rnaseq", "scrnaseq", "sarek")
FIELDS = {
    "type",
    "minimum",
    "maximum",
    "oneOf",
    "enum",
    "pattern",
    "minLength",
    "maxLength",
}


def pinned_revision(skill_dir: Path) -> str:
    tree = ast.parse((skill_dir / "schemas.py").read_text(encoding="utf-8"))
    for statement in tree.body:
        if isinstance(statement, ast.Assign) and any(
            isinstance(target, ast.Name) and target.id == "DEFAULT_PIPELINE_VERSION"
            for target in statement.targets
        ):
            return ast.literal_eval(statement.value)
    raise ValueError(f"No pinned revision in {skill_dir}")


_IGNORE_PARAMS = re.compile(r"defaultIgnoreParams\s*=\s*\[(?P<items>[^\]]*)\]")


def ignored_parameters(config_text: str) -> list[str]:
    """validation.defaultIgnoreParams from nextflow.config, in upstream order."""
    found = _IGNORE_PARAMS.search(config_text)
    return re.findall(r"[\"']([^\"']+)[\"']", found.group("items")) if found else []


def generate_contract(name: str, revision: str, source_dir: Path) -> str:
    raw = (source_dir / "nextflow_schema.json").read_bytes()
    config = (source_dir / "nextflow.config").read_bytes()
    schema = json.loads(raw)
    parameters = {
        key: {field: value for field, value in spec.items() if field in FIELDS}
        for group in schema.get("$defs", schema.get("definitions", {})).values()
        for key, spec in group.get("properties", {}).items()
    }
    if not parameters:
        raise ValueError(f"No parameter definitions in {source_dir}")
    metadata = {
        "pipeline": f"nf-core/{name}",
        "revision": revision,
        "schema_url": f"https://raw.githubusercontent.com/nf-core/{name}/{revision}/nextflow_schema.json",
        "schema_sha256": hashlib.sha256(raw).hexdigest(),
        "config_url": f"https://raw.githubusercontent.com/nf-core/{name}/{revision}/nextflow.config",
        "config_sha256": hashlib.sha256(config).hexdigest(),
        "ignored_parameters": ignored_parameters(config.decode("utf-8")),
    }
    entries = [
        json.dumps(key) + ": " + json.dumps(value, sort_keys=True)
        for key, value in sorted(parameters.items())
    ]
    return (
        json.dumps(metadata, indent=2)[:-2]
        + ',\n  "parameters": {\n    '
        + ",\n    ".join(entries)
        + "\n  }\n}\n"
    )


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--upstream-dir",
        required=True,
        type=Path,
        help="Contains rnaseq-<tag>, scrnaseq-<tag>, and sarek-<tag> source trees",
    )
    parser.add_argument(
        "--check", action="store_true", help="Report drift without changing contracts"
    )
    args = parser.parse_args(argv)
    failed = False
    for name in PIPELINES:
        directory = ROOT / "skills" / f"nfcore-{name}-wrapper"
        revision = pinned_revision(directory)
        expected = generate_contract(
            name,
            revision,
            args.upstream_dir / f"{name}-{revision}",
        )
        target = directory / "reproducibility/upstream_scalar_contract.json"
        if args.check:
            matches = (
                target.is_file() and target.read_text(encoding="utf-8") == expected
            )
            failed |= not matches
            print(f"{name} {revision}: {'matches' if matches else 'DRIFT'}")
        else:
            target.write_text(expected, encoding="utf-8", newline="\n")
            print(f"{name} {revision}: generated")
    return int(failed)


if __name__ == "__main__":
    raise SystemExit(main())
