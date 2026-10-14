#!/usr/bin/env python3
"""Offline donor/cell holdout audit of explicitly pinned observation metadata.

This does not fetch an atlas, inspect expression matrices, run a model or certify
unknown model-training data. PASS concerns the declared metadata only.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import re
import shlex
import sys
from pathlib import Path

VERSION = "0.1.0"
UNKNOWN = {"", "unknown", "na", "n/a", "nan", "none", "null", "unassigned", "not reported"}
FLOATING = {"latest", "stable", "current", "main", "master", "head"}
DISCLAIMER = (
    "ClawBio is a research and educational tool. It is not a medical device "
    "and does not provide clinical diagnoses. Consult a healthcare professional "
    "before making any medical decisions."
)
LIMITATIONS = [
    "PASS means no donor/cell overlap found in the declared metadata, not a clean model-training history.",
    "Identifiers must be stable and globally unique within the same study namespace across releases. Aliases are not resolved.",
    "Metadata checksums pin the inspected bytes. Source release labels and expected cell counts are caller declarations.",
    "Expression matrices, export completeness against H5AD, labels and upstream preprocessing leakage are not verified.",
    "A combined atlas release may contain earlier samples. A newer version is not evidence of independent holdout data.",
    "Donor metadata is descriptive. Missing ancestry is unknown, never evidence of representativeness or inferred ancestry.",
]


class InputError(ValueError):
    def __init__(self, code, message):
        super().__init__(message)
        self.code = code


def _pairs(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise InputError("INVALID_INPUT", "Duplicate JSON key is not permitted.")
        result[key] = value
    return result


def _reject_constant(value):
    raise InputError("INVALID_INPUT", "Non-finite JSON values are not permitted.")


def _json(data):
    return json.loads(data, object_pairs_hook=_pairs, parse_constant=_reject_constant)


def _text(value, label):
    if (not isinstance(value, str) or value != value.strip()
            or value.casefold() in UNKNOWN or any(ord(c) < 32 for c in value)):
        raise InputError("INVALID_INPUT", f"{label} must be a known, non-empty, whitespace-trimmed string.")
    return value


def _finding(findings, code, severity, message, split=None):
    item = {"code": code, "severity": severity, "message": message}
    if split:
        item["split"] = split
    findings.append(item)


def _validate_spec(spec):
    if not isinstance(spec, dict):
        raise InputError("INVALID_INPUT", "Each split must be an object.")
    for key in ("dataset_id", "dataset_version", "source_url", "identity_namespace",
                "observations_path", "observations_sha256"):
        _text(spec.get(key), key)
    if spec["dataset_version"].casefold() in FLOATING:
        raise InputError("INVALID_INPUT", "Pin a dataset release, not a floating version alias.")
    if not re.fullmatch(r"https?://[^\s/]+(?:/[^\s]*)?", spec["source_url"]):
        raise InputError("INVALID_INPUT", "source_url must be an HTTP(S) citation URL; it is never fetched.")
    if "://" in spec["observations_path"]:
        raise InputError("INVALID_INPUT", "observations_path must name a local JSONL file, not a URL.")
    if not re.fullmatch(r"[0-9a-fA-F]{64}", spec["observations_sha256"]):
        raise InputError("INVALID_INPUT", "observations_sha256 must be a SHA-256 hex digest.")
    if type(spec.get("expected_cells")) is not int or spec["expected_cells"] <= 0:
        raise InputError("INVALID_INPUT", "expected_cells must be a positive integer.")


def _demographics(donors, records, field, findings, split):
    distribution = {}
    known = conflict = 0
    for donor in donors:
        values = records.get(donor, set())
        if len(values) == 1:
            known += 1
            value = next(iter(values))
            distribution[value] = distribution.get(value, 0) + 1
        elif len(values) > 1:
            conflict += 1
    unknown = len(donors) - known - conflict
    if unknown:
        _finding(findings, f"{field.upper()}_COVERAGE_INCOMPLETE", "warning",
                 f"{unknown} donors have no reported {field}; no values were inferred.", split)
    if conflict:
        _finding(findings, f"{field.upper()}_METADATA_CONFLICT", "warning",
                 f"{conflict} donors have conflicting {field} values; excluded from the distribution.", split)
    return {"known_donors": known, "unknown_donors": unknown,
            "conflicting_donors": conflict, "distribution": dict(sorted(distribution.items()))}


def _read_split(name, spec, parent, findings):
    path = (parent / spec["observations_path"]).resolve()
    digest = hashlib.sha256()
    cells, donors = set(), set()
    demographics = {"sex": {}, "ancestry": {}}
    count = duplicates = 0
    # Hash exactly the bytes parsed in this single read, avoiding a separate hash/read race.
    with path.open("rb") as handle:
        for line_number, line in enumerate(handle, 1):
            digest.update(line)
            if not line.strip():
                continue
            try:
                row = _json(line)
            except (ValueError, UnicodeError) as error:
                raise InputError("INVALID_INPUT", f"Invalid JSONL in {name}, line {line_number}.") from error
            if not isinstance(row, dict):
                raise InputError("INVALID_INPUT", f"{name} line {line_number} must be an object.")
            cell = _text(row.get("cell_id"), f"{name} cell_id at line {line_number}")
            donor = _text(row.get("donor_id"), f"{name} donor_id at line {line_number}")
            duplicates += int(cell in cells)
            cells.add(cell)
            donors.add(donor)
            count += 1
            for field, groups in demographics.items():
                value = row.get(field)
                if value is None:
                    continue
                if not isinstance(value, str):
                    raise InputError("INVALID_INPUT", f"Optional {field} must be a string or null.")
                value = value.strip()
                if value.casefold() not in UNKNOWN:
                    groups.setdefault(donor, set()).add(value)
    actual_hash = digest.hexdigest()
    if actual_hash != spec["observations_sha256"].lower():
        raise InputError("CHECKSUM_MISMATCH", f"{name} observations do not match the declared SHA-256.")
    if count != spec["expected_cells"]:
        raise InputError("CELL_COUNT_MISMATCH", f"{name}: read {count} cells; expected {spec['expected_cells']}.")
    if duplicates:
        _finding(findings, "DUPLICATE_CELL_ID", "failure",
                 f"{duplicates} repeated cell IDs within the split.", name)
    summary = {
        "dataset_id": spec["dataset_id"], "dataset_version": spec["dataset_version"],
        "source_url": spec["source_url"], "identity_namespace": spec["identity_namespace"],
        "observations_path": str(path), "observations_sha256": actual_hash,
        "cells": count, "unique_cells": len(cells), "donors": len(donors),
        "duplicate_cell_rows": duplicates,
        "donor_metadata": {field: _demographics(donors, groups, field, findings, name)
                           for field, groups in demographics.items()},
    }
    return summary, cells, donors


def _audit(path, raw):
    result = {
        "tool": "atlas-benchmark-preflight", "tool_version": VERSION,
        "status": "ERROR", "scope": "declared_observation_metadata_only",
        "model_training_contamination": "not_assessed", "synthetic_demo": False,
        "splits": {}, "overlap": None, "findings": [], "limitations": LIMITATIONS.copy(),
        "provenance": {"manifest_path": str(path),
                       "manifest_sha256": hashlib.sha256(raw).hexdigest()},
    }
    try:
        manifest = _json(raw)
        if not isinstance(manifest, dict) or type(manifest.get("schema_version")) is not int or manifest["schema_version"] != 1:
            raise InputError("INVALID_INPUT", "Expected manifest schema_version 1.")
        splits = manifest.get("splits")
        if not isinstance(splits, dict) or set(splits) != {"train", "test"}:
            raise InputError("INVALID_INPUT", "Exactly train and test splits are required.")
        for spec in splits.values():
            _validate_spec(spec)
        result["synthetic_demo"] = manifest.get("synthetic_demo") is True
        if splits["train"]["identity_namespace"] != splits["test"]["identity_namespace"]:
            raise InputError("IDENTITY_NAMESPACE_MISMATCH", "Cannot compare IDs from different identity namespaces.")
        identities = {}
        for name in ("train", "test"):
            summary, cells, donors = _read_split(name, splits[name], path.parent, result["findings"])
            result["splits"][name] = summary
            identities[name] = (cells, donors)
        cell_overlap = len(identities["train"][0] & identities["test"][0])
        donor_overlap = len(identities["train"][1] & identities["test"][1])
        result["overlap"] = {"cells": cell_overlap, "donors": donor_overlap}
        for count, code in ((cell_overlap, "CELL_OVERLAP"), (donor_overlap, "DONOR_OVERLAP")):
            if count:
                _finding(result["findings"], code, "failure",
                         f"{count} shared {code.split('_')[0].lower()} identities between train and test.")
        result["status"] = "FAIL" if any(f["severity"] == "failure" for f in result["findings"]) else "PASS"
    except (InputError, ValueError, OSError, UnicodeError) as error:
        _finding(result["findings"], getattr(error, "code", "INVALID_INPUT"), "error", str(error))
    return result


def audit_manifest(path):
    """Audit pinned local JSONL observation exports; never use network or write files."""
    path = Path(path).resolve()
    try:
        raw = path.read_bytes()
    except OSError as error:
        result = _audit(path, b"")
        result["provenance"]["manifest_sha256"] = None
        result["findings"] = [{"code": "INVALID_INPUT", "severity": "error", "message": str(error)}]
        return result
    return _audit(path, raw)


def _write_bundle(out, result, raw, argv):
    # Refuse existing directories, files and symlinks, including concurrent creation.
    out.parent.mkdir(parents=True, exist_ok=True)
    out.mkdir(exist_ok=False)
    repro = out / "reproducibility"
    repro.mkdir()
    (out / "result.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    lines = ["# Atlas benchmark preflight", "", f"Status: **{result['status']}**", "",
             "Scope: declared observation metadata only. Model-training contamination was not assessed.", ""]
    if result["synthetic_demo"]:
        lines += ["**SYNTHETIC DEMO: these are toy records, not Tabula Sapiens observations.**", ""]
    lines += ["| Split | Cells | Unique donors |", "|---|---:|---:|"]
    for name, split in result["splits"].items():
        lines.append(f"| {name} | {split['cells']} | {split['donors']} |")
    if result["overlap"] is not None:
        lines += ["", f"Shared cells: {result['overlap']['cells']}. Shared donors: {result['overlap']['donors']}."]
    lines += ["", "## Findings", ""]
    lines.extend(f"- {f['severity'].upper()}: {f['code']}" + (f" ({f['split']})" if "split" in f else "")
                 for f in result["findings"])
    if not result["findings"]:
        lines.append("No overlap or metadata warnings found within the declared scope.")
    lines += ["", "See result.json for structured findings and donor-level metadata coverage.", "",
              "## Interpretation boundaries", ""] + [f"- {line}" for line in LIMITATIONS]
    lines += ["", DISCLAIMER, ""]
    (out / "report.md").write_text("\n".join(lines), encoding="utf-8")
    if raw is not None:
        (repro / "manifest.json").write_bytes(raw)
    script = Path(__file__).resolve()
    (repro / "environment.json").write_text(json.dumps({
        "python": platform.python_version(), "platform": platform.platform(),
        "tool_version": VERSION, "script_sha256": hashlib.sha256(script.read_bytes()).hexdigest(),
        "dependencies": "Python standard library only", "network_access": "none",
        "working_directory": str(Path.cwd()),
    }, indent=2) + "\n", encoding="utf-8")
    (repro / "commands.sh").write_text(
        "#!/usr/bin/env bash\n# Original invocation. Choose a fresh --output directory for a rerun.\n"
        + "cd " + shlex.quote(str(Path.cwd())) + "\n"
        + shlex.join([sys.executable, str(script), *argv]) + "\n", encoding="utf-8")
    hashes = []
    for path in sorted(out.rglob("*")):
        if path.is_file():
            hashes.append(f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.relative_to(out).as_posix()}")
    (repro / "checksums.sha256").write_text("\n".join(hashes) + "\n", encoding="utf-8")


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    parser = argparse.ArgumentParser(description=__doc__)
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument("--input", type=Path, help="Pinned manifest for local JSONL metadata exports.")
    inputs.add_argument("--demo", action="store_true", help="Synthetic leaking split, expected FAIL and exit 1.")
    parser.add_argument("--output", type=Path, required=True, help="New output directory; existing paths are refused.")
    args = parser.parse_args(argv)
    path = (Path(__file__).parent / "fixtures" / "atlas_preflight" / "demo_manifest.json"
            if args.demo else args.input).resolve()
    try:
        raw = path.read_bytes()
        result = _audit(path, raw)
    except OSError:
        raw = None
        result = audit_manifest(path)
    try:
        _write_bundle(args.output.absolute(), result, raw, argv)
    except OSError as error:
        print(f"Refusing to overwrite, or unable to write, output: {error}", file=sys.stderr)
        return 2
    print(json.dumps({"status": result["status"], "scope": result["scope"],
                      "output": str(args.output.absolute()), "overlap": result["overlap"]}))
    return {"PASS": 0, "FAIL": 1, "ERROR": 2}[result["status"]]


if __name__ == "__main__":
    raise SystemExit(main())
