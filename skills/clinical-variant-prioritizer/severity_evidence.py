"""Offline, fail-closed validation of supplied disease-severity evidence.

Aggregation follows T0hid/hpo-classification-agent at UPSTREAM_COMMIT.
Evidence completeness and citation presence are checked, not biological truth.
No patient information or evidence is sent to external services.
"""
from __future__ import annotations

from collections import Counter
from datetime import date
from fractions import Fraction
import hashlib
import json
import re

UPSTREAM_COMMIT = "03bc5f6be6456a3ca2c5206e399f0a38f879ce57"
UNKNOWN = "severity_unknown"
LABELS = {"Profound", "Severe", "Moderate", "Mild"}
FLOATING = {"latest", "main", "master", "head", "unknown", "unversioned"}


def aggregate(counts: dict) -> str:
    """Apply published tier-count rules to already validated, unique HPO terms."""
    first, second, third = (counts.get(str(i), 0) for i in (1, 2, 3))
    if first > 1:
        return "Profound"
    if first == 1:
        return "Severe"
    if second:
        return "Severe" if second + third >= 4 else "Moderate"
    return "Moderate" if third else "Mild"


def _text(value) -> bool:
    return isinstance(value, str) and bool(value.strip())


def _release(value) -> bool:
    if not isinstance(value, str) or not re.fullmatch(r"\d{4}-\d{2}-\d{2}", value):
        return False
    try:
        date.fromisoformat(value)
        return True
    except ValueError:
        return False


def _frequency(value) -> bool:
    """Return eligibility. Reject missing, nonfinite and ambiguous encodings."""
    if isinstance(value, bool) or value is None:
        raise ValueError("invalid_frequency")
    if isinstance(value, str):
        value = value.strip()
        if value in {"HP:0040280", "HP:0040281", "HP:0040282"}:
            return True
        if value in {"HP:0040283", "HP:0040284", "HP:0040285"}:
            return False
    try:
        if isinstance(value, str) and value.endswith("%"):
            number = Fraction(value[:-1]) / 100
        else:
            number = Fraction(str(value))
    except (ValueError, ZeroDivisionError, TypeError):
        raise ValueError("invalid_frequency") from None
    if not 0 <= number <= 1:
        raise ValueError("invalid_frequency")
    return number >= Fraction(3, 10)


def evaluate(bundle, gene, disease_id, declared_severity=None) -> dict:
    """Return a separate research annotation; never revise pathogenicity.

    Input is a JSON-shaped object. Any invalid record invalidates the bundle:
    silently discarding malformed evidence could bias the severity downwards.
    The caller must supply a curated gene + MONDO identity; no fuzzy joins.
    """
    result = {"severity": UNKNOWN, "warnings": [], "reasons": [], "evidence": [],
              "hpo_release": None, "upstream_commit": UPSTREAM_COMMIT,
              "evidence_sha256": None, "tier_counts": {}, "excluded_count": 0,
              "validation_scope": "structure_and_provenance_presence_only"}
    try:
        serialised = json.dumps(bundle, sort_keys=True, separators=(",", ":"), allow_nan=False)
        result["evidence_sha256"] = hashlib.sha256(serialised.encode()).hexdigest()
    except (TypeError, ValueError):
        result["reasons"] = ["not_finite_json"]
        return result
    try:
        if not isinstance(bundle, dict) or bundle.get("schema_version") != "1.0":
            raise ValueError("unsupported_schema")
        if bundle.get("upstream_commit") != UPSTREAM_COMMIT:
            raise ValueError("incompatible_upstream_commit")
        if not _release(bundle.get("hpo_release")):
            raise ValueError("missing_or_invalid_hpo_release")
        result["hpo_release"] = bundle["hpo_release"]
        records = bundle.get("records")
        if not isinstance(records, list) or not records:
            raise ValueError("missing_records")
        checked, seen = [], set()
        for row in records:
            if not isinstance(row, dict):
                raise ValueError("invalid_record")
            for field in ("gene", "disease_id", "disease", "hpo_id", "citation",
                          "source_span", "source_version", "hpo_release"):
                if not _text(row.get(field)):
                    raise ValueError("missing_" + field)
            if not re.fullmatch(r"MONDO:\d{7}", row["disease_id"]):
                raise ValueError("invalid_disease_id")
            if not re.fullmatch(r"HP:\d{7}", row["hpo_id"]):
                raise ValueError("invalid_hpo_id")
            if row["hpo_release"] != bundle["hpo_release"]:
                raise ValueError("incompatible_hpo_release")
            if row["source_version"].strip().lower() in FLOATING:
                raise ValueError("unversioned_source")
            if not re.fullmatch(r"(?:PMID:[1-9]\d*|https?://[^\s]+|doi:10\.\d{4,9}/[^\s]+)", row["citation"]):
                raise ValueError("invalid_citation")
            if type(row.get("tier")) is not int or row["tier"] not in (1, 2, 3, 4):
                raise ValueError("invalid_tier")
            qualifier = row.get("qualifier")
            if not isinstance(qualifier, str) or qualifier.strip().upper() not in ("", "NOT"):
                raise ValueError("invalid_qualifier")
            eligible = _frequency(row.get("frequency"))
            key = (row["gene"], row["disease_id"], row["hpo_id"])
            if key in seen:
                raise ValueError("duplicate_hpo_evidence")
            seen.add(key)
            checked.append((row, eligible and qualifier.strip().upper() != "NOT"))
        if not _text(gene) or not isinstance(disease_id, str) or not re.fullmatch(r"MONDO:\d{7}", disease_id):
            raise ValueError("missing_gene_disease_identity")
        matched = [(r, keep) for r, keep in checked
                   if r["gene"] == gene and r["disease_id"] == disease_id]
        if not matched:
            raise ValueError("no_matching_evidence")
        if len({r["disease"] for r, _ in matched}) != 1:
            raise ValueError("inconsistent_disease_labels")
        kept = [r for r, keep in matched if keep]
        result["excluded_count"] = len(matched) - len(kept)
        # Retain supplied source spans for included and excluded records.
        result["evidence"] = [dict(r, included=keep) for r, keep in matched]
        if not kept:
            raise ValueError("no_eligible_evidence")
        counts = Counter(str(r["tier"]) for r in kept)
        result["tier_counts"] = {str(i): counts[str(i)] for i in range(1, 5)}
        result["severity"] = aggregate(counts)
        if _text(declared_severity):
            declared = declared_severity.strip().capitalize()
            if declared in LABELS and declared != result["severity"]:
                result["warnings"].append("severity_disagreement")
            elif declared not in LABELS:
                result["warnings"].append("legacy_severity_not_comparable")
    except ValueError as exc:
        result["reasons"].append(str(exc))
    return result


def annotate(findings: list[dict], panel: list[dict], bundle) -> dict:
    """Attach evidence by the finding's exact panel identifier, gene and MONDO."""
    entries = {entry["id"]: entry for entry in panel}
    identifiers = Counter(entry["id"] for entry in panel)
    for finding in findings:
        entry = entries[finding["id"]]
        if identifiers[finding["id"]] != 1:
            unknown = evaluate(bundle, None, None)
            unknown["reasons"].append("ambiguous_panel_identity")
            finding["disease_severity"] = unknown
            continue
        finding["disease_severity"] = evaluate(
            bundle, entry.get("gene"), entry.get("disease_id"), entry.get("condition_severity"))
    return {"supplied": True, "upstream_commit": UPSTREAM_COMMIT,
            "scope": "research_annotation_only",
            "classified": sum(f["disease_severity"]["severity"] != UNKNOWN for f in findings),
            "unknown": sum(f["disease_severity"]["severity"] == UNKNOWN for f in findings)}
