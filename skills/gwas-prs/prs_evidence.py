"""Deterministic research-percentile gate for exact scoring files.

Evidence is caller-curated and hash-bound, not independently verified scientific
truth. No network, ancestry inference, clinical classification or causal inference.
"""
from __future__ import annotations

import hashlib
import math
from pathlib import Path

POLICY_VERSION = "1.0.0"


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _text(value):
    return isinstance(value, str) and bool(value.strip())


def _number(value):
    return type(value) in (int, float) and math.isfinite(value)


def _positive_int(value):
    return type(value) is int and value > 0


def _mapping(value):
    return value if isinstance(value, dict) else {}


def assess(*, score_id, scoring_path, input_path, build, variants, genotypes,
           raw_score, evidence=None, curated=False):
    """Return a stable decision; absent/incomplete evidence withholds interpretation.

    Supported means the declared evidence satisfies this policy for a *research
    percentile*. It does not certify external evidence or individual disease risk.
    """
    reasons = []
    def flag(code, action, incompatible=False):
        reasons.append({"code": code, "action": action, "incompatible": incompatible})

    result = {
        "policy_version": POLICY_VERSION, "status": "not_established",
        "scope": "research_percentile", "evidence_basis": "caller_supplied",
        "score_id": score_id, "scoring_sha256": digest(scoring_path),
        "input_sha256": digest(input_path), "reasons": reasons,
        "percentile": None, "z_score": None, "risk_category": None,
        "method": "withheld_evidence", "reference_population": None,
        "predictive_performance": None, "heritability": None,
        "causal_evidence": None, "causal_claim_supported": False,
    }
    if curated:
        flag("ILLUSTRATIVE_PANEL", "Use an exact published score with validated weights.")
    if evidence is None:
        flag("EVIDENCE_MISSING", "Supply a hash-bound evidence JSON; ancestry resemblance is insufficient.")
        return result
    if (not isinstance(evidence, dict) or type(evidence.get("schema_version")) is not int
            or evidence["schema_version"] != 1
            or not isinstance(evidence.get("input"), dict)
            or not isinstance(evidence.get("scores"), dict)):
        flag("EVIDENCE_SCHEMA_INVALID", "Use evidence schema version 1 with input and scores objects.")
        return result
    inp = evidence["input"]
    card = _mapping(evidence["scores"].get(score_id))
    if not card:
        flag("SCORE_EVIDENCE_MISSING", "Provide evidence for this exact score identifier.")
        return result
    result["score_version"] = card.get("version")
    result["score_source"] = card.get("source")
    if card.get("sha256") != result["scoring_sha256"]:
        flag("SCORE_IDENTITY_MISMATCH", "Verify the scoring-file hash against the cited score version.", True)
    if inp.get("sha256") != result["input_sha256"]:
        flag("INPUT_IDENTITY_MISMATCH", "Provide an input manifest for the exact genotype file.", True)
    if not all(_text(card.get(k)) for k in ("version", "source")):
        flag("SCORE_PROVENANCE_INCOMPLETE", "Record the score version and source.")
    if not all(_text(inp.get(k)) for k in ("population", "context", "source")):
        flag("TARGET_CONTEXT_MISSING", "Document target population, validation context and input source.")
    if inp.get("build") not in ("GRCh37", "GRCh38") or card.get("build") not in ("GRCh37", "GRCh38"):
        flag("BUILD_UNKNOWN", "Document the input and score genome builds.")
    elif not inp["build"] == card["build"] == build:
        flag("BUILD_MISMATCH", "Harmonise builds using a validated upstream workflow.", True)

    if inp.get("strand") != "forward":
        flag("STRAND_UNKNOWN", "Supply a source-backed forward-strand input declaration.")

    ids = [v.get("rsid") for v in variants]
    if len(ids) != len(set(ids)):
        flag("DUPLICATE_VARIANTS", "Resolve duplicate scoring variants upstream.", True)
    if not _positive_int(card.get("variant_count")) or card["variant_count"] != len(variants):
        flag("VARIANT_COUNT_MISMATCH", "Reconcile parsed variants with the exact published score.", True)
    if not variants or any(v.get("rsid") not in genotypes for v in variants):
        flag("INCOMPLETE_VARIANTS", "This policy requires the complete score; partial sums cannot use its reference.")
    bad_alleles = False
    for v in variants:
        ea, oa = v.get("effect_allele"), v.get("other_allele")
        gt = genotypes.get(v.get("rsid"))
        # v1 deliberately supports diploid biallelic SNVs only, with explicit
        # other alleles. It does not guess strand flips or impute missing calls.
        if (ea not in ("A", "C", "G", "T") or oa not in ("A", "C", "G", "T") or ea == oa
                or (gt is not None and (not isinstance(gt, str) or len(gt) != 2
                                       or not set(gt).issubset({ea, oa})))):
            bad_alleles = True
    if bad_alleles:
        flag("ALLELE_INCOMPATIBLE", "Supply harmonised diploid SNVs with explicit effect/other alleles and valid calls.", True)
    if not _number(raw_score) or any(not _number(v.get("effect_weight")) for v in variants):
        flag("NONFINITE_SCORE", "Resolve non-finite weights or raw score.", True)

    if (not bad_alleles and all(v.get("rsid") in genotypes for v in variants)
            and _number(raw_score) and all(_number(v.get("effect_weight")) for v in variants)):
        expected = sum(genotypes[v["rsid"]].count(v["effect_allele"]) * v["effect_weight"]
                       for v in variants)
        if not math.isclose(raw_score, expected, rel_tol=1e-12, abs_tol=1e-12):
            flag("RAW_SCORE_MISMATCH", "Recalculate the raw sum from these exact compatible dosages.", True)

    validation = _mapping(card.get("validation"))
    ref = _mapping(card.get("reference"))
    result["predictive_performance"] = validation or None
    # These fields are deliberately never substituted for prediction evidence.
    result["heritability"] = card.get("heritability")
    result["causal_evidence"] = card.get("causal_evidence")
    metric = validation.get("metric")
    lo, est, hi = (validation.get(k) for k in ("ci_lower", "estimate", "ci_upper"))
    valid_metric = (metric in ("r2", "auc", "correlation")
                    and all(_number(x) for x in (lo, est, hi)))
    if valid_metric:
        lower_bound = -1 if metric == "correlation" else 0
        valid_metric = lower_bound <= lo <= est <= hi <= 1
    if (validation.get("independent") is not True or not _positive_int(validation.get("n"))
            or not _text(validation.get("source")) or not valid_metric):
        flag("VALIDATION_INCOMPLETE", "Provide independent validation, sample size, source, metric and ordered confidence interval.")
    if (not all(_text(validation.get(k)) for k in ("population", "context"))
            or any(validation.get(k) != inp.get(k) for k in ("population", "context"))):
        flag("VALIDATION_CONTEXT_MISMATCH", "Obtain validation for the specified target population and context.")
    if (not _text(ref.get("source")) or not _positive_int(ref.get("n"))
            or not _number(ref.get("mean")) or not _number(ref.get("sd")) or ref.get("sd", 0) <= 0
            or ref.get("distribution") != "normal" or ref.get("missingness") != "complete_only"):
        flag("REFERENCE_INCOMPLETE", "Provide a documented normal reference distribution for the complete score.")
    if (not all(_text(ref.get(k)) for k in ("population", "context"))
            or any(ref.get(k) != inp.get(k) for k in ("population", "context"))):
        flag("REFERENCE_CONTEXT_MISMATCH", "Provide a reference distribution for this population and context.")

    if reasons:
        result["status"] = "incompatible" if any(r["incompatible"] for r in reasons) else "not_established"
    else:
        z = (raw_score - ref["mean"]) / ref["sd"]
        if not math.isfinite(z):
            flag("REFERENCE_NUMERIC_OVERFLOW", "Check raw score and reference scale.", True)
            result["status"] = "incompatible"
            return result
        result.update(status="supported", z_score=z,
                      percentile=50 * (1 + math.erf(z / math.sqrt(2))),
                      method="evidence_normal_reference", reference_population=ref["population"])
    return result
