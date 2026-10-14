#!/usr/bin/env python3
"""prs-applicability-gate: may this already-computed PGS result be interpreted for this person?

Input: ONE candidate score plus normalised evidence about it and the person (schema
``prs-applicability-gate.input.v2``, documented in SKILL.md). Output: SUPPORTED / RAW_ONLY / ABSTAIN with reason
codes, a full rule trace, the evidence used and missing, and what would change the result.

The gate is deliberately narrow: it does not search the PGS Catalog, rank scores, infer ancestry, compute scores
or choose a primary score. It is a pure, deterministic function of (input, calibration config).

Originally developed in PRSGuard (ClawBio Boston Hackathon, Challenge 3):
https://github.com/rinatrizvanov/prsguard-demo

    python prs_applicability_gate.py --input gate_input.json --output out_dir
    python prs_applicability_gate.py --demo --output /tmp/prs_gate_demo
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import shlex
import sys
from collections.abc import Callable
from pathlib import Path
from typing import Any

import yaml

SKILL_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SKILL_DIR.parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from clawbio.common.report import DISCLAIMER as CLAWBIO_DISCLAIMER  # noqa: E402
from clawbio.common.reproducibility import (  # noqa: E402
    write_checksums,
    write_commands_sh,
    write_environment_yml,
)

SKILL_NAME = "prs-applicability-gate"
GATE_VERSION = "2.2.0"
INPUT_SCHEMA = "prs-applicability-gate.input.v2"
OUTPUT_SCHEMA = "prs-applicability-gate.output.v2"
DEFAULT_CONFIG = SKILL_DIR / "config" / "calibration.yaml"
DEMO_DIR = SKILL_DIR / "examples"
SUPPORTED, RAW_ONLY, ABSTAIN = "SUPPORTED", "RAW_ONLY", "ABSTAIN"
SEVERITY = {SUPPORTED: 0, RAW_ONLY: 1, ABSTAIN: 2}
GROUPS = ("AFR", "AMR", "EAS", "EUR", "SAS")
# The assemblies PGS Catalog releases harmonised scoring files on, and the only ones upstream harmonisation
# trusts positions for. Anything else (NCBI36, T2T, an unnormalised alias such as hg19) fails G3.
SUPPORTED_BUILDS = ("GRCh37", "GRCh38")

REASON_CODES = {
    "INVALID_GATE_INPUT": "The gate input is missing required fields, has wrong types, or fails its digest.",
    "UNSUPPORTED_SCORE_FORMAT": "The scoring file is not a plain additive log-scale score (ratio weights, "
                                "dosage-specific weights, haplotype/interaction terms, or unreadable rows).",
    "SCORE_FORMAT_UNVERIFIED": "The scoring file's format was not reported (ratio_weight_type not a boolean or "
                               "unsupported_features not a list), so it cannot be confirmed as a plain additive "
                               "log-scale score.",
    "BUILD_UNRESOLVED": "The genotype file's genome build was not established as a supported assembly "
                        "(GRCh37 or GRCh38); positions cannot be trusted.",
    "ALLELE_HARMONIZATION_FAILED": "Too many located variants have alleles that cannot be reconciled with the "
                                   "scoring file (wrong strand/build/file suspected).",
    "LOW_SCOREABILITY": "The score computable from this genotype correlates too weakly with the published score.",
    "SCOREABILITY_UNVERIFIED": "No reference panel or allele frequencies to measure how well the computable "
                               "score represents the published score.",
    "PALINDROMIC_VARIANT_UNRESOLVED": "Excluded strand-ambiguous (A/T, C/G) variants are the largest loss.",
    "VARIANTS_MISSING": "Variants absent from the genotype file are the largest loss.",
    "DUPLICATE_OR_CONFLICTING_VARIANTS": "Duplicated or position-conflicting variants are the largest loss.",
    "SEX_POPULATION_MISMATCH": "The score or all of its relevant evaluations are specific to the other sex.",
    "SEX_NOT_PROVIDED": "The score is sex-specific but the person's sex was not provided.",
    "SEX_EVALUATION_UNVERIFIED": "No informative evaluation is known to include the person's sex (sex composition "
                                 "not reported, or the person's sex not provided and no evaluation includes both "
                                 "sexes).",
    "METADATA_CONTRADICTION": "PGS Catalog metadata for this score failed a blocking consistency check.",
    "EVALUATION_METADATA_UNAVAILABLE": "PGS Catalog metadata for this score could not be retrieved.",
    "TARGET_REFERENCE_UNRESOLVED": "The person could not be placed stably inside one reference group "
                                   "(intermediate/admixed, unstable, or too few sites).",
    "NO_RELEVANT_EVALUATION": "No single-ancestry evaluation in the person's reference group reports a metric.",
    "EVALUATION_NOT_INFORMATIVE": "Relevant evaluations exist but none shows evidence of association in the score's "
                                  "direction (no metric's 95% CI lies entirely above the null).",
    "REFERENCE_DISTRIBUTION_UNAVAILABLE": "No usable reference distribution on the person's matched variant set "
                                          "(none, another group's, or no positive reference size and variant "
                                          "intersection).",
    "REFERENCE_SENSITIVE": "The percentile depends on which reference population is chosen within the group.",
    "REFERENCE_SENSITIVITY_UNVERIFIED": "Whether the percentile depends on the reference population chosen was not "
                                        "reported (reference_sensitive is not a boolean).",
}
LOSS_CODES = {"palindromic_excluded": "PALINDROMIC_VARIANT_UNRESOLVED", "missing": "VARIANTS_MISSING",
              "no_call": "VARIANTS_MISSING", "allele_mismatch": "ALLELE_HARMONIZATION_FAILED",
              "duplicate_excluded": "DUPLICATE_OR_CONFLICTING_VARIANTS",
              "position_conflict": "DUPLICATE_OR_CONFLICTING_VARIANTS", "build_unresolved": "BUILD_UNRESOLVED"}
DISCLAIMER = ("Research software, not a medical device. SUPPORTED is a research-prototype reportability state, not a "
              "clinical recommendation: it requires evidence of association in a relevant evaluation group, which is "
              "not evidence of clinically useful discrimination or calibration. A polygenic score is not a diagnosis, "
              "and this gate never converts a score into absolute risk. Genetic reference placement is not ethnicity "
              "or identity. " + CLAWBIO_DISCLAIMER)
STATUS_MEANING = {SUPPORTED: "raw score, standardised score and reference percentile may be shown",
                  RAW_ONLY: "only the raw score may be shown; no relative statement is supported",
                  ABSTAIN: "nothing may be reported for this score and person"}


# ---------------------------------------------------------------------------
# Config and input validation
# ---------------------------------------------------------------------------


def _sha256_bytes(data: bytes) -> str:
    return "sha256:" + hashlib.sha256(data).hexdigest()


def canonical_digest(obj: Any) -> str:
    return _sha256_bytes(json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=True,
                                    default=str).encode())


def load_config(path: Path = DEFAULT_CONFIG) -> dict:
    raw = Path(path).read_bytes()
    cfg = yaml.safe_load(raw)
    problems = []

    def num(section, key, lo, hi):
        v = (cfg.get(section) or {}).get(key) if isinstance(cfg, dict) else None
        if not isinstance(v, (int, float)) or isinstance(v, bool) or not lo <= v <= hi:
            problems.append(f"{section}.{key}={v!r} not in [{lo}, {hi}]")

    if not isinstance(cfg, dict) or not isinstance(cfg.get("calibration_version"), str):
        raise ValueError("config: calibration_version missing")
    num("scoreability", "r_min", 0.5, 1.0)
    num("allele_harmonisation", "max_mismatch_fraction", 0.0, 0.5)
    num("placement", "min_stability", 0.5, 1.0)
    min_sites = (cfg.get("placement") or {}).get("min_sites")
    if not isinstance(min_sites, int) or isinstance(min_sites, bool) or min_sites < 1:
        problems.append(f"placement.min_sites={min_sites!r} must be a positive integer")
    num("reference_distribution", "interval", 0.5, 0.999)
    if (cfg.get("evaluation") or {}).get("require_ci_excluding_null") is not True:
        problems.append("evaluation.require_ci_excluding_null must be true")
    if problems:
        raise ValueError("config invalid: " + "; ".join(problems))
    cfg["_sha256"] = _sha256_bytes(raw)
    # A non-shipped calibration is allowed for calibration research (direct CLI only) but is always labelled.
    cfg["_canonical"] = cfg["_sha256"] == _sha256_bytes(DEFAULT_CONFIG.read_bytes())
    return cfg


def validate_input(gi: Any) -> list[str]:
    p = []
    if not isinstance(gi, dict):
        return ["input is not a JSON object"]
    if gi.get("schema") != INPUT_SCHEMA:
        p.append(f"schema is {gi.get('schema')!r}, expected {INPUT_SCHEMA}")
    for block in ("candidate", "score_file", "genotype", "person", "harmonisation", "scoreability", "placement",
                  "catalog_metadata", "evaluation"):
        if not isinstance(gi.get(block), dict):
            p.append(f"{block}: missing or not an object")
    if p:
        return p
    if not isinstance(gi["candidate"].get("pgs_id"), str):
        p.append("candidate.pgs_id missing")
    h = gi["harmonisation"]
    for k in ("n_variants", "n_matched", "n_located"):
        if not isinstance(h.get(k), int) or h[k] < 0:
            p.append(f"harmonisation.{k} must be a non-negative integer")
    if isinstance(h.get("n_variants"), int) and isinstance(h.get("n_matched"), int) \
            and h["n_matched"] > h["n_variants"]:
        p.append("harmonisation.n_matched exceeds n_variants")
    # Numeric domains, each checked on its own. allele_mismatch_fraction is null exactly when no scoring variant
    # was located (n_located == 0) and a fraction in [0, 1] otherwise.
    mm = h.get("allele_mismatch_fraction")
    if isinstance(h.get("n_located"), int) and h["n_located"] == 0:
        if mm is not None:
            p.append("harmonisation.allele_mismatch_fraction must be null when n_located == 0")
    elif not _number_in(mm, 0.0, 1.0):
        p.append(f"harmonisation.allele_mismatch_fraction={mm!r} must be a number in [0, 1] when variants were located")
    for k, v, lo, hi in (("scoreability.r", gi["scoreability"].get("r"), -1.0, 1.0),
                         ("placement.placement_stability", gi["placement"].get("placement_stability"), 0.0, 1.0)):
        if v is not None and not _number_in(v, lo, hi):
            p.append(f"{k}={v!r} must be null or a number in [{lo:g}, {hi:g}]")
    loss = h.get("weight_loss_by_status")
    if loss is not None and (not isinstance(loss, dict) or not all(_number_in(x, 0.0, 1.0) for x in loss.values())):
        p.append("harmonisation.weight_loss_by_status must map statuses to fractions in [0, 1]")
    sf = gi["score_file"]
    if not isinstance(sf.get("unsupported_features") or [], list):
        p.append("score_file.unsupported_features must be a list")
    for k in ("n_parse_problems", "variants_interactions"):
        if sf.get(k) is not None and (not isinstance(sf[k], int) or isinstance(sf[k], bool) or sf[k] < 0):
            p.append(f"score_file.{k} must be null or a non-negative integer")
    if gi["candidate"].get("sex_specific") not in (None, "female", "male"):
        p.append("candidate.sex_specific must be female, male or null")
    if not isinstance(gi["catalog_metadata"].get("status"), str):
        p.append("catalog_metadata.status must be a string")
    rd = gi.get("reference_distribution")
    if rd is not None and not isinstance(rd, dict):
        p.append("reference_distribution must be an object or null")
    if gi["person"].get("sex") not in (None, "female", "male"):
        p.append("person.sex must be female, male or null")
    if gi["placement"].get("status") not in ("RESOLVED", "INTERMEDIATE", "UNSTABLE", "UNRESOLVED"):
        p.append("placement.status invalid")
    if gi["placement"].get("status") == "RESOLVED" and gi["placement"].get("placement") not in GROUPS:
        p.append("placement.placement must be a reference group when RESOLVED")
    units = gi["evaluation"].get("units")
    if not isinstance(units, list):
        p.append("evaluation.units must be a list")
    else:
        for i, u in enumerate(units):
            if not isinstance(u, dict) or not isinstance(u.get("metrics", []), list) \
                    or not all(isinstance(m, dict) for m in u.get("metrics", [])):
                p.append(f"evaluation.units[{i}] must be an object whose metrics are a list of objects")
                break
    if "input_digest" in gi:
        body = {k: v for k, v in gi.items() if k != "input_digest"}
        if canonical_digest(body) != gi["input_digest"]:
            p.append("input_digest does not match the input (modified after it was built)")
    return p


def _finite_number(v: Any) -> bool:
    """A real number: not a bool, NaN or +/-inf. Ints are never converted to float, so a huge int cannot overflow."""
    if isinstance(v, bool) or not isinstance(v, (int, float)):
        return False
    return isinstance(v, int) or math.isfinite(v)


def _number_in(v: Any, lo: float, hi: float) -> bool:
    return _finite_number(v) and lo <= v <= hi


def _positive_int(v: Any) -> bool:
    return isinstance(v, int) and not isinstance(v, bool) and v > 0


def _is_bool(v: Any) -> bool:
    return isinstance(v, bool)


def _is_list(v: Any) -> bool:
    return isinstance(v, list)


def _supported_build(v: Any) -> bool:
    return isinstance(v, str) and v in SUPPORTED_BUILDS


# ---------------------------------------------------------------------------
# Rules
# ---------------------------------------------------------------------------


G1_QUESTION = "Is the gate input well-formed and unmodified?"
G1_CHANGE = (f"a well-formed gate input conforming to schema {INPUT_SCHEMA} (SKILL.md, Input Formats), "
             "built from the actual evidence")
G9_QUESTION = ("Is there evidence of association (95% CI above the null) in an evaluation of the person's group? "
               "(Association only: not clinically useful discrimination or calibration.)")


class Trace:
    def __init__(self):
        self.rules: list[dict] = []
        self.used: dict[str, Any] = {}
        self.missing: list[str] = []

    def use(self, path: str, value: Any, needed: bool = True, *,
            usable: Callable[[Any], bool] | None = None) -> Any:
        """Record a field a rule read. A needed field that is null, or fails `usable`, is missing evidence."""
        self.used[path] = value
        unusable = value is None or (usable is not None and not usable(value))
        if unusable and needed and path not in self.missing:
            self.missing.append(path)
        return value

    def add(self, rule: str, name: str, question: str, outcome: str, effect: str | None = None,
            codes: list[str] | None = None, detail: str = "", change: str | None = None):
        self.rules.append({"rule": rule, "name": name, "question": question, "outcome": outcome,
                           "effect": effect if outcome == "fail" else None, "codes": codes or [],
                           "detail": detail, "what_would_change": change if outcome == "fail" else None})


def _metric_ok(m: dict) -> bool | None:
    """True if the 95% CI lies entirely above the null (association in the direction the score is built for:
    OR/HR/RR > 1, beta > 0, AUROC/C-index > 0.5, r > 0); None if it cannot be assessed.

    An inverse association (CI entirely below the null, e.g. a case-only subtype comparison) is not evidence
    for interpreting the score. Evidence of association is NOT evidence of clinically useful discrimination or
    calibration."""
    lo, hi, null = m.get("ci_lower"), m.get("ci_upper"), m.get("null")
    if not all(_finite_number(x) for x in (lo, hi, null)) or lo > hi:
        return None
    return lo > null


def _includes_sex(percent_male: float, sex: str | None) -> bool:
    """Whether a unit with this known percent_male includes participants of `sex` (both sexes when sex is None)."""
    if sex == "female":
        return percent_male < 100
    if sex == "male":
        return percent_male > 0
    return 0 < percent_male < 100


def _largest_losses(loss: dict) -> list[str]:
    """The status(es) with the largest weight loss, whatever the input's key order: every status tied at the
    maximum, in name order; none when nothing was lost. Values are validated fractions in [0, 1]."""
    most = max(loss.values(), default=0)
    return sorted(k for k, v in loss.items() if v == most) if most > 0 else []


def invalid_input(gi: Any, cfg: dict, detail: str) -> dict:
    """ABSTAIN with INVALID_GATE_INPUT: the fail-closed decision for anything that is not a usable gate input."""
    t = Trace()
    t.add("G1", "INPUT_VALID", G1_QUESTION, "fail", ABSTAIN, ["INVALID_GATE_INPUT"], detail, G1_CHANGE)
    return _finish(gi if isinstance(gi, dict) else {}, cfg, t)


def evaluate(gi: Any, cfg: dict) -> dict:
    """Deterministic decision for one gate input. Never raises on malformed input: anything the validator did not
    anticipate still fails closed as INVALID_GATE_INPUT (ABSTAIN)."""
    try:
        return _evaluate(gi, cfg)
    except (TypeError, KeyError, AttributeError, ValueError, IndexError, RecursionError) as exc:
        return invalid_input(gi, cfg, f"malformed input ({type(exc).__name__} while evaluating rules)")


def _evaluate(gi: dict, cfg: dict) -> dict:
    t = Trace()
    problems = validate_input(gi)
    if problems:
        return invalid_input(gi, cfg, "; ".join(problems))
    t.add("G1", "INPUT_VALID", G1_QUESTION, "pass", detail="schema v2")
    sf, gt, h, sc = gi["score_file"], gi["genotype"], gi["harmonisation"], gi["scoreability"]
    pl, md, person = gi["placement"], gi["catalog_metadata"], gi["person"]
    rd = gi.get("reference_distribution") or {}

    # G2 score format ------------------------------------------------------------------------------------
    # The weight scale and the feature list must be reported: unknown is not evidence of a plain additive score.
    features = t.use("score_file.unsupported_features", sf.get("unsupported_features"), usable=_is_list)
    ratio = t.use("score_file.ratio_weight_type", sf.get("ratio_weight_type"), usable=_is_bool)
    unverified = [f"{k} {v!r} not reported as {kind}" for k, v, kind, ok in (
        ("ratio_weight_type", ratio, "true/false", _is_bool(ratio)),
        ("unsupported_features", features, "a list", _is_list(features))) if not ok]
    unsup = list(features) if _is_list(features) else []
    if ratio is True:
        unsup.append(f"ratio-scale weight_type {sf.get('weight_type')!r}")
    inter = t.use("score_file.variants_interactions", sf.get("variants_interactions"), False)
    if isinstance(inter, int) and inter > 0:
        unsup.append(f"{inter} interaction terms")
    bad_rows = t.use("score_file.n_parse_problems", sf.get("n_parse_problems"), False)
    if isinstance(bad_rows, int) and bad_rows > 0:
        unsup.append(f"{bad_rows} unreadable scoring rows")
    codes = (["UNSUPPORTED_SCORE_FORMAT"] if unsup else []) + (["SCORE_FORMAT_UNVERIFIED"] if unverified else [])
    changes = (["a scoring file without the unsupported features"] if unsup else []) + (
        [("an evidence builder that reports score_file.ratio_weight_type (true/false) and "
          "score_file.unsupported_features (a list)")] if unverified else [])
    t.add("G2", "SCORE_FORMAT", "Is the scoring file a plain additive score on a log scale?",
          "fail" if codes else "pass", ABSTAIN, codes,
          "; ".join(unsup + unverified) or f"additive; weight_type {sf.get('weight_type')!r}",
          "; ".join(changes))

    # G3 build ---------------------------------------------------------------------------------------------
    build = t.use("genotype.build", gt.get("build"), usable=_supported_build)
    build_ok = _supported_build(build)
    method = (gt.get("build_evidence") or {}).get("method")
    t.add("G3", "BUILD", "Is the genotype file's genome build established (never assumed)?",
          "pass" if build_ok else "fail", ABSTAIN, ["BUILD_UNRESOLVED"],
          f"build {build}; method {method}" if build_ok else
          f"build {build!r} is not established as {' or '.join(SUPPORTED_BUILDS)}; method {method}",
          f"a genotype file on {' or '.join(SUPPORTED_BUILDS)} whose build is declared in its header or "
          "verifiable from rsID positions")

    # G4 alleles ---------------------------------------------------------------------------------------------
    mm = t.use("harmonisation.allele_mismatch_fraction", h.get("allele_mismatch_fraction"), False)
    cut = cfg["allele_harmonisation"]["max_mismatch_fraction"]
    if h["n_located"] == 0:
        t.add("G4", "ALLELES", "Are located variants' alleles consistent with the scoring file?", "not_applicable",
              detail="no scoring variant located in the genotype file")
    else:
        t.add("G4", "ALLELES", "Are located variants' alleles consistent with the scoring file?",
              "fail" if mm > cut else "pass", ABSTAIN, ["ALLELE_HARMONIZATION_FAILED"],
              f"{mm:.1%} of {h['n_located']} located variants unreconcilable (max {cut:.0%})",
              "a genotype file on the forward strand of the declared build")

    # G5 scoreability ------------------------------------------------------------------------------------------
    r = t.use("scoreability.r", sc.get("r"))
    r_min = cfg["scoreability"]["r_min"]
    loss = h.get("weight_loss_by_status") or {}
    t.use("harmonisation.weight_loss_by_status", loss, False)
    if r is None:
        t.add("G5", "SCOREABILITY", "Does the computable score represent the published score (r >= r_min)?",
              "fail", ABSTAIN, ["SCOREABILITY_UNVERIFIED"], "no reference panel scores and no allele frequencies",
              "reference-panel genotypes (or reported allele frequencies) for this score's variants")
    else:
        codes = ["LOW_SCOREABILITY"]
        top = _largest_losses(loss)
        if r < r_min:
            codes += list(dict.fromkeys(LOSS_CODES[s] for s in top if s in LOSS_CODES))
        t.add("G5", "SCOREABILITY", "Does the computable score represent the published score (r >= r_min)?",
              "fail" if r < r_min else "pass", ABSTAIN, codes,
              f"r = {r:.3f} ({sc.get('method')}); r_min {r_min}; {h['n_matched']}/{h['n_variants']} variants "
              f"matched; largest loss: {', '.join(top) or 'none'}",
              f"genotype data covering the score's heavily weighted variants so that r >= {r_min} "
              "(e.g. sequencing or imputation instead of a sparse array)")

    # G6 score sex-specificity --------------------------------------------------------------------------------
    sex = t.use("person.sex", person.get("sex"), False)
    specific = t.use("candidate.sex_specific", gi["candidate"].get("sex_specific"), False)
    if specific and sex and specific != sex:
        t.add("G6", "SEX_SCORE", "Is a sex-specific score applied to a person of that sex?", "fail", ABSTAIN,
              ["SEX_POPULATION_MISMATCH"], f"{specific}-specific score, person is {sex}",
              "a score not restricted to the other sex")
    elif specific and not sex:
        t.add("G6", "SEX_SCORE", "Is a sex-specific score applied to a person of that sex?", "fail", RAW_ONLY,
              ["SEX_NOT_PROVIDED"], f"{specific}-specific score; person's sex not provided",
              "the person's sex (only needed because this score is sex-specific)")
    else:
        t.add("G6", "SEX_SCORE", "Is a sex-specific score applied to a person of that sex?",
              "pass" if specific else "not_applicable", detail=f"score specificity {specific or 'none'}")

    # G7 catalog metadata -------------------------------------------------------------------------------------
    status = t.use("catalog_metadata.status", md.get("status"))
    if status == "resolved":
        t.add("G7", "METADATA", "Is the PGS Catalog record resolved and internally consistent?", "pass",
              detail="all blocking consistency checks passed")
    else:
        code = "METADATA_CONTRADICTION" if status == "contradictory" else "EVALUATION_METADATA_UNAVAILABLE"
        t.add("G7", "METADATA", "Is the PGS Catalog record resolved and internally consistent?", "fail", RAW_ONLY,
              [code], f"{status}: {md.get('detail')}",
              "a consistent PGS Catalog record (the discrepancy reported to the Catalog and corrected)")

    # G8 placement -------------------------------------------------------------------------------------------
    # A RESOLVED status is not taken on trust: it must rest on enough placement sites and a stable placement, at
    # the calibrated thresholds. Low-marker over-calls are the commonest ancestry-placement failure.
    pst = t.use("placement.status", pl.get("status"))
    claimed = pst == "RESOLVED"
    group = t.use("placement.placement", pl.get("placement"), claimed)
    sites = t.use("placement.n_sites_used", pl.get("n_sites_used"), claimed)
    stability = t.use("placement.placement_stability", pl.get("placement_stability"), claimed)
    min_sites, min_stability = cfg["placement"]["min_sites"], cfg["placement"]["min_stability"]
    shortfalls = []
    if claimed:
        if not (isinstance(sites, int) and not isinstance(sites, bool) and sites >= min_sites):
            shortfalls.append(f"n_sites_used={sites!r}, needs an integer >= min_sites {min_sites}")
        if not (_number_in(stability, 0.0, 1.0) and stability >= min_stability):
            shortfalls.append(f"placement_stability={stability!r}, needs >= min_stability {min_stability}")
    resolved = claimed and not shortfalls
    t.add("G8", "PLACEMENT", "Is the person placed stably inside one reference group?",
          "pass" if resolved else "fail", RAW_ONLY, ["TARGET_REFERENCE_UNRESOLVED"],
          f"{pst}: {pl.get('detail')}" + (f"; RESOLVED not supported: {'; '.join(shortfalls)}" if shortfalls else ""),
          f"a placement RESOLVED from at least {min_sites} sites with stability >= {min_stability}, against a "
          "reference panel that represents this person's genetic background (e.g. admixed references with "
          "local-ancestry-aware scoring); not something the gate can relax")

    # G9 evaluation ------------------------------------------------------------------------------------------
    units = gi["evaluation"]["units"]
    if not resolved:
        t.add("G9", "EVALUATION", G9_QUESTION,
              "not_applicable", detail="no resolved reference group to match evaluations against")
        t.add("G10", "SEX_EVALUATION", "Do the informative evaluations include the person's sex?",
              "not_applicable", detail="no resolved reference group")
    else:
        rel = [u for u in units if u.get("code") == group and u.get("pooled") is False]
        with_metric = [u for u in rel if u.get("metrics")]
        informative = [u for u in with_metric if any(_metric_ok(m) for m in u["metrics"])]
        t.use("evaluation.relevant_units", len(rel), False)
        t.use("evaluation.relevant_units_with_metric", len(with_metric), False)
        t.use("evaluation.relevant_units_informative", len(informative), False)
        n_rel = sum(u.get("n") or 0 for u in rel)
        if not with_metric:
            t.add("G9", "EVALUATION", G9_QUESTION,
                  "fail", RAW_ONLY, ["NO_RELEVANT_EVALUATION"],
                  f"0 single-ancestry {group} evaluation units with a reported metric "
                  f"({len(units)} units in total)",
                  f"a published evaluation of this score in a {group} cohort, deposited in the PGS Catalog")
        elif not informative:
            t.add("G9", "EVALUATION", G9_QUESTION,
                  "fail", RAW_ONLY, ["EVALUATION_NOT_INFORMATIVE"],
                  f"{len(with_metric)} {group} units with metrics, none with a 95% CI entirely above the null",
                  f"a {group} evaluation showing association in the score's direction")
        else:
            t.add("G9", "EVALUATION", G9_QUESTION,
                  "pass", detail=f"{len(informative)} {group} evaluation units with a 95% CI above the null "
                                 f"({len(with_metric)} with metrics; {n_rel:,} individuals)")
        # G10 evaluation sex -----------------------------------------------------------------------------
        # Unknown sex composition never counts as including the person, and a person whose sex is not provided
        # is covered only by an evaluation that includes both sexes.
        if informative:
            known = [p for p in (u.get("percent_male") for u in informative) if _number_in(p, 0.0, 100.0)]
            n_unknown = t.use("evaluation.informative_units_sex_unknown", len(informative) - len(known), False)
            n_cover = sum(_includes_sex(p, sex) for p in known)
            who = f"{sex} participants" if sex else "both sexes"
            if n_cover:
                t.add("G10", "SEX_EVALUATION", "Do the informative evaluations include the person's sex?", "pass",
                      detail=f"{n_cover} of {len(informative)} informative {group} evaluations include {who}; "
                             f"{n_unknown} of unknown sex composition")
            elif sex and not n_unknown:
                t.add("G10", "SEX_EVALUATION", "Do the informative evaluations include the person's sex?", "fail",
                      RAW_ONLY, ["SEX_POPULATION_MISMATCH"],
                      f"all {len(informative)} informative {group} evaluations contain only the other sex",
                      f"an evaluation including {sex} participants")
            else:
                if not sex:
                    t.use("person.sex", sex)
                t.add("G10", "SEX_EVALUATION", "Do the informative evaluations include the person's sex?", "fail",
                      RAW_ONLY, ["SEX_EVALUATION_UNVERIFIED"],
                      f"no informative {group} evaluation is known to include {who} ({n_unknown} of "
                      f"{len(informative)} of unknown sex composition{'' if sex else '; person sex not provided'})",
                      f"an informative {group} evaluation that reports its sex composition and includes {who}"
                      + ("" if sex else ", or the person's sex"))
        else:
            t.add("G10", "SEX_EVALUATION", "Do the informative evaluations include the person's sex?",
                  "not_applicable", detail="no informative evaluation")

    # G11/G12 reference distribution -----------------------------------------------------------------------------
    if not resolved:
        t.add("G11", "REFERENCE_DISTRIBUTION", "Is there a reference distribution on the matched variants?",
              "not_applicable", detail="no resolved reference group")
        t.add("G12", "REFERENCE_SENSITIVITY", "Is the percentile robust to the reference population chosen?",
              "not_applicable", detail="no resolved reference group")
    else:
        # Structural validity only (the distribution exists for this group, on a non-empty variant set); no
        # sample-size threshold is applied here.
        avail = t.use("reference_distribution.available", rd.get("available"))
        claimed = avail is True
        ref_group = t.use("reference_distribution.reference_group", rd.get("reference_group"), claimed)
        ref_n = t.use("reference_distribution.reference_n", rd.get("reference_n"), claimed, usable=_positive_int)
        n_inter = t.use("reference_distribution.n_intersection", rd.get("n_intersection"), claimed,
                        usable=_positive_int)
        if not claimed:
            problems = [f"unavailable: {rd.get('detail')}"]
        else:
            problems = [msg for bad, msg in (
                (ref_group != group, f"reference group {ref_group!r} is not the placed group {group}"),
                (not _positive_int(ref_n), f"reference_n {ref_n!r} is not a positive integer"),
                (not _positive_int(n_inter), f"n_intersection {n_inter!r} is not a positive integer")) if bad]
        ok = not problems
        t.add("G11", "REFERENCE_DISTRIBUTION", "Is there a reference distribution on the matched variants?",
              "pass" if ok else "fail", RAW_ONLY, ["REFERENCE_DISTRIBUTION_UNAVAILABLE"],
              f"{ref_n} {group} reference individuals scored on {n_inter} matched variants" if ok else
              "; ".join(problems),
              "reference-panel genotypes at this score's variants")
        sens = t.use("reference_distribution.reference_sensitive", rd.get("reference_sensitive"), ok,
                     usable=_is_bool)
        if ok and sens is True:
            pairs = rd.get("reference_sensitive_pairs") or []
            t.add("G12", "REFERENCE_SENSITIVITY", "Is the percentile robust to the reference population chosen?",
                  "fail", RAW_ONLY, ["REFERENCE_SENSITIVE"],
                  "disjoint 95% intervals between " + ", ".join("/".join(p) for p in pairs),
                  "a reference distribution matched more finely to the person (not a choice the agent may make "
                  "to obtain a preferred percentile)")
        elif ok and sens is not False:
            t.add("G12", "REFERENCE_SENSITIVITY", "Is the percentile robust to the reference population chosen?",
                  "fail", RAW_ONLY, ["REFERENCE_SENSITIVITY_UNVERIFIED"],
                  f"reference_sensitive {sens!r} is not reported as true/false",
                  "a reference distribution that reports whether the percentile changes with the reference "
                  "population (reference_sensitive true or false)")
        elif ok:
            t.add("G12", "REFERENCE_SENSITIVITY", "Is the percentile robust to the reference population chosen?",
                  "pass", detail=f"percentile intervals overlap across {group} reference populations")
        else:
            t.add("G12", "REFERENCE_SENSITIVITY", "Is the percentile robust to the reference population chosen?",
                  "not_applicable", detail="no reference distribution")
    return _finish(gi, cfg, t)


def _finish(gi: Any, cfg: dict, t: Trace) -> dict:
    failed = [r for r in t.rules if r["outcome"] == "fail"]
    status = max((r["effect"] for r in failed), key=SEVERITY.get, default=SUPPORTED)
    primary = next((r for r in failed if r["effect"] == status), None)
    codes: list[str] = []
    for r in failed:
        for c in r["codes"]:
            if c not in codes:
                codes.append(c)
    allowed = {"raw_score": status in (SUPPORTED, RAW_ONLY), "standardized_score": status == SUPPORTED,
               "percentile": status == SUPPORTED, "absolute_risk": False}
    body = gi if isinstance(gi, dict) else {}
    return {
        "schema": OUTPUT_SCHEMA,
        "pgs_id": (body.get("candidate") or {}).get("pgs_id") if isinstance(body.get("candidate"), dict) else None,
        "status": status,
        "primary_reason": None if primary is None else {"code": primary["codes"][0], "rule": primary["rule"],
                                                        "detail": primary["detail"],
                                                        "meaning": REASON_CODES[primary["codes"][0]]},
        "reason_codes": codes,
        "allowed_claims": allowed,
        "evidence_used": t.used,
        "evidence_missing": t.missing,
        "calibration_version": cfg.get("calibration_version"),
        "rule_trace": t.rules,
        "what_would_change_result": [{"rule": r["rule"], "code": r["codes"][0], "change": r["what_would_change"]}
                                     for r in failed],
        "provenance": {"gate": SKILL_NAME, "gate_version": GATE_VERSION,
                       "config_sha256": cfg.get("_sha256"),
                       "config_canonical": cfg.get("_canonical"),
                       "input_digest": canonical_digest({k: v for k, v in body.items() if k != "input_digest"}),
                       "deterministic": True},
        "disclaimer": DISCLAIMER,
    }


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

NON_CANONICAL_NOTE = ("**NON-CANONICAL CALIBRATION**: this result used a configuration other than the shipped "
                      "config/calibration.yaml (calibration research only; not a canonical gate result).")


def _cell(value: Any) -> str:
    """Text from the input must not break a markdown table row: collapse whitespace, escape pipes."""
    return " ".join(str(value).split()).replace("|", "\\|")


def render_report(result: dict, level: int = 1, label: str | None = None, disclaimer: bool = True) -> str:
    """Markdown for one decision: status, allowed claims, rule trace and what would change the result."""
    h = "#" * level
    title = label or result["pgs_id"] or "unidentified candidate"
    lines = [f"{h} PRS applicability gate: {_cell(title)}", "",
             f"**Status: {result['status']}**" + (f" ({result['primary_reason']['code']})"
                                                  if result["primary_reason"] else ""), ""]
    if result["primary_reason"]:
        lines += [result["primary_reason"]["meaning"], "", f"Detail: {_cell(result['primary_reason']['detail'])}",
                  ""]
    claims = result["allowed_claims"]
    lines += [f"{h}# What may be reported", "",
              f"- Raw score: {'yes' if claims['raw_score'] else 'no'}",
              f"- Standardised score and percentile: {'yes' if claims['percentile'] else 'no'}",
              "- Absolute risk: never", "", f"{h}# Rule trace", "", "| Rule | Question | Outcome | Detail |",
              "|---|---|---|---|"]
    for r in result["rule_trace"]:
        lines.append(f"| {r['rule']} {r['name']} | {r['question']} | {r['outcome']}"
                     f"{' -> ' + r['effect'] if r['effect'] else ''} | {_cell(r['detail'])} |")
    if result["what_would_change_result"]:
        lines += ["", f"{h}# What would change the result", ""]
        lines += [f"- {w['code']}: {_cell(w['change'])}" for w in result["what_would_change_result"]]
    if result["provenance"].get("config_canonical") is False:
        lines += ["", NON_CANONICAL_NOTE]
    lines += ["", f"Calibration {result['calibration_version']}; gate {result['provenance']['gate_version']}; "
                  f"input {result['provenance']['input_digest']}", ""]
    if disclaimer:
        lines += [f"*{result['disclaimer']}*", ""]
    return "\n".join(lines)


def render_run_report(rows: list[dict], decisions: list[dict], cfg: dict, demo: bool) -> str:
    """The run's report.md: a summary table, then the full rule trace of every decision."""
    canonical = "canonical" if cfg["_canonical"] else "NON-CANONICAL"
    if demo:
        source = (f"**Mode**: demo, {len(rows)} SYNTHETIC gate inputs bundled with the skill (invented values; "
                  "no person, no genotype, no real PGS Catalog record)")
    else:
        source = f"**Input**: `{_cell(rows[0]['input'])}`"
    lines = ["# PRS Applicability Gate Report", "", source,
             f"**Calibration**: `{cfg['calibration_version']}` ({canonical}; config {cfg['_sha256']}) · "
             f"**Gate**: {GATE_VERSION}", ""]
    if not cfg["_canonical"]:
        lines += [NON_CANONICAL_NOTE, ""]
    lines += ["| Input | Score | Status | Primary reason | Reason codes |", "|---|---|---|---|---|"]
    for row in rows:
        lines.append(f"| {_cell(row['input'])} | {_cell(row['pgs_id'] or '-')} | **{row['status']}** | "
                     f"{row['primary_reason'] or '-'} | {', '.join(row['reason_codes']) or '-'} |")
    lines += ["", "; ".join(f"**{s}**: {m}" for s, m in STATUS_MEANING.items()) + ".", ""]
    for row, decision in zip(rows, decisions):
        lines += [render_report(decision, level=2, label=f"{decision['pgs_id'] or row['input']}", disclaimer=False)]
    lines += ["---", "", f"*{DISCLAIMER}*", ""]
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def read_gate_input(path: Path, cfg: dict) -> dict:
    """Evaluate one gate input file. Content that is not UTF-8 JSON fails closed (ABSTAIN), it never crashes."""
    try:
        gi = json.loads(Path(path).read_text(encoding="utf-8-sig"))
    except (ValueError, RecursionError) as exc:  # JSONDecodeError and UnicodeDecodeError are ValueErrors
        return invalid_input(None, cfg, f"input file is not valid UTF-8 JSON ({type(exc).__name__})")
    return evaluate(gi, cfg)


def _write_bundle(out: Path, args: argparse.Namespace, cfg: dict) -> None:
    command = ["python", "skills/prs-applicability-gate/prs_applicability_gate.py"]
    command += ["--demo"] if args.demo else ["--input", str(args.input.resolve())]
    command += ["--output", str(out.resolve())]
    if not cfg["_canonical"]:
        command += ["--config", str(Path(args.config).resolve())]
    write_commands_sh(out, "# Run from the ClawBio repository root.\n" + shlex.join(command))
    env = write_environment_yml(out, env_name="clawbio-prs-applicability-gate", pip_deps=["pyyaml>=6.0"])
    write_checksums([out / "report.md", out / "result.json", out / "reproducibility" / "commands.sh", env],
                    out, anchor=out)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    source = ap.add_mutually_exclusive_group(required=True)
    source.add_argument("--input", type=Path, help=f"one gate input JSON file (schema {INPUT_SCHEMA})")
    source.add_argument("--demo", action="store_true", help="evaluate the bundled SYNTHETIC example inputs")
    ap.add_argument("--output", type=Path, required=True, help="output directory")
    ap.add_argument("--config", type=Path, default=DEFAULT_CONFIG,
                    help="EXPERT/DEBUG ONLY: an alternative calibration file for calibration research. Results are "
                         "labelled NON-CANONICAL CALIBRATION (config_canonical: false). Not forwarded by "
                         "`clawbio.py run prs-gate`, which always uses the shipped config/calibration.yaml.")
    args = ap.parse_args(argv)

    try:
        cfg = load_config(args.config)
    except (OSError, ValueError, yaml.YAMLError) as exc:
        ap.error(f"cannot use calibration config {args.config}: {exc}")  # exits with status 2
        return 2  # not reached; makes the exit explicit, so cfg is provably bound below
    if not cfg["_canonical"]:
        print(f"WARNING: NON-CANONICAL CALIBRATION: {args.config.name} is not the shipped config/calibration.yaml; "
              "every result is labelled config_canonical: false (calibration research only).", file=sys.stderr)
    if args.demo:
        inputs = sorted(DEMO_DIR.glob("*.gate_input.json"))
    else:
        if not args.input.is_file():
            ap.error(f"input file not found: {args.input}")
        inputs = [args.input]
    out = args.output
    if out.exists() and not out.is_dir():
        ap.error(f"output path exists and is not a directory: {out}")
    out.mkdir(parents=True, exist_ok=True)
    existing = [name for name in ("report.md", "result.json") if (out / name).exists()]
    if existing:
        print(f"WARNING: overwriting existing {', '.join(existing)} in {out}", file=sys.stderr)

    decisions, rows, datasets = [], [], {}
    for path in inputs:
        decision = read_gate_input(path, cfg)
        decisions.append(decision)
        datasets[path.name] = _sha256_bytes(path.read_bytes())
        rows.append({"input": path.name, "pgs_id": decision["pgs_id"], "status": decision["status"],
                     "primary_reason": (decision["primary_reason"] or {}).get("code"),
                     "reason_codes": decision["reason_codes"]})
        print(f"{path.name}: {decision['status']} {decision['reason_codes']}")

    summary = {"mode": "demo" if args.demo else "input", "synthetic_demo": bool(args.demo),
               "n_inputs": len(rows),
               "status_counts": {s: sum(r["status"] == s for r in rows) for s in (SUPPORTED, RAW_ONLY, ABSTAIN)},
               "calibration_version": cfg["calibration_version"], "gate_version": GATE_VERSION,
               "config_sha256": cfg["_sha256"], "config_canonical": cfg["_canonical"], "decisions": rows}
    # The ClawBio result.json envelope, without a completion timestamp: identical runs give byte-identical output.
    result = {"skill": SKILL_NAME, "version": GATE_VERSION,
              "input_checksum": "" if args.demo else datasets[inputs[0].name],
              "datasets": datasets, "summary": summary, "data": {"decisions": decisions},
              "chat_summary_lines": [f"{r['pgs_id'] or r['input']}: {r['status']}"
                                     + (f" ({r['primary_reason']})" if r["primary_reason"] else "") for r in rows],
              "disclaimer": DISCLAIMER}
    (out / "report.md").write_text(render_run_report(rows, decisions, cfg, args.demo), encoding="utf-8")
    (out / "result.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    _write_bundle(out, args, cfg)
    return 0


if __name__ == "__main__":
    sys.exit(main())
