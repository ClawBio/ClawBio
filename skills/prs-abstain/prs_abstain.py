#!/usr/bin/env python3
"""prs-abstain — gate polygenic score percentiles on ancestry transferability.

One job: decide whether a PRS percentile is interpretable for a given individual,
and refuse to report it when it is not. This skill does not infer ancestry
(see `ancestry-risk-profiler`, `claw-ancestry-pca`) and does not compute polygenic
scores (see `gwas-prs`). It sits between them.

ClawBio is a research and educational tool. It is not a medical device.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import statistics as st
import sys
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Sequence

__version__ = "0.3.0"

SKILL_DIR = Path(__file__).resolve().parent
EXAMPLES = SKILL_DIR / "examples"
DEFAULT_PCS = ("PC1", "PC2", "PC3", "PC4")

DISCLAIMER = (
    "ClawBio is a research and educational tool. It is not a medical device and "
    "does not provide clinical diagnoses. Consult a healthcare professional "
    "before making any medical decisions."
)
NOT_REASSURANCE = (
    "A withheld percentile is **not evidence of low risk**. It means this score "
    "cannot be interpreted for this individual, not that their risk is average or low."
)
# Kosoy et al. 2009 (Hum Mutat 30:69-78): lower bound of markers for reliable
# continental assignment. Same bound used by ancestry-risk-profiler.
DEFAULT_MIN_MARKERS = 30


class CalibrationError(ValueError):
    """Raised when a reference population cannot support a defensible threshold."""


@dataclass
class Individual:
    sample_id: str
    population: str | None
    pcs: Sequence[float] | None
    n_markers_shared: int | None
    sex: str | None = None


@dataclass
class Calibration:
    reference_population: str
    pcs_used: tuple[str, ...]
    centroid: list[float]
    n: int
    mean: float
    sd: float
    k_sd: float
    threshold: float
    within_max: float
    nearest_other: float | None
    other_populations: dict[str, float] = field(default_factory=dict)

    @property
    def threshold_exceeds_nearest_other(self) -> bool:
        """True when the threshold is wide enough to admit non-reference individuals."""
        return self.nearest_other is not None and self.threshold > self.nearest_other


@dataclass
class Decision:
    sample_id: str
    verdict: str
    distance: float | None
    threshold: float
    reason: str
    remedy: str
    n_markers_shared: int | None
    declared_population: str | None



# ══════════════════════════════════════════════════════════════════════════════
# v0.2 — Applicability, score integrity, and per-variant allele-frequency audit
#
# Mechanism this section exists to expose: for every curated score in the demo,
# the EUR reference mean equals Sum(2 * AF_EUR * w) to within rounding. The
# percentile is therefore centred on a European allele-frequency calculation.
# Re-centring on another population's frequencies shifts the mean by
# Sum(2 * (AF_ref - AF_pop) * w), and that sum decomposes variant by variant.
# ══════════════════════════════════════════════════════════════════════════════

PALINDROMIC = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}

# Traits whose scores are derived in, and clinically meaningful for, one sex only.
SEX_SPECIFIC = {
    "female": ("breast", "ovarian", "cervical", "endometrial", "uterine"),
    "male": ("prostate", "testicular"),
}


@dataclass
class ScoreDefinition:
    pgs_id: str
    trait: str
    build: str | None
    variants: list[dict[str, Any]]


@dataclass
class Applicability:
    applicable: bool
    reason: str


@dataclass
class ScoreAudit:
    pgs_id: str
    n_total: int
    n_matched: int
    weight_total: float
    weight_covered: float
    weight_coverage: float
    weight_at_risk: float
    top1_share: float
    effective_n: float
    palindromic_n: int
    palindromic_share: float
    missing_top: list[dict[str, Any]] = field(default_factory=list)


@dataclass
class IntegrityVerdict:
    passed: bool
    reasons: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)


@dataclass
class AFShift:
    pgs_id: str
    population: str
    n_variants_with_af: int
    coverage: float
    shift_raw: float
    shift_sd: float
    per_variant: list[dict[str, Any]] = field(default_factory=list)


def load_score_definitions(path: Path) -> dict[str, ScoreDefinition]:
    """Load PGS Catalog scoring files from a directory."""
    out: dict[str, ScoreDefinition] = {}
    for fp in sorted(Path(path).glob("*.txt")):
        hdr: dict[str, str] = {}
        variants: list[dict[str, Any]] = []
        for line in fp.read_text().splitlines():
            if line.startswith("#"):
                if "=" in line:
                    k, v = line[1:].split("=", 1)
                    hdr[k.strip()] = v.strip()
                continue
            if line.startswith("rsID"):
                continue
            parts = line.split("\t")
            if len(parts) < 6:
                continue
            af = None
            if len(parts) > 6 and parts[6].strip():
                try:
                    af = float(parts[6])
                except ValueError:
                    af = None
            variants.append({
                "rsid": parts[0], "chr": parts[1], "pos": parts[2],
                "effect_allele": parts[3].strip().upper(),
                "other_allele": parts[4].strip().upper(),
                "weight": float(parts[5]), "af_reference": af,
            })
        if variants:
            pid = hdr.get("pgs_id", fp.stem.split("_")[0])
            out[pid] = ScoreDefinition(pid, hdr.get("trait_reported", "unknown"),
                                       hdr.get("genome_build"), variants)
    return out


def load_genotype(path: Path) -> dict[str, str]:
    """Load a 23andMe/AncestryDNA-style genotype file into {rsid: genotype}."""
    out: dict[str, str] = {}
    for line in Path(path).read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) >= 4:
            out[parts[0]] = parts[3].strip().upper()
    return out


def load_population_af(path: Path) -> dict[str, dict[str, float]]:
    """Load per-population effect-allele frequencies: {rsid: {pop: af}}."""
    out: dict[str, dict[str, float]] = {}
    for line in Path(path).read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 3 or parts[0].lower() == "rsid":
            continue
        try:
            out.setdefault(parts[0], {})[parts[1].strip().upper()] = float(parts[2])
        except ValueError:
            continue
    return out


def check_applicability(score: dict[str, Any] | ScoreDefinition, sex: str | None) -> Applicability:
    """Refuse sex-specific scores for the wrong or unstated sex."""
    trait = (score.trait if isinstance(score, ScoreDefinition) else score.get("trait", "")) or ""
    t = trait.lower()
    sex_norm = (sex or "").strip().lower() or None
    if sex_norm in ("unknown", "unspecified", "na", ""):
        sex_norm = None
    for required_sex, keywords in SEX_SPECIFIC.items():
        if any(k in t for k in keywords):
            if sex_norm is None:
                return Applicability(False, (
                    f"{trait} is a sex-specific trait and no sex was recorded for this "
                    f"individual. The score cannot be shown without knowing it applies."))
            if sex_norm != required_sex:
                return Applicability(False, (
                    f"{trait} is derived in {required_sex} cohorts, but this individual is "
                    f"recorded as {sex_norm}. The score does not apply."))
            return Applicability(True, f"{trait} applies to {sex_norm} individuals.")
    return Applicability(True, "Trait is not sex-specific.")


def expected_mean(score: ScoreDefinition, af_key: str = "af_reference",
                  af_table: dict[str, dict[str, float]] | None = None,
                  population: str | None = None) -> float | None:
    """Sum(2 * AF * w) — the allele-frequency expectation of the raw score."""
    total, seen = 0.0, 0
    for v in score.variants:
        if af_table is not None and population is not None:
            af = af_table.get(v["rsid"], {}).get(population)
        else:
            af = v.get(af_key)
        if af is None:
            continue
        total += 2 * af * v["weight"]
        seen += 1
    return total if seen else None


def audit_score(score: ScoreDefinition, genotype: dict[str, str]) -> ScoreAudit:
    """Variant-level integrity of one score against one genotype."""
    weights = [abs(v["weight"]) for v in score.variants]
    w_total = sum(weights) or 1e-12
    matched = [v for v in score.variants if v["rsid"] in genotype]
    w_cov = sum(abs(v["weight"]) for v in matched)
    missing = sorted((v for v in score.variants if v["rsid"] not in genotype),
                     key=lambda v: -abs(v["weight"]))
    pal = [v for v in score.variants if (v["effect_allele"], v["other_allele"]) in PALINDROMIC]
    eff_n = (w_total ** 2) / sum(w * w for w in weights) if weights else 0.0
    return ScoreAudit(
        pgs_id=score.pgs_id, n_total=len(score.variants), n_matched=len(matched),
        weight_total=round(w_total, 6), weight_covered=round(w_cov, 6),
        weight_coverage=round(w_cov / w_total, 6),
        weight_at_risk=round(w_total - w_cov, 6),
        top1_share=round(max(weights) / w_total, 6) if weights else 0.0,
        effective_n=round(eff_n, 2), palindromic_n=len(pal),
        palindromic_share=round(len(pal) / len(score.variants), 6) if score.variants else 0.0,
        missing_top=[{"rsid": v["rsid"], "weight": v["weight"]} for v in missing[:5]],
    )



@dataclass
class LDAudit:
    """Physical-proximity proxy for linkage disequilibrium within a score.

    Real LD needs a haplotype reference panel. Physical clustering is a weaker but
    always-available proxy: variants tens of kb apart on the same chromosome are
    usually correlated, so treating them as independent contributions overstates
    how many separate pieces of evidence a score carries.
    """
    pgs_id: str
    window_kb: float
    n_variants: int
    clusters: list[dict[str, Any]]
    n_clusters_multi: int
    clustered_weight_share: float
    effective_n_independent: float
    effective_n_ld: float
    duplicate_positions: list[dict[str, Any]] = field(default_factory=list)


def _effective_n(weights: Sequence[float]) -> float:
    w = [abs(x) for x in weights]
    denom = sum(x * x for x in w)
    return (sum(w) ** 2) / denom if denom else 0.0


def ld_audit(score: ScoreDefinition, window_kb: float = 250.0) -> LDAudit:
    """Single-linkage clustering of score variants by physical position."""
    by_chr: dict[str, list[dict[str, Any]]] = {}
    for v in score.variants:
        try:
            pos = int(v["pos"])
        except (TypeError, ValueError):
            continue
        by_chr.setdefault(str(v["chr"]), []).append({**v, "pos_int": pos})

    window = window_kb * 1000.0
    clusters: list[dict[str, Any]] = []
    duplicates: list[dict[str, Any]] = []

    for chrom, variants in by_chr.items():
        variants.sort(key=lambda v: v["pos_int"])
        seen_pos: dict[int, list[str]] = {}
        for v in variants:
            seen_pos.setdefault(v["pos_int"], []).append(v["rsid"])
        for pos, rsids in seen_pos.items():
            if len(rsids) > 1:
                duplicates.append({"chr": chrom, "pos": pos, "rsids": rsids})

        current = [variants[0]]
        for v in variants[1:]:
            if v["pos_int"] - current[-1]["pos_int"] <= window:
                current.append(v)
            else:
                clusters.append(_mk_cluster(chrom, current))
                current = [v]
        clusters.append(_mk_cluster(chrom, current))

    total_w = sum(abs(v["weight"]) for v in score.variants) or 1e-12
    multi = [c for c in clusters if len(c["rsids"]) > 1]
    clustered_w = sum(c["weight_sum"] for c in multi)

    return LDAudit(
        pgs_id=score.pgs_id, window_kb=window_kb, n_variants=len(score.variants),
        clusters=sorted(clusters, key=lambda c: -c["weight_sum"]),
        n_clusters_multi=len(multi),
        clustered_weight_share=round(clustered_w / total_w, 6),
        effective_n_independent=round(_effective_n([v["weight"] for v in score.variants]), 2),
        effective_n_ld=round(_effective_n([c["weight_sum"] for c in clusters]), 2),
        duplicate_positions=duplicates,
    )


def _mk_cluster(chrom: str, members: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "chr": chrom,
        "start": members[0]["pos_int"],
        "end": members[-1]["pos_int"],
        "span_kb": round((members[-1]["pos_int"] - members[0]["pos_int"]) / 1000.0, 1),
        "rsids": [m["rsid"] for m in members],
        "weight_sum": sum(abs(m["weight"]) for m in members),
    }


def integrity_verdict(audit: ScoreAudit, min_weight_coverage: float = 0.90,
                      min_effective_n: float = 10.0,
                      max_palindromic_share: float = 0.10,
                      ld: "LDAudit | None" = None,
                      max_clustered_weight_share: float = 0.30) -> IntegrityVerdict:
    """A score computed on too little of itself is not that score."""
    reasons, warnings = [], []
    if ld is not None and ld.duplicate_positions:
        pairs = "; ".join(
            f"chr{d['chr']}:{d['pos']} ({', '.join(d['rsids'])})" for d in ld.duplicate_positions[:3])
        reasons.append(
            f"Data integrity: {len(ld.duplicate_positions)} genomic position(s) carry more than "
            f"one scored variant ({pairs}). This is either a coordinate error or the same locus "
            f"counted twice, and both corrupt the sum. Verify the scoring file before use.")
    if audit.weight_coverage < min_weight_coverage:
        reasons.append(
            f"Only {audit.weight_coverage:.1%} of this score's total effect weight was "
            f"genotyped (minimum {min_weight_coverage:.0%}). {audit.weight_at_risk:.3f} of "
            f"{audit.weight_total:.3f} weight is missing, so the raw score is biased downward "
            f"by an unknown amount rather than merely noisy.")
    eff = ld.effective_n_ld if ld is not None else audit.effective_n
    if eff < min_effective_n:
        detail = (f"{eff:.1f} after grouping variants that sit within {ld.window_kb:g} kb of each "
                  f"other (was {audit.effective_n:.1f} assuming independence)"
                  if ld is not None else f"{eff:.1f}")
        warnings.append(
            f"Effective number of independent contributions is {detail}. The top variant carries "
            f"{audit.top1_share:.1%} of the weight. A score this concentrated moves sharply with "
            f"a single allele-frequency or correlation difference between populations.")
    if ld is not None and ld.clustered_weight_share > max_clustered_weight_share:
        top = next((c for c in ld.clusters if len(c["rsids"]) > 1), None)
        where = (f" The largest such group is {', '.join(top['rsids'])} on chr{top['chr']} "
                 f"spanning {top['span_kb']:g} kb." if top else "")
        warnings.append(
            f"{ld.clustered_weight_share:.0%} of this score's weight sits in groups of variants "
            f"close enough together to be correlated rather than independent.{where} Correlation "
            f"between such variants differs between populations, so this part of the score does "
            f"not transfer at face value.")
    if audit.palindromic_share > max_palindromic_share:
        warnings.append(
            f"{audit.palindromic_n} of {audit.n_total} variants "
            f"({audit.palindromic_share:.1%}) are strand-ambiguous (A/T or C/G). Effect-allele "
            f"orientation cannot be verified from the genotype file alone.")
    return IntegrityVerdict(passed=not reasons, reasons=reasons, warnings=warnings)


def af_shift(score: ScoreDefinition, af_table: dict[str, dict[str, float]],
             sd: float, population: str = "AFR") -> AFShift | None:
    """Decompose the percentile error into per-variant contributions.

    The reference mean is Sum(2*AF_ref*w). Re-centring on `population` moves it by
    Sum(2*(AF_pop - AF_ref)*w). Dividing by the reference sd gives the shift in sd
    units, which maps directly onto a percentile error.
    """
    per: list[dict[str, Any]] = []
    for v in score.variants:
        af_ref = v.get("af_reference")
        af_pop = af_table.get(v["rsid"], {}).get(population.upper())
        if af_ref is None or af_pop is None:
            continue
        delta = 2 * (af_pop - af_ref) * v["weight"]
        per.append({
            "rsid": v["rsid"], "weight": v["weight"], "af_reference": af_ref,
            "af_population": af_pop, "af_delta": round(af_pop - af_ref, 4),
            "delta_mean": delta,
        })
    if not per:
        return None
    shift_raw = sum(p["delta_mean"] for p in per)
    per.sort(key=lambda p: -abs(p["delta_mean"]))
    return AFShift(
        pgs_id=score.pgs_id, population=population.upper(), n_variants_with_af=len(per),
        coverage=round(len(per) / len(score.variants), 4),
        shift_raw=shift_raw, shift_sd=(shift_raw / sd) if sd else float("nan"),
        per_variant=per,
    )


# ── Loading ───────────────────────────────────────────────────────────────────

def _to_float(value: str | None) -> float | None:
    if value is None:
        return None
    value = value.strip()
    if value == "":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def load_reference_panel(path: Path, pcs: Sequence[str] = DEFAULT_PCS) -> list[Individual]:
    """Load labelled reference individuals with PC coordinates."""
    out: list[Individual] = []
    with Path(path).open(newline="") as handle:
        reader = csv.DictReader(handle)
        missing = [c for c in ("sample_id", "population", *pcs) if c not in (reader.fieldnames or [])]
        if missing:
            raise CalibrationError(f"reference panel missing columns: {missing}")
        for row in reader:
            coords = [_to_float(row.get(p)) for p in pcs]
            if any(c is None for c in coords):
                continue  # a reference individual without coordinates cannot anchor a centroid
            out.append(Individual(row["sample_id"], (row.get("population") or "").strip() or None, coords, None))
    if not out:
        raise CalibrationError(f"no usable reference individuals in {path}")
    return out


def load_query_individuals(path: Path, pcs: Sequence[str] = DEFAULT_PCS) -> list[Individual]:
    out: list[Individual] = []
    with Path(path).open(newline="") as handle:
        for row in csv.DictReader(handle):
            coords = [_to_float(row.get(p)) for p in pcs]
            coords = None if any(c is None for c in coords) else coords
            markers = row.get("n_markers_shared")
            n = int(markers) if (markers or "").strip().isdigit() else None
            out.append(Individual(row["sample_id"], (row.get("population") or "").strip() or None,
                                  coords, n, (row.get("sex") or "").strip() or None))
    return out


def load_prs_results(path: Path) -> list[dict[str, Any]]:
    data = json.loads(Path(path).read_text())
    if isinstance(data, dict):
        data = data.get("scores") or data.get("results") or []
    if not isinstance(data, list):
        raise ValueError(f"unrecognised PRS results structure in {path}")
    return data


# ── Calibration ───────────────────────────────────────────────────────────────

def _distance(pcs: Sequence[float], centroid: Sequence[float]) -> float:
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(pcs, centroid)))


def calibrate(
    panel: Iterable[Individual],
    ref_pop: str,
    k_sd: float = 3.0,
    pcs: Sequence[str] = DEFAULT_PCS,
    min_reference_n: int = 10,
) -> Calibration:
    """Derive the abstention threshold from the spread of the reference population."""
    panel = list(panel)
    ref = [s for s in panel if s.population == ref_pop and s.pcs]
    if not ref:
        raise CalibrationError(
            f"reference population {ref_pop!r} not present in panel; "
            f"available: {sorted({s.population for s in panel if s.population})}"
        )
    if len(ref) < min_reference_n:
        raise CalibrationError(
            f"reference population {ref_pop!r} has n={len(ref)}, below min_reference_n="
            f"{min_reference_n}. A threshold from this few individuals is not defensible."
        )

    centroid = [st.mean(s.pcs[i] for s in ref) for i in range(len(pcs))]
    within = sorted(_distance(s.pcs, centroid) for s in ref)
    mean, sd = st.mean(within), st.stdev(within)
    threshold = mean + k_sd * sd

    others = [s for s in panel if s.population != ref_pop and s.pcs]
    by_pop: dict[str, float] = {}
    for s in others:
        d = _distance(s.pcs, centroid)
        by_pop[s.population or "UNLABELLED"] = min(by_pop.get(s.population or "UNLABELLED", d), d)
    nearest_other = min((_distance(s.pcs, centroid) for s in others), default=None)

    return Calibration(
        reference_population=ref_pop,
        pcs_used=tuple(pcs),
        centroid=[round(c, 6) for c in centroid],
        n=len(ref),
        mean=round(mean, 4),
        sd=round(sd, 4),
        k_sd=k_sd,
        threshold=round(threshold, 4),
        within_max=round(within[-1], 4),
        nearest_other=round(nearest_other, 4) if nearest_other is not None else None,
        other_populations={k: round(v, 4) for k, v in sorted(by_pop.items())},
    )


# ── Decision ──────────────────────────────────────────────────────────────────

def decide(individual: Individual, cal: Calibration, min_markers: int = DEFAULT_MIN_MARKERS) -> Decision:
    """Three outcomes. Marker sufficiency is checked before distance, always."""
    common = dict(
        sample_id=individual.sample_id,
        threshold=cal.threshold,
        n_markers_shared=individual.n_markers_shared,
        declared_population=individual.population,
    )

    n = individual.n_markers_shared
    if n is None or n < min_markers:
        return Decision(
            verdict="REFUSE_UNDETERMINABLE",
            distance=None,
            reason=(
                f"Ancestry could not be determined: {n if n is not None else 'unknown'} markers "
                f"shared with the reference panel, below the minimum of {min_markers} "
                f"(Kosoy et al. 2009). Placement in PC space is not possible, so distance to the "
                f"{cal.reference_population} reference was never computed."
            ),
            remedy=(
                "Genotype more ancestry-informative markers, or supply PC coordinates derived "
                "from a panel that overlaps this individual's markers."
            ),
            **common,
        )

    if not individual.pcs:
        return Decision(
            verdict="REFUSE_UNDETERMINABLE",
            distance=None,
            reason=(
                "No PC coordinates supplied for this individual. Missing coordinates are treated "
                "as unplaceable, never as the origin."
            ),
            remedy="Run claw-ancestry-pca (or ancestry-risk-profiler) and supply the coordinates.",
            **common,
        )

    dist = _distance(individual.pcs, cal.centroid)
    if dist > cal.threshold:
        return Decision(
            verdict="REFUSE_DISTANT",
            distance=round(dist, 4),
            reason=(
                f"Distance to the {cal.reference_population} centroid is {dist:.2f}, beyond the "
                f"threshold of {cal.threshold:.2f} (mean {cal.mean:.2f} + {cal.k_sd:g}x sd "
                f"{cal.sd:.2f} of within-{cal.reference_population} spread). The percentile is "
                f"derived from a {cal.reference_population} reference distribution, so for this "
                f"individual it would be uninterpretable rather than merely imprecise."
            ),
            remedy=(
                f"Use a score with a reference distribution matched to this individual's ancestry, "
                f"or report the raw score with no percentile. Non-genetic risk assessment remains valid."
            ),
            **common,
        )

    return Decision(
        verdict="REPORT",
        distance=round(dist, 4),
        reason=(
            f"Distance to the {cal.reference_population} centroid is {dist:.2f}, within the "
            f"threshold of {cal.threshold:.2f}. The percentile's reference distribution is "
            f"applicable to this individual."
        ),
        remedy="None required. Percentile reported with its reference population disclosed.",
        **common,
    )


def gate_scores(scores: Iterable[dict[str, Any]], decision: Decision, cal: Calibration,
                sex: str | None = None,
                audits: dict[str, ScoreAudit] | None = None,
                integrity: dict[str, IntegrityVerdict] | None = None,
                shifts: dict[str, AFShift] | None = None) -> list[dict[str, Any]]:
    """Retain raw scores always; withhold percentile unless every tier passes.

    Tiers, in order: applicability, score integrity, ancestry transferability,
    reference-population match. The first failure decides, and the reason is kept.
    """
    gated: list[dict[str, Any]] = []
    for s_ in scores:
        pid = s_.get("pgs_id")
        score_pop = s_.get("reference_population")
        reasons: list[str] = []
        warnings: list[str] = []

        app = check_applicability(s_, sex)
        if not app.applicable:
            reasons.append(f"Not applicable: {app.reason}")

        integ = (integrity or {}).get(pid)
        if integ is not None:
            reasons.extend(f"Score integrity: {r}" for r in integ.reasons)
            warnings.extend(integ.warnings)

        if decision.verdict != "REPORT":
            reasons.append(f"Ancestry gate: {decision.verdict}.")

        if score_pop is not None and score_pop != cal.reference_population:
            reasons.append(
                f"Reference mismatch: this score is centred on {score_pop} but the gate was "
                f"calibrated against {cal.reference_population}; transferability to "
                f"{score_pop} was not evaluated.")

        allow = not reasons
        shift = (shifts or {}).get(pid)
        if shift is not None and abs(shift.shift_sd) > 0.5:
            warnings.append(
                f"Re-centring on {shift.population} allele frequencies would move the "
                f"reference mean by {shift.shift_sd:+.2f} sd, i.e. this percentile is off by "
                f"roughly that much before any individual-level error.")

        note = ("Reported against the "
                f"{score_pop or cal.reference_population} reference distribution."
                if allow else " ".join(reasons))
        if warnings:
            note = (note + " | Caveats: " + " ".join(warnings)).strip()

        audit = (audits or {}).get(pid)
        gated.append({
            "pgs_id": pid, "trait": s_.get("trait"), "raw_score": s_.get("raw_score"),
            "variants_used": s_.get("variants_used"), "reference_population": score_pop,
            "percentile": s_.get("percentile") if allow else None,
            "risk_category": s_.get("risk_category") if allow else None,
            "z_score": s_.get("z_score") if allow else None,
            "withheld_reasons": reasons, "caveats": warnings, "note": note,
            "weight_coverage": audit.weight_coverage if audit else None,
            "effective_n": audit.effective_n if audit else None,
            "af_shift_sd": round(shift.shift_sd, 4) if shift else None,
        })
    return gated


# ── Figures ───────────────────────────────────────────────────────────────────

def _write_figures(outdir: Path, panel: list[Individual], cal: Calibration,
                   decisions: list[tuple[Individual, Decision]]) -> list[str]:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.patches import Circle
    except Exception:  # pragma: no cover - graceful degradation
        return []

    figdir = outdir / "figures"
    figdir.mkdir(parents=True, exist_ok=True)
    written: list[str] = []
    colours = {"EUR": "#E8B923", "AFR": "#E07B39", "EAS": "#2E8B7A", "SAS": "#2C5F9E", "AMR": "#7EC8E3"}
    VERDICT_COLOUR = {"REPORT": "#2E8B7A", "REFUSE_DISTANT": "#C1443C", "REFUSE_UNDETERMINABLE": "#6B6B6B"}

    def _base(ax):
        for pop in sorted({s.population for s in panel if s.population}):
            pts = [s.pcs for s in panel if s.population == pop]
            ax.scatter([p[0] for p in pts], [p[1] for p in pts], s=28,
                       c=colours.get(pop, "#999999"), label=pop, alpha=0.75, edgecolors="none")
        ax.scatter([cal.centroid[0]], [cal.centroid[1]], marker="X", s=170, c="black",
                   zorder=5, label=f"{cal.reference_population} centroid")
        ax.add_patch(Circle((cal.centroid[0], cal.centroid[1]), cal.threshold, fill=False,
                            ls="--", lw=1.8, ec="#C1443C", zorder=4))
        ax.set_aspect("equal", adjustable="datalim")
        ax.set_xlabel(f"{cal.pcs_used[0]}")
        ax.set_ylabel(f"{cal.pcs_used[1]}")
        ax.grid(alpha=0.25)

    # Overview
    fig, ax = plt.subplots(figsize=(8, 6.5))
    _base(ax)
    for ind, dec in decisions:
        if ind.pcs:
            ax.scatter([ind.pcs[0]], [ind.pcs[1]], marker="*", s=460,
                       c=VERDICT_COLOUR[dec.verdict], edgecolors="black", linewidths=1.1, zorder=6)
            ax.annotate(ind.sample_id, (ind.pcs[0], ind.pcs[1]), textcoords="offset points",
                        xytext=(9, 7), fontsize=8, weight="bold")
    unplaceable = [i.sample_id for i, d in decisions if not i.pcs]
    title = (f"Abstention gate — threshold {cal.threshold:.2f} "
             f"(dashed circle, radius in {cal.pcs_used[0]}-{cal.pcs_used[-1]} space)")
    if unplaceable:
        title += f"\nUnplaceable, not shown: {', '.join(unplaceable)}"
    ax.set_title(title, fontsize=10)
    ax.legend(fontsize=8, loc="best")
    fig.tight_layout()
    fig.savefig(figdir / "gate_overview.png", dpi=150)
    plt.close(fig)
    written.append("figures/gate_overview.png")

    # Per individual
    for ind, dec in decisions:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5),
                                       gridspec_kw={"width_ratios": [1.25, 1]})
        _base(ax1)
        if ind.pcs:
            ax1.scatter([ind.pcs[0]], [ind.pcs[1]], marker="*", s=520,
                        c=VERDICT_COLOUR[dec.verdict], edgecolors="black", linewidths=1.2, zorder=6)
            ax1.plot([cal.centroid[0], ind.pcs[0]], [cal.centroid[1], ind.pcs[1]],
                     ls=":", c="black", lw=1.4, zorder=5)
            if dec.distance is None:
                # Placeable, but the marker check refused before distance was computed.
                ax1.set_title(f"{ind.sample_id}: coordinates supplied but not used\n"
                              f"(refused at the marker check, before distance)", fontsize=10)
            else:
                ax1.set_title(f"{ind.sample_id}: distance {dec.distance:.2f} vs threshold "
                              f"{cal.threshold:.2f}", fontsize=10)
        else:
            ax1.set_title(f"{ind.sample_id}: cannot be placed in PC space", fontsize=10)
        ax1.legend(fontsize=7, loc="best")

        # Distance ruler
        ax2.axis("off")
        far = max([cal.threshold * 1.6] + ([dec.distance * 1.15] if dec.distance else []))
        ax2.set_xlim(0, far)
        ax2.set_ylim(0, 1)
        ax2.hlines(0.62, 0, far, color="#CCCCCC", lw=7)
        ax2.hlines(0.62, 0, min(cal.threshold, far), color="#2E8B7A", lw=7)
        ax2.vlines(cal.threshold, 0.5, 0.74, color="#C1443C", lw=2.5)
        ax2.text(cal.threshold, 0.79, f"threshold {cal.threshold:.2f}", ha="center",
                 fontsize=9, color="#C1443C", weight="bold")
        ax2.text(0, 0.79, f"{cal.reference_population} centroid", fontsize=8, color="#555555")
        if dec.distance is not None:
            ax2.plot([dec.distance], [0.62], marker="*", ms=22,
                     color=VERDICT_COLOUR[dec.verdict], markeredgecolor="black")
            ax2.text(dec.distance, 0.45, f"{dec.distance:.2f}", ha="center", fontsize=9, weight="bold")
        ax2.text(0, 0.28, f"VERDICT: {dec.verdict}", fontsize=11, weight="bold",
                 color=VERDICT_COLOUR[dec.verdict])
        checks = [
            (f"markers >= {DEFAULT_MIN_MARKERS}",
             (ind.n_markers_shared or 0) >= DEFAULT_MIN_MARKERS),
            ("placeable in PC space", bool(ind.pcs)),
            (f"distance <= {cal.threshold:.2f}",
             dec.distance is not None and dec.distance <= cal.threshold),
        ]
        for i, (label, ok) in enumerate(checks):
            ax2.text(0, 0.18 - i * 0.07, f"{'PASS' if ok else 'FAIL'}  {label}",
                     fontsize=9, family="monospace",
                     color="#2E8B7A" if ok else "#C1443C")
        fig.tight_layout()
        name = f"figures/individual_{ind.sample_id}.png"
        fig.savefig(outdir / name, dpi=150)
        plt.close(fig)
        written.append(name)
    return written


# ── Report ────────────────────────────────────────────────────────────────────

def _write_report(outdir: Path, cal: Calibration, results: list[dict[str, Any]],
                  figures: list[str], args_line: str, min_markers: int) -> None:
    L: list[str] = []
    a = L.append
    a("# PRS Abstention Gate Report\n")
    a(f"**Generated**: {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}  ")
    a(f"**Skill**: prs-abstain v{__version__}  ")
    a(f"**Individuals assessed**: {len(results)}\n")

    a("## Decision criteria\n")
    a("Applied in order. The first failing check decides the verdict.\n")
    a("| # | Check | Rule | Source |")
    a("|---|-------|------|--------|")
    a(f"| 1 | Marker sufficiency | shared markers >= {min_markers} | Kosoy et al. 2009 |")
    a(f"| 2 | Placeable | PC coordinates present, never imputed as 0 | this skill |")
    a(f"| 3 | Transferability | distance to {cal.reference_population} centroid <= {cal.threshold:.2f} | calibrated below |")
    a(f"| 4 | Score match | score reference population == {cal.reference_population} | per-score field |\n")

    a("## Threshold calibration\n")
    a(f"- Reference population: **{cal.reference_population}** (n={cal.n})")
    a(f"- PCs used: {', '.join(cal.pcs_used)}")
    a(f"- Within-{cal.reference_population} distance: mean {cal.mean:.2f}, sd {cal.sd:.2f}, max {cal.within_max:.2f}")
    a(f"- Threshold: **{cal.threshold:.2f}** = mean + {cal.k_sd:g} x sd")
    if cal.nearest_other is not None:
        a(f"- Nearest non-{cal.reference_population} individual: {cal.nearest_other:.2f}")
        a(f"- Empirical separation gap: {cal.within_max:.2f} to {cal.nearest_other:.2f} (no individual falls inside)")
    if cal.other_populations:
        a("\n| Population | Closest individual to reference centroid |")
        a("|------------|------------------------------------------|")
        for pop, d in cal.other_populations.items():
            a(f"| {pop} | {d:.2f} |")
    if cal.threshold_exceeds_nearest_other:
        a(f"\n> **THRESHOLD OVERREACH.** The threshold {cal.threshold:.2f} is larger than the "
          f"distance to the nearest non-{cal.reference_population} individual "
          f"({cal.nearest_other:.2f}). At k={cal.k_sd:g} this rule admits individuals from other "
          f"populations and is no longer an abstention rule. Any `REPORT` verdict below should "
          f"be treated as uncalibrated.\n")
    a("\n> The threshold is a **policy choice placed inside an empirically empty gap**, "
      "not a value derived from theory. Admixed individuals fall between clusters and are "
      "the unhandled case: this rule will classify them by distance alone, which is a "
      "known limitation rather than a validated behaviour.\n")

    a("## Verdicts\n")
    a("| Individual | Verdict | Distance | Threshold | Markers |")
    a("|------------|---------|----------|-----------|---------|")
    for r in results:
        d = r["decision"]
        dist = "n/a" if d["distance"] is None else f"{d['distance']:.2f}"
        a(f"| {d['sample_id']} | **{d['verdict']}** | {dist} | {cal.threshold:.2f} | "
          f"{d['n_markers_shared'] if d['n_markers_shared'] is not None else 'unknown'} |")
    a("")

    for r in results:
        d, scores = r["decision"], r["scores"]
        a(f"### {d['sample_id']} — {d['verdict']}\n")
        a(f"**Why**: {d['reason']}\n")
        a(f"**What would change this**: {d['remedy']}\n")
        if d["verdict"] != "REPORT":
            a(f"> {NOT_REASSURANCE}\n")
        a("| PGS ID | Trait | Raw score | Percentile | Note |")
        a("|--------|-------|-----------|------------|------|")
        for s in scores:
            pct = "**WITHHELD**" if s["percentile"] is None else f"{s['percentile']:.1f}%"
            raw = "n/a" if s["raw_score"] is None else f"{s['raw_score']:.4f}"
            a(f"| {s['pgs_id']} | {s['trait']} | {raw} | {pct} | {s['note']} |")
        a("")
        fig = f"figures/individual_{d['sample_id']}.png"
        if fig in figures:
            a(f"![{d['sample_id']}]({fig})\n")

    if "figures/gate_overview.png" in figures:
        a("## Panel overview\n")
        a("![Gate overview](figures/gate_overview.png)\n")

    a("## Reproducibility\n")
    a(f"```bash\n{args_line}\n```\n")
    a("---\n")
    a(f"*{DISCLAIMER}*")
    (outdir / "report.md").write_text("\n".join(L))



def _verdict_plain(verdict: str) -> str:
    return {
        "REPORT": "Result released",
        "REFUSE_DISTANT": "Result withheld: ancestry is too far from the score's reference group",
        "REFUSE_UNDETERMINABLE": "Result withheld: ancestry could not be established",
    }.get(verdict, verdict)


def _write_clinician_report(outdir: Path, cal: Calibration, results: list[dict[str, Any]],
                            figures: list[str], min_markers: int,
                            lds: dict[str, LDAudit] | None = None) -> None:
    """Plain-language report. No jargon, no matrix algebra, no unexplained acronyms."""
    L: list[str] = []
    a = L.append
    a("# Polygenic score review — summary for clinicians\n")
    a(f"*Generated {datetime.now(timezone.utc).strftime('%d %B %Y, %H:%M UTC')} "
      f"by prs-abstain v{__version__}*\n")

    a("## What this document is\n")
    a("A polygenic score adds up many small genetic effects to place someone on a risk "
      "scale. That placement only means something when it is compared against a group the "
      "person actually resembles. Every score reviewed here was built and calibrated in "
      "people of European ancestry.\n")
    a("This tool checks, for each person and each score, whether the comparison holds. "
      "Where it does not, the percentile is withheld and the reason is given. The "
      "underlying genetic sum is still shown, because that part is always valid.\n")

    a("## What this means\n")
    a("- A **released** result can be read normally, with the usual caveats of any risk score.\n"
      "- A **withheld** result is not a borderline or a low number. It means the comparison "
      "group does not fit this person, so no percentile can be honestly quoted.\n"
      "- A withheld percentile is **not evidence of low risk**. Risk may be high, low or "
      "average; this test simply cannot say which. Clinical assessment, family history and "
      "non-genetic risk factors remain valid and should be used.\n")

    a("## What was checked, and why\n")
    a("| Check | The question it answers | Why it matters |")
    a("|---|---|---|")
    a("| Does the score apply | Is this trait relevant to this person at all? | "
      "Some scores are only meaningful in one sex. |")
    a(f"| Enough of the score measured | Were enough of the score's genetic markers actually "
      f"read? | Missing markers pull the total down and look like lower risk. |")
    a("| Is the score fragile | Does one marker dominate the total? | "
      "A score resting on a few markers swings sharply between populations. |")
    a(f"| Ancestry established | Were at least {min_markers} ancestry markers available? | "
      "Without them the person cannot be compared to any group. |")
    a("| Ancestry close enough | Is this person similar to the group the score was built in? | "
      "Outside that group the percentile stops being interpretable. |\n")

    if lds:
        worst = max(lds.values(), key=lambda l: l.clustered_weight_share)
        if worst.clustered_weight_share > 0.30:
            # Prefer a small, legible cluster for the worked example: a pair explains
            # the idea better than an eight-variant block.
            candidates = [(l, c) for l in lds.values() for c in l.clusters
                          if 2 <= len(c["rsids"]) <= 3]
            if candidates:
                worst, top = max(candidates, key=lambda t: t[1]["weight_sum"])
            else:
                top = next((c for c in worst.clusters if len(c["rsids"]) > 1), None)
            a("## One further caution about how these scores are built\n")
            a("A score is meant to add up many independent pieces of genetic evidence. Some of "
              "the markers in these scores sit in the same stretch of DNA, so they are inherited "
              "together and largely repeat the same information rather than adding new "
              "information.\n")
            if top:
                a(f"In the {worst.pgs_id} score, {', '.join(top['rsids'])} lie within "
                  f"{top['span_kb']:g} thousand DNA letters of each other, so they behave as "
                  f"one signal rather than {len(top['rsids'])} separate ones.\n")
            a(f"That matters here because how strongly such markers travel together differs "
              f"between ancestral groups. A pair that acts as one signal in the group the score "
              f"was built in may act differently in someone else, which is a second reason a "
              f"percentile can mislead even when everything else checks out.\n")
    a("## Results\n")
    a("| Person | Outcome | Scores released | Scores withheld |")
    a("|---|---|---|---|")
    for r in results:
        d, sc = r["decision"], r["scores"]
        rel = sum(1 for x in sc if x["percentile"] is not None)
        a(f"| {d['sample_id']} | {_verdict_plain(d['verdict'])} | {rel} | {len(sc) - rel} |")
    a("")

    for r in results:
        d, sc = r["decision"], r["scores"]
        a(f"### {d['sample_id']}\n")
        a(f"**Outcome:** {_verdict_plain(d['verdict'])}\n")
        a(f"**In plain terms:** {_plain_reason(d, cal)}\n")
        a(f"**What would change this:** {d['remedy']}\n")
        a("| Trait | Percentile | Status |")
        a("|---|---|---|")
        for x in sc:
            if x["percentile"] is not None:
                status = "Released"
                pct = f"{x['percentile']:.0f}th"
            else:
                status = _plain_withheld(x)
                pct = "withheld"
            a(f"| {x['trait']} | {pct} | {status} |")
        a("")
        if any(x["percentile"] is None for x in sc):
            a("> A withheld percentile is **not evidence of low risk**.\n")
        fig = f"figures/individual_{d['sample_id']}.png"
        if fig in figures:
            a(f"![{d['sample_id']}]({fig})\n")

    a("## Limitations you should know before using this\n")
    a("- The ancestry comparison uses a small demonstration reference panel. It is not a "
      "clinical-grade ancestry test.\n"
      "- People of mixed ancestry sit between the reference groups. This tool judges them by "
      "distance alone, which is not a validated approach for them.\n"
      "- A released percentile is a position in a distribution, not a probability of disease "
      "and not a diagnosis.\n"
      "- The tool cannot detect an error in the underlying laboratory data.\n")

    a("---\n")
    a(f"*{DISCLAIMER}*")
    (outdir / "report_clinician.md").write_text("\n".join(L))


def _plain_reason(d: dict[str, Any], cal: Calibration) -> str:
    if d["verdict"] == "REPORT":
        return ("This person closely resembles the group the scores were built in, so the "
                "percentiles can be compared meaningfully.")
    if d["verdict"] == "REFUSE_DISTANT":
        return ("This person's genetic background sits well outside the group these scores "
                "were built in. The scores can still be added up, but the resulting position "
                "on the scale would not mean what it appears to mean.")
    return ("Too few ancestry markers were available to establish which group this person "
            "should be compared against, so no comparison was attempted at all.")


def _plain_withheld(x: dict[str, Any]) -> str:
    reasons = x.get("withheld_reasons") or []
    joined = " ".join(reasons).lower()
    if "not applicable" in joined:
        return "Withheld: trait does not apply to this person"
    if "score integrity" in joined:
        return "Withheld: too little of the score was measured"
    if "ancestry gate" in joined:
        return "Withheld: ancestry comparison does not hold"
    if "reference mismatch" in joined:
        return "Withheld: score built in a different group"
    return "Withheld"


def _write_technical_report(outdir: Path, cal: Calibration, results: list[dict[str, Any]],
                            audits: dict[str, ScoreAudit], shifts: dict[str, AFShift],
                            figures: list[str], args_line: str, min_markers: int,
                            af_population: str | None,
                            lds: dict[str, LDAudit] | None = None) -> None:
    L: list[str] = []
    a = L.append
    a("# PRS abstention gate — technical report\n")
    a(f"**Skill** prs-abstain v{__version__}  |  **Generated** "
      f"{datetime.now(timezone.utc).isoformat(timespec='seconds')}\n")

    a("## 1. Mechanism\n")
    a("For every curated score in this dataset the reference mean satisfies\n")
    a("```\nmean_ref  =  Sum_i  2 * AF_ref(i) * w(i)\n```\n")
    a("to within rounding. The percentile is therefore not merely *derived from* a European "
      "cohort; its centre **is** a European allele-frequency calculation. Re-centring on a "
      "population with frequencies AF_pop moves the mean by\n")
    a("```\ndelta_mean  =  Sum_i  2 * (AF_pop(i) - AF_ref(i)) * w(i)\n"
      "shift_sd    =  delta_mean / sd_ref\n```\n")
    a("This sum decomposes per variant, so the percentile error is attributable to specific "
      "SNPs rather than asserted at the level of the whole score. That decomposition is the "
      "justification for auditing SNP by SNP rather than stopping at a distance threshold.\n")

    a("## 2. Calibration\n")
    a(f"- Reference population: {cal.reference_population} (n={cal.n})")
    a(f"- PCs used: {', '.join(cal.pcs_used)}; distance metric: Euclidean")
    a(f"- Within-reference distance: mean {cal.mean:.4f}, sd {cal.sd:.4f}, max {cal.within_max:.4f}")
    a(f"- Threshold: {cal.threshold:.4f} = mean + {cal.k_sd:g} x sd")
    a(f"- Nearest non-reference individual: {cal.nearest_other}")
    a(f"- Marker minimum: {min_markers} (Kosoy et al. 2009)")
    if cal.threshold_exceeds_nearest_other:
        a("\n**THRESHOLD OVERREACH**: the threshold exceeds the nearest non-reference "
          "individual. Any REPORT verdict below is uncalibrated.")
    a("")

    a("## 3. Per-score integrity\n")
    a("`effective_n` is the inverse Herfindahl index of |weight|, i.e. the number of equally "
      "weighted variants that would produce the same concentration. Low values mark scores "
      "whose percentile is hostage to a handful of allele frequencies.\n")
    a("| PGS | variants | matched | weight coverage | weight at risk | top1 share | effective_n | palindromic |")
    a("|---|---|---|---|---|---|---|---|")
    for pid, au in sorted(audits.items()):
        a(f"| {pid} | {au.n_total} | {au.n_matched} | {au.weight_coverage:.1%} | "
          f"{au.weight_at_risk:.4f} | {au.top1_share:.1%} | {au.effective_n:.1f} | "
          f"{au.palindromic_n} ({au.palindromic_share:.0%}) |")
    a("")

    if shifts:
        a(f"## 4. Allele-frequency re-centring ({af_population})\n")
        a("| PGS | variants with AF | AF coverage | delta_mean | shift (sd) | top driver |")
        a("|---|---|---|---|---|---|")
        for pid, sh in sorted(shifts.items()):
            top = sh.per_variant[0] if sh.per_variant else None
            drv = (f"{top['rsid']} ({top['delta_mean']:+.4f})" if top else "n/a")
            a(f"| {pid} | {sh.n_variants_with_af} | {sh.coverage:.0%} | {sh.shift_raw:+.4f} | "
              f"{sh.shift_sd:+.3f} | {drv} |")
        a("")
        a("### Top per-variant contributions to the shift\n")
        for pid, sh in sorted(shifts.items()):
            a(f"**{pid}** — five largest contributors\n")
            a("| rsID | weight | AF ref | AF pop | dAF | contribution to mean |")
            a("|---|---|---|---|---|---|")
            for v in sh.per_variant[:5]:
                a(f"| {v['rsid']} | {v['weight']:+.4f} | {v['af_reference']:.3f} | "
                  f"{v['af_population']:.3f} | {v['af_delta']:+.3f} | {v['delta_mean']:+.5f} |")
            a("")
    else:
        a("## 4. Allele-frequency re-centring\n")
        a("No population allele-frequency table supplied, so the percentile shift could not "
          "be quantified. This is a gap in the assessment, not evidence that the shift is "
          "small. Supply `--population-af` with gnomAD or 1000 Genomes frequencies.\n")

    if lds:
        a("## 4b. Linkage disequilibrium\n")
        a("Most GWAS associations are carried by tag variants rather than causal ones. A tag "
          "predicts the causal allele only through the correlation (r^2) between them, and that "
          "correlation is a property of the population it was measured in. LD blocks are shorter "
          "in African-ancestry genomes than in European ones, so a tag chosen in a European "
          "cohort commonly tags the causal variant more weakly elsewhere. The effect is "
          "attenuation: the transferred weight overstates the true effect, and the direction of "
          "the error is not random.\n")
        a("Measuring r^2 needs a haplotype reference panel, which this skill does not bundle. "
          "What it can measure without one is physical clustering, which is a lower bound on the "
          "problem: variants close together are correlated, so an effective-n computed as if "
          "they were independent is too high.\n")
        a("| PGS | variants | correlated groups | clustered weight | effective_n (independent) | effective_n (grouped) | duplicate positions |")
        a("|---|---|---|---|---|---|---|")
        for pid, ld in sorted(lds.items()):
            a(f"| {pid} | {ld.n_variants} | {ld.n_clusters_multi} | "
              f"{ld.clustered_weight_share:.0%} | {ld.effective_n_independent:.1f} | "
              f"{ld.effective_n_ld:.1f} | {len(ld.duplicate_positions)} |")
        a("")
        for pid, ld in sorted(lds.items()):
            multi = [c for c in ld.clusters if len(c["rsids"]) > 1][:3]
            if not multi:
                continue
            a(f"**{pid}** — largest correlated groups (window {ld.window_kb:g} kb)\n")
            a("| chr | span (kb) | variants | summed \\|weight\\| |")
            a("|---|---|---|---|")
            for c in multi:
                a(f"| {c['chr']} | {c['span_kb']:g} | {', '.join(c['rsids'])} | {c['weight_sum']:.4f} |")
            a("")
        dups = [(pid, ld) for pid, ld in sorted(lds.items()) if ld.duplicate_positions]
        if dups:
            a("### Data integrity: duplicated genomic positions\n")
            a("Two scored variants at one coordinate is a coordinate error or a locus counted "
              "twice. Either way the sum is wrong, so these scores are blocked rather than "
              "warned about.\n")
            a("| PGS | position | variants |")
            a("|---|---|---|")
            for pid, ld in dups:
                for d in ld.duplicate_positions:
                    a(f"| {pid} | chr{d['chr']}:{d['pos']} | {', '.join(d['rsids'])} |")
            a("")
    a("## 5. Decisions\n")
    a("| Individual | verdict | distance | threshold | markers | released | withheld |")
    a("|---|---|---|---|---|---|---|")
    for r in results:
        d, sc = r["decision"], r["scores"]
        rel = sum(1 for x in sc if x["percentile"] is not None)
        dist = "n/a" if d["distance"] is None else f"{d['distance']:.4f}"
        a(f"| {d['sample_id']} | {d['verdict']} | {dist} | {cal.threshold:.4f} | "
          f"{d['n_markers_shared']} | {rel} | {len(sc)-rel} |")
    a("")
    for r in results:
        d, sc = r["decision"], r["scores"]
        a(f"### {d['sample_id']}\n")
        a(f"{d['reason']}\n")
        a("| PGS | trait | raw | percentile | weight cov | eff_n | shift sd | note |")
        a("|---|---|---|---|---|---|---|---|")
        for x in sc:
            pct = "WITHHELD" if x["percentile"] is None else f"{x['percentile']:.1f}"
            wc = "n/a" if x["weight_coverage"] is None else f"{x['weight_coverage']:.0%}"
            en = "n/a" if x["effective_n"] is None else f"{x['effective_n']:.0f}"
            sh = "n/a" if x["af_shift_sd"] is None else f"{x['af_shift_sd']:+.2f}"
            a(f"| {x['pgs_id']} | {x['trait']} | {x['raw_score']:.4f} | {pct} | {wc} | "
              f"{en} | {sh} | {x['note'][:160]} |")
        a("")

    a("## 6. Known limitations\n")
    for item in KNOWN_LIMITATIONS:
        a(f"- {item}")
    a("")
    a("## 7. Reproducibility\n")
    a(f"```bash\n{args_line}\n```\n")
    a("---\n")
    a(f"*{DISCLAIMER}*")
    (outdir / "report_technical.md").write_text("\n".join(L))


KNOWN_LIMITATIONS = [
    "Ancestry is taken as supplied PC coordinates. The skill does not verify that the "
    "coordinates came from a panel appropriate to the individual, and a projection computed "
    "against an unsuitable panel will be confidently wrong.",
    "Distance is Euclidean in unscaled PC space. If the reference cluster is elongated, "
    "Mahalanobis distance would be the correct metric; Euclidean over-refuses along the "
    "narrow axis and under-refuses along the broad one.",
    "Admixed individuals are the principal unhandled case. A single centroid and a single "
    "radius cannot express partial membership, and the demo panel contains no admixed "
    "samples to calibrate against.",
    "The threshold is a policy choice placed inside an empirically empty gap. Real cohorts "
    "populate that gap and the clean separation will not reproduce.",
    "Allele-frequency re-centring requires an external population AF source. The bundled "
    "demo AF table is synthetic and must never be used for interpretation.",
    "Linkage disequilibrium is approximated by physical distance, not measured. Real r^2 "
    "requires a haplotype reference panel. Physical clustering is a lower bound: it finds "
    "correlated variants that are close together and misses correlated variants that are far "
    "apart, so the reported effective_n remains an overestimate.",
    "LD-driven attenuation of effect sizes is described but not corrected for. Quantifying it "
    "would require per-population r^2 between each tag and its causal variant, which is "
    "unknown for most GWAS loci.",
    "Strand-ambiguous (A/T, C/G) variants are counted and reported but not resolved. "
    "Resolving them requires strand-aware allele frequencies.",
    "Genome build is read from the score header and not verified against the genotype file.",
    "The reference sd is taken from the score's curated distribution and is assumed to hold "
    "in the target population; in general it does not, so shift_sd is a first-order estimate.",
    "Absolute risk is not computed. A percentile is a position in a distribution, not a "
    "probability, and calibration to absolute risk requires population incidence data.",
]



def _write_pdfs(outdir: Path, results: list[dict[str, Any]]) -> list[str]:
    """Render both reports to PDF. Absence of reportlab is not an error."""
    try:
        sys.path.insert(0, str(SKILL_DIR))
        from pdf_report import markdown_to_pdf
    except Exception as exc:  # pragma: no cover
        print(f"PDF rendering unavailable ({exc}); markdown reports were still written.",
              file=sys.stderr)
        return []

    stamp = datetime.now(timezone.utc).strftime("%d %B %Y, %H:%M UTC")
    n = len(results)
    jobs = [
        ("report_clinician", "Polygenic score review",
         f"Summary for clinicians  |  {n} individuals  |  {stamp}"),
        ("report_technical", "PRS abstention gate: technical report",
         f"prs-abstain v{__version__}  |  {n} individuals  |  {stamp}"),
    ]
    written: list[str] = []
    for stem, title, subtitle in jobs:
        md = outdir / f"{stem}.md"
        if not md.exists():
            continue
        try:
            if markdown_to_pdf(md, outdir / f"{stem}.pdf", title, subtitle,
                               figures_dir=outdir / "figures"):
                written.append(f"{stem}.pdf")
        except Exception as exc:  # pragma: no cover
            print(f"Could not render {stem}.pdf: {exc}", file=sys.stderr)
    return written


# ── CLI ───────────────────────────────────────────────────────────────────────

def main(argv: Sequence[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Gate PRS percentiles on ancestry transferability.")
    p.add_argument("--reference-panel", type=Path, help="CSV: sample_id, population, PC1..PCn")
    p.add_argument("--individuals", type=Path, help="CSV: sample_id, PC1..PCn, n_markers_shared")
    p.add_argument("--prs-results", type=Path, help="gwas-prs prs_results.json")
    p.add_argument("--output", type=Path, required=True, help="Output directory")
    p.add_argument("--ref-pop", default="EUR", help="Reference population to gate against")
    p.add_argument("--k-sd", type=float, default=3.0, help="Threshold = mean + k*sd (default 3.0)")
    p.add_argument("--min-markers", type=int, default=DEFAULT_MIN_MARKERS,
                   help=f"Minimum shared markers (default {DEFAULT_MIN_MARKERS}, Kosoy 2009)")
    p.add_argument("--pcs", default=",".join(DEFAULT_PCS), help="Comma-separated PC columns")
    p.add_argument("--scores", type=Path, help="Directory of PGS Catalog scoring files")
    p.add_argument("--genotype", type=Path, help="23andMe/AncestryDNA genotype file")
    p.add_argument("--population-af", type=Path, help="TSV: rsid, population, effect_allele_frequency")
    p.add_argument("--af-population", default="AFR", help="Population column to re-centre on")
    p.add_argument("--min-weight-coverage", type=float, default=0.90,
                   help="Minimum fraction of a score's total |weight| that must be genotyped")
    p.add_argument("--ld-window-kb", type=float, default=250.0,
                   help="Group score variants within this distance as potentially correlated")
    p.add_argument("--min-effective-n", type=float, default=10.0,
                   help="Warn below this effective number of independent contributions")
    p.add_argument("--demo", action="store_true", help="Run on bundled synthetic demo data")
    p.add_argument("--no-figures", action="store_true")
    p.add_argument("--no-pdf", action="store_true", help="Skip PDF rendering of the two reports")
    args = p.parse_args(argv)

    if args.demo:
        args.reference_panel = EXAMPLES / "demo_reference_pcs.csv"
        args.individuals = EXAMPLES / "demo_query_individuals.csv"
        args.prs_results = EXAMPLES / "demo_prs_results.json"
        args.scores = EXAMPLES / "scores"
        args.genotype = EXAMPLES / "demo_genotype.txt"
        args.population_af = EXAMPLES / "demo_population_af.tsv"
    missing = [n for n, v in [("--reference-panel", args.reference_panel),
                              ("--individuals", args.individuals),
                              ("--prs-results", args.prs_results)] if v is None]
    if missing:
        p.error(f"missing required arguments {missing} (or use --demo)")

    pcs = tuple(c.strip() for c in args.pcs.split(",") if c.strip())
    outdir = Path(args.output)
    if outdir.exists() and any(outdir.iterdir()):
        print(f"Warning: {outdir} is not empty; existing files may be overwritten.", file=sys.stderr)
    (outdir / "tables").mkdir(parents=True, exist_ok=True)
    (outdir / "reproducibility").mkdir(parents=True, exist_ok=True)

    try:
        panel = load_reference_panel(args.reference_panel, pcs)
        cal = calibrate(panel, ref_pop=args.ref_pop, k_sd=args.k_sd, pcs=pcs)
    except CalibrationError as exc:
        print(f"Calibration failed: {exc}", file=sys.stderr)
        return 2
    individuals = load_query_individuals(args.individuals, pcs)
    scores = load_prs_results(args.prs_results)

    score_defs = load_score_definitions(args.scores) if args.scores else {}
    genotype = load_genotype(args.genotype) if args.genotype else {}
    af_table = load_population_af(args.population_af) if args.population_af else {}

    audits: dict[str, ScoreAudit] = {}
    lds: dict[str, LDAudit] = {}
    integrity: dict[str, IntegrityVerdict] = {}
    shifts: dict[str, AFShift] = {}
    curated_sd = {s_.get("pgs_id"): s_ for s_ in scores}
    for pid, sdef in score_defs.items():
        lds[pid] = ld_audit(sdef, window_kb=args.ld_window_kb)
        if genotype:
            audits[pid] = audit_score(sdef, genotype)
            integrity[pid] = integrity_verdict(
                audits[pid], min_weight_coverage=args.min_weight_coverage,
                min_effective_n=args.min_effective_n, ld=lds[pid])
        if af_table:
            rec = curated_sd.get(pid, {})
            raw, z = rec.get("raw_score"), rec.get("z_score")
            sd_ref = abs((raw - expected_mean(sdef)) / z) if (raw is not None and z) else None
            sh = af_shift(sdef, af_table, sd=sd_ref or 1.0, population=args.af_population)
            if sh is not None:
                shifts[pid] = sh

    results, pairs = [], []
    for ind in individuals:
        dec = decide(ind, cal, min_markers=args.min_markers)
        gated = gate_scores(scores, dec, cal, sex=ind.sex,
                            audits=audits, integrity=integrity, shifts=shifts)
        results.append({"decision": dec.__dict__, "scores": gated})
        pairs.append((ind, dec))

    figures = [] if args.no_figures else _write_figures(outdir, panel, cal, pairs)

    with (outdir / "tables" / "decisions.csv").open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["sample_id", "verdict", "distance", "threshold", "n_markers_shared", "reason"])
        for r in results:
            d = r["decision"]
            w.writerow([d["sample_id"], d["verdict"], d["distance"], d["threshold"],
                        d["n_markers_shared"], d["reason"]])
    with (outdir / "tables" / "gated_scores.csv").open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["sample_id", "pgs_id", "trait", "raw_score", "percentile", "note"])
        for r in results:
            for s in r["scores"]:
                w.writerow([r["decision"]["sample_id"], s["pgs_id"], s["trait"],
                            s["raw_score"], s["percentile"], s["note"]])

    cmd = "python " + " ".join([str(Path(__file__).name), *(argv or sys.argv[1:])])
    (outdir / "reproducibility" / "commands.sh").write_text(f"#!/usr/bin/env bash\nset -euo pipefail\n{cmd}\n")
    (outdir / "reproducibility" / "environment.yml").write_text(
        "name: prs-abstain\nchannels:\n  - conda-forge\n  - nodefaults\n"
        "dependencies:\n  - python>=3.10\n  - matplotlib\n"
    )
    json.dump(
        {
            "skill": "prs-abstain",
            "version": __version__,
            "generated": datetime.now(timezone.utc).isoformat(),
            "calibration": {
                "reference_population": cal.reference_population,
                "pcs_used": list(cal.pcs_used),
                "centroid": cal.centroid,
                "n": cal.n,
                "mean": cal.mean,
                "sd": cal.sd,
                "k_sd": cal.k_sd,
                "threshold": cal.threshold,
                "within_max": cal.within_max,
                "nearest_other": cal.nearest_other,
                "threshold_exceeds_nearest_other": cal.threshold_exceeds_nearest_other,
                "min_markers": args.min_markers,
            },
            "decisions": results and [{**r["decision"], "scores": r["scores"]} for r in results],
            "figures": figures,
            "disclaimer": DISCLAIMER,
        },
        (outdir / "result.json").open("w"), indent=2,
    )
    _write_report(outdir, cal, results, figures, cmd, args.min_markers)
    _write_clinician_report(outdir, cal, results, figures, args.min_markers, lds)
    _write_technical_report(outdir, cal, results, audits, shifts, figures, cmd,
                            args.min_markers, args.af_population if shifts else None, lds)
    if audits:
        with (outdir / "tables" / "variant_audit.csv").open("w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["pgs_id", "n_total", "n_matched", "weight_coverage", "weight_at_risk",
                        "top1_share", "effective_n", "palindromic_n"])
            for pid, au in sorted(audits.items()):
                w.writerow([pid, au.n_total, au.n_matched, au.weight_coverage,
                            au.weight_at_risk, au.top1_share, au.effective_n, au.palindromic_n])
    if shifts:
        with (outdir / "tables" / "af_shift_per_variant.csv").open("w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["pgs_id", "rsid", "weight", "af_reference", "af_population",
                        "af_delta", "delta_mean"])
            for pid, sh in sorted(shifts.items()):
                for v in sh.per_variant:
                    w.writerow([pid, v["rsid"], v["weight"], v["af_reference"],
                                v["af_population"], v["af_delta"], v["delta_mean"]])

    pdfs = [] if args.no_pdf else _write_pdfs(outdir, results)

    print(f"\nThreshold {cal.threshold:.2f} ({cal.reference_population}, n={cal.n}, "
          f"mean {cal.mean:.2f} + {cal.k_sd:g}sd)")
    if cal.threshold_exceeds_nearest_other:
        print(f"  WARNING: THRESHOLD OVERREACH - threshold {cal.threshold:.2f} exceeds the "
              f"nearest non-{cal.reference_population} individual at {cal.nearest_other:.2f}. "
              f"k={cal.k_sd:g} is too permissive to be an abstention rule.", file=sys.stderr)
    for r in results:
        d = r["decision"]
        dist = "n/a" if d["distance"] is None else f"{d['distance']:.2f}"
        print(f"  {d['sample_id']:<15} {d['verdict']:<24} distance {dist}")
    print(f"\nReports: {outdir / 'report_clinician.md'}")
    print(f"         {outdir / 'report_technical.md'}")
    for f in pdfs:
        print(f"         {outdir / f}")
    print(f"\n{DISCLAIMER}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
