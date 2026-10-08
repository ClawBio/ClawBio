"""Fail-closed self-audit for the clinical variant reporter.

Deterministic invariants run after ACMG classification and BEFORE a report is
emitted. If any invariant is violated the variant is hard-abstained (its call is
withheld and replaced with ABSTAIN_LABEL plus machine-readable reason codes)
rather than emitted as a confident classification. Safe uncertainty beats
confident hallucination.

The invariants encode failure classes observed in practice:

  IDENTITY_MISMATCH   — the variant resolved at this coordinate is not the one the
                        caller asserted (by gene / HGVS). This catches the
                        wrong-variant-lookup class where an rsID or coordinate was
                        chosen that does not correspond to the intended variant, and
                        the classification is silently produced for a lookalike.
  CONTRADICTORY_EVIDENCE — mutually exclusive computational criteria both fired
                        (PP3 and BP4). Per ClinGen SVI these are exclusive; a single
                        variant cannot have in-silico evidence both for and against
                        pathogenicity applied simultaneously.
  MISSING_PROVENANCE  — a criterion is triggered but carries no evidence source
                        (fired citing "no data"), so its contribution is unauditable.
  ALLELE_IMBALANCE    — the genotype call is not supported by its own read counts:
                        a heterozygous call with an alternate-allele fraction below
                        HET_MIN_ALT_FRACTION, or fewer than MIN_ALT_READS alternate
                        reads (FORMAT/AD). Clusters of low-fraction indels in one
                        window are a classic misalignment artefact; classifying them
                        as loss-of-function (PVS1) produces confident false positives.
                        Only applied when the VCF carries AD.

These are deterministic and reproducible; they cannot hallucinate. Every new failure
class the adversarial review surfaces should be promoted to an invariant here.
"""
from __future__ import annotations
import re
from dataclasses import dataclass, field

ABSTAIN_LABEL = "Abstained (self-audit)"

# criterion pairs that must never both be triggered (same evidence axis, opposite direction)
_MUTUALLY_EXCLUSIVE = [("PP3", "BP4")]
# phrases an evaluator uses when it has no underlying data
_NO_DATA_MARKERS = ("no in silico data", "no data available", "not available", "")
# read-support floors for a called genotype (germline, diploid)
HET_MIN_ALT_FRACTION = 0.20
HOM_MIN_ALT_FRACTION = 0.80
MIN_ALT_READS = 4


def _allele_balance(record) -> AuditViolation | None:
    """Check that the called genotype is supported by its own AD counts."""
    sample = getattr(record, "sample", None) or {}
    gt = (sample.get("GT") or getattr(record, "genotype", "") or "").replace("|", "/")
    ad = sample.get("AD", "")
    if not gt or not ad or "." in ad:
        return None
    try:
        counts = [int(x) for x in ad.split(",")]
        alleles = [int(a) for a in gt.split("/")]
    except ValueError:
        return None
    idx = getattr(record, "alt_index", 1)
    if idx >= len(counts) or idx not in alleles:
        return None
    total = sum(counts)
    if total <= 0:
        return None
    alt_reads = counts[idx]
    frac = alt_reads / total
    het = len(set(alleles)) > 1
    if alt_reads < MIN_ALT_READS:
        return AuditViolation("ALLELE_IMBALANCE",
                              f"only {alt_reads} of {total} reads support the alternate allele (GT {gt})")
    if het and frac < HET_MIN_ALT_FRACTION:
        return AuditViolation("ALLELE_IMBALANCE",
                              f"heterozygous call with alternate-allele fraction {frac:.0%} "
                              f"({alt_reads}/{total} reads), below {HET_MIN_ALT_FRACTION:.0%}")
    if not het and frac < HOM_MIN_ALT_FRACTION:
        return AuditViolation("ALLELE_IMBALANCE",
                              f"homozygous-alternate call with alternate-allele fraction {frac:.0%} "
                              f"({alt_reads}/{total} reads), below {HOM_MIN_ALT_FRACTION:.0%}")
    return None


@dataclass
class AuditViolation:
    code: str
    detail: str


@dataclass
class AuditResult:
    passed: bool
    violations: list[AuditViolation] = field(default_factory=list)

    def reason_codes(self) -> list[str]:
        return [v.code for v in self.violations]


def _norm_protein(hgvs: str) -> str:
    """Normalise an HGVS protein string to its p.change for comparison.

    'ENSP00000418447.2:p.Ser403Ile' -> 'ser403ile'; 'p.Lys1638Ter' -> 'lys1638ter'.
    Maps '*' to 'ter' and strips parentheses/whitespace so representations compare.
    """
    if not hgvs:
        return ""
    s = hgvs.split(":")[-1].strip()
    s = re.sub(r"^p\.", "", s, flags=re.IGNORECASE)
    s = s.replace("(", "").replace(")", "").replace("*", "Ter").strip()
    return s.lower()


def _norm_coding(hgvs: str) -> str:
    if not hgvs:
        return ""
    s = hgvs.split(":")[-1].strip()
    return re.sub(r"^c\.", "", s, flags=re.IGNORECASE).replace(" ", "").lower()


def expected_from_record(record) -> dict:
    """Derive the asserted identity from a VcfRecord's ID and INFO.

    Recognised INFO keys: GENE, EXPECTED_HGVSP (p.change) and EXPECTED_HGVSC (c.change).
    The ID column contributes an rsID assertion when present.
    """
    if record is None:
        return {}
    info = getattr(record, "info", {}) or {}
    rsid = getattr(record, "id", "") or ""
    return {
        "rsid": rsid if rsid.lower().startswith("rs") else "",
        "gene": info.get("GENE", ""),
        "hgvsp": info.get("EXPECTED_HGVSP", ""),
        "hgvsc": info.get("EXPECTED_HGVSC", ""),
    }


def audit_classified(classified, expected: dict | None = None, record=None) -> AuditResult:
    """Run the fail-closed invariants against a ClassifiedVariant."""
    ev = classified.evidence
    expected = expected or {}
    violations: list[AuditViolation] = []

    # 1. IDENTITY_MISMATCH — asserted identity vs VEP-resolved identity
    exp_gene = (expected.get("gene") or "").strip()
    if exp_gene and ev.gene and exp_gene.upper() != ev.gene.upper():
        violations.append(AuditViolation(
            "IDENTITY_MISMATCH",
            f"asserted gene {exp_gene} but coordinate resolves to {ev.gene}"))
    exp_p = _norm_protein(expected.get("hgvsp", ""))
    got_p = _norm_protein(ev.hgvsp)
    if exp_p and got_p and exp_p != got_p:
        violations.append(AuditViolation(
            "IDENTITY_MISMATCH",
            f"asserted protein change p.{expected['hgvsp'].split('p.')[-1]} "
            f"but coordinate resolves to {ev.hgvsp}"))
    exp_c = _norm_coding(expected.get("hgvsc", ""))
    got_c = _norm_coding(ev.hgvsc)
    if exp_c and got_c and exp_c != got_c:
        violations.append(AuditViolation(
            "IDENTITY_MISMATCH",
            f"asserted coding change c.{exp_c} but coordinate resolves to {ev.hgvsc}"))

    triggered = {c.code: c for c in classified.criteria if c.triggered}

    # 2. CONTRADICTORY_EVIDENCE — mutually exclusive computational criteria co-fired
    for a, b in _MUTUALLY_EXCLUSIVE:
        if a in triggered and b in triggered:
            violations.append(AuditViolation(
                "CONTRADICTORY_EVIDENCE",
                f"{a} and {b} both triggered; they are mutually exclusive (ClinGen SVI)"))

    # 3. MISSING_PROVENANCE — a triggered criterion cites no evidence source
    for code, crit in triggered.items():
        src = (getattr(crit, "source", "") or "").strip().lower()
        if any(marker and marker in src for marker in _NO_DATA_MARKERS if marker) or src == "":
            violations.append(AuditViolation(
                "MISSING_PROVENANCE",
                f"{code} triggered but its evidence source is empty/absent"))

    # 4. ALLELE_IMBALANCE — the call is not supported by its own read counts
    if record is not None:
        imbalance = _allele_balance(record)
        if imbalance:
            violations.append(imbalance)

    return AuditResult(passed=not violations, violations=violations)
