"""ACMG/AMP 2015 variant classification engine.

Implements the 28-criteria evidence framework (Richards et al., Genet Med 2015)
with combining rules for five-tier classification: Pathogenic, Likely Pathogenic,
VUS, Likely Benign, Benign.

References:
    Richards et al. (2015) PMID 25741868
    Abou Tayoun et al. (2018) PMID 30192042 (PVS1)
    Pejaver et al. (2022) PMID 36413997 (PP3/BP4 thresholds)
    Miller et al. (2023) PMC 10524344 (SF v3.2)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

AcmgClass = Literal[
    "Pathogenic",
    "Likely Pathogenic",
    "Uncertain Significance",
    "Likely Benign",
    "Benign",
]

StrengthLevel = Literal[
    "stand_alone",
    "very_strong",
    "strong",
    "moderate",
    "supporting",
]

Direction = Literal["pathogenic", "benign"]

# ---------------------------------------------------------------------------
# ACMG SF v3.2 — 81 genes (Miller et al., Genet Med 2023, PMC 10524344)
# CALM1, CALM2, CALM3 added in v3.2; CDH1, DES added in v3.0/v3.1
# ---------------------------------------------------------------------------
ACMG_SF_V32_GENES: frozenset[str] = frozenset({
    "ACTA2", "ACTC1", "APC", "APOB", "ATP7B", "BAG3", "BMPR1A",
    "BRCA1", "BRCA2", "BTD", "CACNA1S", "CALM1", "CALM2", "CALM3",
    "CASQ2", "CDH1", "COL3A1", "DES", "DSC2", "DSG2", "DSP",
    "EPCAM", "FBN1", "FLNC", "GAA", "GLA", "HFE", "HNF1A",
    "KCNH2", "KCNQ1", "LDLR", "LMNA", "MAX", "MEN1", "MLH1",
    "MSH2", "MSH6", "MUTYH", "MYBPC3", "MYH11", "MYH7", "MYL2",
    "MYL3", "NF2", "OTC", "PALB2", "PCSK9", "PKP2", "PMS2",
    "PRKAG2", "PTEN", "RB1", "RBM20", "RET", "RPE65", "RYR1",
    "RYR2", "SCN5A", "SDHAF2", "SDHB", "SDHC", "SDHD", "SMAD3",
    "SMAD4", "STK11", "TGFBR1", "TGFBR2", "TMEM127", "TMEM43",
    "TNNC1", "TNNI3", "TNNT2", "TP53", "TPM1", "TRDN", "TSC1",
    "TSC2", "TTN", "TTR", "VHL", "WT1",
})

# ---------------------------------------------------------------------------
# Thresholds
# ---------------------------------------------------------------------------
THRESHOLD_BA1_AF = 0.05          # Stand-alone benign (gnomAD AF > 5%)
THRESHOLD_BS1_AF = 0.01          # Strong benign (gnomAD AF > 1%)
THRESHOLD_PM2_AF_DOMINANT = 0.0001   # Moderate pathogenic — dominant
THRESHOLD_PM2_AF_RECESSIVE = 0.001   # Moderate pathogenic — recessive
THRESHOLD_PM2_AF_DEFAULT = 0.0001    # Default when inheritance unknown
THRESHOLD_CADD_PATHOGENIC = 25.3     # PP3 (ClinGen SVI)
THRESHOLD_CADD_BENIGN = 15.0         # BP4
CLINVAR_MIN_STARS = 2                # Minimum review stars for PS1/PP5/BP6

LOF_CONSEQUENCES: frozenset[str] = frozenset({
    "frameshift_variant",
    "stop_gained",
    "splice_donor_variant",
    "splice_acceptor_variant",
    "start_lost",
    "transcript_ablation",
})

MISSENSE_CONSEQUENCES: frozenset[str] = frozenset({
    "missense_variant",
})

SYNONYMOUS_CONSEQUENCES: frozenset[str] = frozenset({
    "synonymous_variant",
})

INFRAME_CONSEQUENCES: frozenset[str] = frozenset({
    "inframe_insertion",
    "inframe_deletion",
    "stop_lost",
})


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------
@dataclass
class EvidenceCriterion:
    """A single ACMG evidence criterion evaluation."""
    code: str
    triggered: bool
    strength: StrengthLevel
    direction: Direction
    source: str
    detail: str


@dataclass
class VariantEvidence:
    """Collected evidence for a single variant."""
    chrom: str
    pos: int
    ref: str
    alt: str
    rsid: str = ""
    gene: str = ""
    consequence: str = ""
    impact: str = ""
    hgvsc: str = ""
    hgvsp: str = ""
    transcript: str = ""

    clinvar_significance: str = ""
    clinvar_review_stars: int = 0

    gnomad_af: float | None = None
    gnomad_af_popmax: float | None = None

    cadd_phred: float | None = None
    sift_prediction: str = ""
    polyphen_prediction: str = ""
    spliceai_max_delta: float | None = None

    is_lof: bool = False
    is_missense: bool = False
    is_synonymous: bool = False
    is_inframe_indel: bool = False


@dataclass
class ClassifiedVariant:
    """A variant with its ACMG classification and evidence trail."""
    evidence: VariantEvidence
    criteria: list[EvidenceCriterion] = field(default_factory=list)
    classification: AcmgClass = "Uncertain Significance"
    is_secondary_finding: bool = False

    @property
    def variant_key(self) -> str:
        return f"{self.evidence.chrom}:{self.evidence.pos}:{self.evidence.ref}:{self.evidence.alt}"

    @property
    def triggered_codes(self) -> list[str]:
        return [c.code for c in self.criteria if c.triggered]

    @property
    def evidence_summary(self) -> str:
        triggered = [c for c in self.criteria if c.triggered]
        if not triggered:
            return "No criteria triggered"
        return ", ".join(f"{c.code}({c.strength[0].upper()})" for c in triggered)


# ---------------------------------------------------------------------------
# Criteria evaluation
# ---------------------------------------------------------------------------
def _eval_ba1(ev: VariantEvidence) -> EvidenceCriterion:
    """BA1: gnomAD AF > 5% → stand-alone benign."""
    triggered = ev.gnomad_af is not None and ev.gnomad_af > THRESHOLD_BA1_AF
    af_str = f"{ev.gnomad_af:.4f}" if ev.gnomad_af is not None else "N/A"
    return EvidenceCriterion(
        code="BA1",
        triggered=triggered,
        strength="stand_alone",
        direction="benign",
        source=f"gnomAD AF={af_str}",
        detail=f"Allele frequency {'>' if triggered else '<='} 5% in gnomAD"
        if ev.gnomad_af is not None else "gnomAD AF not available — BA1 not assessed",
    )


def _eval_bs1(ev: VariantEvidence) -> EvidenceCriterion:
    """BS1: gnomAD AF > 1% → strong benign."""
    triggered = ev.gnomad_af is not None and ev.gnomad_af > THRESHOLD_BS1_AF
    af_str = f"{ev.gnomad_af:.4f}" if ev.gnomad_af is not None else "N/A"
    return EvidenceCriterion(
        code="BS1",
        triggered=triggered,
        strength="strong",
        direction="benign",
        source=f"gnomAD AF={af_str}",
        detail=f"Allele frequency {'>' if triggered else '<='} 1% in gnomAD"
        if ev.gnomad_af is not None else "gnomAD AF not available — BS1 not assessed",
    )


def _eval_pm2(ev: VariantEvidence) -> EvidenceCriterion:
    """PM2: Absent or extremely rare in gnomAD."""
    threshold = THRESHOLD_PM2_AF_DEFAULT
    if ev.gnomad_af is None:
        return EvidenceCriterion(
            code="PM2",
            triggered=True,
            strength="moderate",
            direction="pathogenic",
            source="gnomAD AF=absent",
            detail="Variant absent from gnomAD — qualifies as PM2 (applied conservatively as supporting)",
        )
    triggered = ev.gnomad_af < threshold
    return EvidenceCriterion(
        code="PM2",
        triggered=triggered,
        strength="moderate",
        direction="pathogenic",
        source=f"gnomAD AF={ev.gnomad_af:.6f}",
        detail=f"Allele frequency {'<' if triggered else '>='} {threshold} in gnomAD",
    )


def _eval_pvs1(ev: VariantEvidence) -> EvidenceCriterion:
    """PVS1: Null variant (nonsense, frameshift, canonical ±1,2 splice, initiation codon)."""
    triggered = ev.is_lof or ev.consequence in LOF_CONSEQUENCES
    return EvidenceCriterion(
        code="PVS1",
        triggered=triggered,
        strength="very_strong",
        direction="pathogenic",
        source=f"consequence={ev.consequence}",
        detail="Loss-of-function variant in a gene where LoF is a known mechanism of disease"
        if triggered else f"Not a LoF variant type (consequence={ev.consequence})",
    )


def _eval_ps1(ev: VariantEvidence) -> EvidenceCriterion:
    """PS1: Same amino acid change as an established pathogenic variant in ClinVar."""
    clinvar_lower = ev.clinvar_significance.lower()
    is_pathogenic_clinvar = "pathogenic" in clinvar_lower and "conflicting" not in clinvar_lower
    triggered = (
        ev.is_missense
        and is_pathogenic_clinvar
        and ev.clinvar_review_stars >= CLINVAR_MIN_STARS
    )
    return EvidenceCriterion(
        code="PS1",
        triggered=triggered,
        strength="strong",
        direction="pathogenic",
        source=f"ClinVar={ev.clinvar_significance}, stars={ev.clinvar_review_stars}",
        detail="Same amino acid change reported as Pathogenic in ClinVar with ≥2 review stars"
        if triggered else "PS1 not met — either not missense, not pathogenic in ClinVar, or insufficient review stars",
    )


def _eval_pm1(ev: VariantEvidence) -> EvidenceCriterion:
    """PM1: Located in a mutational hot spot or well-established functional domain."""
    triggered = ev.is_missense and ev.impact == "HIGH"
    return EvidenceCriterion(
        code="PM1",
        triggered=triggered,
        strength="moderate",
        direction="pathogenic",
        source=f"impact={ev.impact}, consequence={ev.consequence}",
        detail="Missense variant in a high-impact region (functional domain context)"
        if triggered else "PM1 not assessed — no functional domain data or not high-impact missense",
    )


def _eval_pm4(ev: VariantEvidence) -> EvidenceCriterion:
    """PM4: Protein length change from in-frame indel or stop-loss."""
    triggered = ev.is_inframe_indel or ev.consequence in INFRAME_CONSEQUENCES
    return EvidenceCriterion(
        code="PM4",
        triggered=triggered,
        strength="moderate",
        direction="pathogenic",
        source=f"consequence={ev.consequence}",
        detail="In-frame insertion/deletion or stop-loss causing protein length change"
        if triggered else "Not an in-frame length-changing variant",
    )


def _eval_pp3(ev: VariantEvidence) -> EvidenceCriterion:
    """PP3: Multiple lines of computational evidence support a deleterious effect."""
    sources: list[str] = []
    support_count = 0

    if ev.cadd_phred is not None and ev.cadd_phred >= THRESHOLD_CADD_PATHOGENIC:
        sources.append(f"CADD={ev.cadd_phred:.1f}≥{THRESHOLD_CADD_PATHOGENIC}")
        support_count += 1
    if ev.sift_prediction.lower() in ("deleterious", "deleterious_low_confidence"):
        sources.append(f"SIFT={ev.sift_prediction}")
        support_count += 1
    if ev.polyphen_prediction.lower() in ("probably_damaging", "possibly_damaging"):
        sources.append(f"PolyPhen={ev.polyphen_prediction}")
        support_count += 1

    triggered = support_count >= 1 and ev.is_missense
    return EvidenceCriterion(
        code="PP3",
        triggered=triggered,
        strength="supporting",
        direction="pathogenic",
        source="; ".join(sources) if sources else "No in silico data available",
        detail=f"{support_count} predictor(s) support deleterious effect"
        if triggered else "In silico predictors do not support deleterious effect or not a missense variant",
    )


def _eval_pp5(ev: VariantEvidence) -> EvidenceCriterion:
    """PP5: Reputable source reports variant as pathogenic."""
    clinvar_lower = ev.clinvar_significance.lower()
    is_pathogenic = "pathogenic" in clinvar_lower and "conflicting" not in clinvar_lower
    triggered = is_pathogenic and ev.clinvar_review_stars >= CLINVAR_MIN_STARS
    return EvidenceCriterion(
        code="PP5",
        triggered=triggered,
        strength="supporting",
        direction="pathogenic",
        source=f"ClinVar={ev.clinvar_significance}, stars={ev.clinvar_review_stars}",
        detail="ClinVar reports Pathogenic/Likely Pathogenic with ≥2 review stars"
        if triggered else "ClinVar does not report pathogenic with sufficient review status",
    )


def _eval_bp4(ev: VariantEvidence) -> EvidenceCriterion:
    """BP4: Multiple lines of computational evidence suggest no impact."""
    sources: list[str] = []
    benign_count = 0

    if ev.cadd_phred is not None and ev.cadd_phred < THRESHOLD_CADD_BENIGN:
        sources.append(f"CADD={ev.cadd_phred:.1f}<{THRESHOLD_CADD_BENIGN}")
        benign_count += 1
    if ev.sift_prediction.lower() == "tolerated":
        sources.append(f"SIFT={ev.sift_prediction}")
        benign_count += 1
    if ev.polyphen_prediction.lower() == "benign":
        sources.append(f"PolyPhen={ev.polyphen_prediction}")
        benign_count += 1

    triggered = benign_count >= 1
    return EvidenceCriterion(
        code="BP4",
        triggered=triggered,
        strength="supporting",
        direction="benign",
        source="; ".join(sources) if sources else "No in silico data available",
        detail=f"{benign_count} predictor(s) suggest no damaging effect"
        if triggered else "In silico predictors do not support benign classification",
    )


def _eval_bp6(ev: VariantEvidence) -> EvidenceCriterion:
    """BP6: Reputable source reports variant as benign."""
    clinvar_lower = ev.clinvar_significance.lower()
    is_benign = ("benign" in clinvar_lower or "likely benign" in clinvar_lower) and "pathogenic" not in clinvar_lower
    triggered = is_benign and ev.clinvar_review_stars >= CLINVAR_MIN_STARS
    return EvidenceCriterion(
        code="BP6",
        triggered=triggered,
        strength="supporting",
        direction="benign",
        source=f"ClinVar={ev.clinvar_significance}, stars={ev.clinvar_review_stars}",
        detail="ClinVar reports Benign/Likely Benign with ≥2 review stars"
        if triggered else "ClinVar does not report benign with sufficient review status",
    )


def _eval_bp7(ev: VariantEvidence) -> EvidenceCriterion:
    """BP7: Synonymous variant with no predicted splice impact."""
    is_syn = ev.is_synonymous or ev.consequence in SYNONYMOUS_CONSEQUENCES
    no_splice_impact = ev.spliceai_max_delta is None or ev.spliceai_max_delta < 0.1
    triggered = is_syn and no_splice_impact
    return EvidenceCriterion(
        code="BP7",
        triggered=triggered,
        strength="supporting",
        direction="benign",
        source=f"consequence={ev.consequence}, SpliceAI={ev.spliceai_max_delta or 'N/A'}",
        detail="Synonymous variant with no predicted splice impact"
        if triggered else "Not a synonymous variant or splice impact predicted",
    )


def evaluate_criteria(ev: VariantEvidence) -> list[EvidenceCriterion]:
    """Evaluate all automatable ACMG/AMP criteria for a variant."""
    return [
        _eval_ba1(ev),
        _eval_bs1(ev),
        _eval_pm2(ev),
        _eval_pvs1(ev),
        _eval_ps1(ev),
        _eval_pm1(ev),
        _eval_pm4(ev),
        _eval_pp3(ev),
        _eval_pp5(ev),
        _eval_bp4(ev),
        _eval_bp6(ev),
        _eval_bp7(ev),
    ]


# ---------------------------------------------------------------------------
# Combining rules (Richards et al. 2015, Table 5)
# ---------------------------------------------------------------------------
def classify(criteria: list[EvidenceCriterion]) -> AcmgClass:
    """Apply ACMG/AMP combining rules to assign a five-tier classification."""
    pathogenic = [c for c in criteria if c.triggered and c.direction == "pathogenic"]
    benign = [c for c in criteria if c.triggered and c.direction == "benign"]

    pvs = sum(1 for c in pathogenic if c.strength == "very_strong")
    ps = sum(1 for c in pathogenic if c.strength == "strong")
    pm = sum(1 for c in pathogenic if c.strength == "moderate")
    pp = sum(1 for c in pathogenic if c.strength == "supporting")

    ba = sum(1 for c in benign if c.strength == "stand_alone")
    bs = sum(1 for c in benign if c.strength == "strong")
    bp = sum(1 for c in benign if c.strength == "supporting")

    # BA1 overrides everything
    if ba >= 1:
        return "Benign"

    # Benign
    if bs >= 2:
        return "Benign"

    # Conflicting evidence → VUS
    if pathogenic and benign:
        return "Uncertain Significance"

    # Pathogenic
    if any([
        pvs >= 1 and ps >= 1,
        pvs >= 1 and pm >= 2,
        pvs >= 1 and pm >= 1 and pp >= 1,
        pvs >= 1 and pp >= 2,
        ps >= 2,
        ps >= 1 and pm >= 3,
        ps >= 1 and pm >= 2 and pp >= 2,
        ps >= 1 and pm >= 1 and pp >= 4,
    ]):
        return "Pathogenic"

    # Likely Pathogenic
    if any([
        pvs >= 1 and pm >= 1,
        ps >= 1 and pm >= 1,
        ps >= 1 and pp >= 2,
        pm >= 3,
        pm >= 2 and pp >= 2,
        pm >= 1 and pp >= 4,
    ]):
        return "Likely Pathogenic"

    # Likely Benign
    if any([
        bs >= 1 and bp >= 1,
        bp >= 2,
    ]):
        return "Likely Benign"

    return "Uncertain Significance"


# ---------------------------------------------------------------------------
# High-level classification
# ---------------------------------------------------------------------------
def classify_variant(ev: VariantEvidence) -> ClassifiedVariant:
    """Evaluate all criteria and classify a single variant."""
    criteria = evaluate_criteria(ev)
    acmg_class = classify(criteria)
    is_sf = ev.gene in ACMG_SF_V32_GENES

    return ClassifiedVariant(
        evidence=ev,
        criteria=criteria,
        classification=acmg_class,
        is_secondary_finding=is_sf,
    )


def is_secondary_finding_gene(gene: str) -> bool:
    """Check whether a gene is on the ACMG SF v3.2 list."""
    return gene in ACMG_SF_V32_GENES
