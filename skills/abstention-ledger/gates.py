"""Entitlement gates: machine-checked reasons to withhold a variant from ranking.

Design
------
Two kinds of finding, deliberately kept apart.

**Cohort abstentions** are properties of the dataset, stated once. They are the
reason no output from this file is an answer about a patient, no matter how the
individual records score. Collapsing these into per-variant rows would let a
reader page past them.

**Variant gates** are per-record. Each one is a function of evidence that is
present in, or absent from, the input — never an assertion. A gate that fires
attaches the value that made it fire, so a reader can disagree with the check
rather than with the conclusion.

A variant with no gate firing is ``REVIEWABLE``: it survived the checks we were
able to run. That is a statement about our evidence, not about the variant.
"""

from __future__ import annotations

import pathlib
import sys
from dataclasses import dataclass, field

from legacy_eff import EffSummary

# --------------------------------------------------------------------------
# Reuse the library's own ACMG secondary-findings list rather than retyping it.
#
# skills/clinical-variant-reporter/acmg_engine.py imports only dataclasses and
# typing: no network, no key, no VEP call. Importing it is cheap and it means the
# gene set we screen against is the one the library already ships, not a copy of
# ours that can drift from it.
#
# Failure is explicit. If the import does not resolve we set the list to None and
# the gate reports SF_LIST_UNAVAILABLE instead of quietly passing every record,
# because "we did not check" and "nothing matched" are different claims.
# --------------------------------------------------------------------------
_CVR = pathlib.Path(__file__).resolve().parents[1] / "clinical-variant-reporter"
ACMG_SF_GENES: frozenset[str] | None
ACMG_SF_SOURCE: str

try:
    if str(_CVR) not in sys.path:
        sys.path.append(str(_CVR))
    from acmg_engine import ACMG_SF_V32_GENES as _SF  # type: ignore[import-not-found]

    ACMG_SF_GENES = _SF
    ACMG_SF_SOURCE = (
        "ACMG_SF_V32_GENES from skills/clinical-variant-reporter/acmg_engine.py "
        f"({len(_SF)} genes, ACMG SF v3.2)"
    )
except Exception as _exc:  # noqa: BLE001 - any import failure must be reported, not hidden
    ACMG_SF_GENES = None
    ACMG_SF_SOURCE = f"unavailable: {type(_exc).__name__}: {_exc}"

# --------------------------------------------------------------------------
# Cohort-level abstentions
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class CohortAbstention:
    code: str
    claim_blocked: str
    because: str
    evidence: str


def cohort_abstentions(records: list[dict], *, roles: dict[str, str]) -> list[CohortAbstention]:
    """Facts about the dataset that bound every downstream claim."""
    n = len(records)
    son_states = {r["genotypes"]["SON"] for r in records}
    single_parent = sum(
        1
        for r in records
        if (r["genotypes"]["FATHER"] != "0/0") != (r["genotypes"]["MOTHER"] != "0/0")
    )

    return [
        CohortAbstention(
            code="NO_PHENOTYPE",
            claim_blocked="Any statement linking a variant to this individual's clinical picture.",
            because=(
                "The input carries no phenotype and no HPO terms. Segregation without a "
                "phenotype is arithmetic over genotypes; it cannot be evidence for or against "
                "a condition nobody has described. Gene-to-condition links we could look up "
                "would be untested against this person."
            ),
            evidence="Input columns: no phenotype field, no HPO field. Confirmed by the data page.",
        ),
        CohortAbstention(
            code="UNPHASED",
            claim_blocked="Any statement about two variants being on the same or opposite copies.",
            because=(
                "Parent-of-origin here is Mendelian deduction from genotypes, not molecular "
                "phase. Without phase, the two-hit question about a gene is not merely "
                "unanswered — it is unaskable from this file."
            ),
            evidence=(
                f"All {n} records have son genotype in {sorted(son_states)}; origin is inferred "
                "from which parent carries, and the data page labels the column UNPHASED."
            ),
        ),
        CohortAbstention(
            code="SELECTION_BIAS",
            claim_blocked="Any count, burden or per-gene total presented as describing this pedigree.",
            because=(
                "The file contains only records where the son carries and exactly one parent "
                "carries. Every variant where both parents carry, and every variant the son "
                "does not carry, was removed upstream. Totals computed here describe the "
                "filter, not the family."
            ),
            evidence=f"{single_parent}/{n} records have exactly one carrier parent, by construction.",
        ),
        CohortAbstention(
            code="HISTORICAL_ANNOTATION",
            claim_blocked="Any claim that the supplied effect labels reflect current evidence.",
            because=(
                "The effect annotation is legacy SnpEff against versioned RefSeq transcripts. "
                "Transcript sets and consequence rules have both changed since. The supplied "
                "labels are a historical record of what a tool said, not current annotation."
            ),
            evidence="EFF field is classic SnpEff format with pinned NM_ accessions.",
        ),
        CohortAbstention(
            code="NON_PROBAND_DISCLOSURE",
            claim_blocked=(
                "Any per-variant report treated as being only about the designated proband."
            ),
            because=(
                "The son is the teaching proband, but the file states three other people's "
                "genotypes at every position. Each parent-of-origin label is a statement about "
                "a parent; each shared call is a statement about the sister. Publishing a "
                "ranking of the son publishes a partial genotype of three non-probands who are "
                "not the subject of the analysis."
            ),
            evidence=(
                "Sample roles resolved from genotypes: "
                + ", ".join(f"{role}={sid}" for role, sid in roles.items())
            ),
        ),
    ]


# --------------------------------------------------------------------------
# Per-variant gates
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class GateHit:
    code: str
    evidence: str
    source: str


@dataclass
class VariantVerdict:
    variant_id: str
    genes: list[str]
    supplied_impact: str
    hits: list[GateHit] = field(default_factory=list)

    @property
    def verdict(self) -> str:
        return "WITHHELD" if self.hits else "REVIEWABLE"

    @property
    def codes(self) -> list[str]:
        return [h.code for h in self.hits]


# Gene families where a HIGH-effect call on short reads is unreliable often
# enough that the call itself, not the variant, is the thing in question.
# Each entry carries the reason it is listed; the report prints them.
LOW_COMPLEXITY_FAMILIES: dict[str, str] = {
    "OR": "olfactory receptor family; largest human paralogous family, high cross-mapping",
    "MUC": "mucin; long VNTR domains, assembly and alignment both unstable",
    "PRAMEF": "PRAME family cluster at 1p36, segmental duplication",
    "FLG": "filaggrin repeat domains; copy-number variable, repeat-length artefacts",
    "HRNR": "hornerin; repeat-rich, adjacent to the filaggrin cluster",
    "PCDHA": "clustered protocadherin locus, tandem near-identical cassettes",
    "PCDHB": "clustered protocadherin locus, tandem near-identical cassettes",
    "PCDHG": "clustered protocadherin locus, tandem near-identical cassettes",
    "NBPF": "neuroblastoma breakpoint family; DUF1220 tandem repeats",
    "GOLGA": "golgin family; recurrent segmental duplication",
    "KRTAP": "keratin-associated proteins; clustered, highly similar paralogues",
    "LCE": "late cornified envelope cluster at 1q21, tandem paralogues",
    "SPRR": "small proline-rich cluster at 1q21, tandem paralogues",
    "TBC1D3": "segmental duplication family, variable copy number",
    "HLA": "MHC; extreme polymorphism and reference-allele bias",
    "PSORS1C": "within the 6p21 MHC block; reference bias",
    "CCHCR1": "within the 6p21 MHC block; reference bias",
    "TCF19": "within the 6p21 MHC block; reference bias",
    "POLR3C": "1q21 duplicated segment, paralogous cross-mapping",
}

# GRCh37 coordinates of the extended MHC.
MHC_REGION = ("6", 28_477_797, 33_448_354)


def _low_complexity_reason(genes: list[str], chrom: str, pos: int) -> str | None:
    for gene in genes:
        for prefix, reason in LOW_COMPLEXITY_FAMILIES.items():
            if gene.upper().startswith(prefix):
                return f"{gene}: {reason}"
    c, lo, hi = MHC_REGION
    if chrom == c and lo <= pos <= hi:
        return f"{chrom}:{pos} lies in the extended MHC ({c}:{lo}-{hi}, GRCh37); reference-allele bias"
    return None


def gate_transcript_artifact(summary: EffSummary) -> GateHit | None:
    """Fire when the record's severity depends on which transcript you read.

    This is the gate the legacy-EFF parser exists to make possible. A pipeline
    that reduces EFF to ``max(impact)`` reports these as HIGH with no hint that
    a different, equally annotated transcript of the same gene disagrees.
    """
    if not summary.is_impact_transcript_dependent:
        return None
    top = summary.max_impact
    details = []
    for gene in summary.genes:
        per_tx = [
            (a.transcript_id or "-", a.impact, a.effect)
            for a in summary.annotations
            if a.gene == gene
        ]
        if len({p[1] for p in per_tx}) > 1:
            details.append(
                gene + ": " + "; ".join(f"{tx}={imp}({eff})" for tx, imp, eff in per_tx)
            )
    return GateHit(
        code="TRANSCRIPT_ARTIFACT",
        evidence=(
            f"supplied max impact {top} is not reproduced across transcripts of the same gene — "
            + " | ".join(details)
        ),
        source="legacy EFF field of the input record, parsed per transcript",
    )


def gate_low_complexity(genes: list[str], chrom: str, pos: int) -> GateHit | None:
    reason = _low_complexity_reason(genes, chrom, pos)
    if reason is None:
        return None
    return GateHit(
        code="LOW_COMPLEXITY_LOCUS",
        evidence=reason,
        source="curated gene-family list in gates.LOW_COMPLEXITY_FAMILIES",
    )


def non_proband_carriers(genotypes: dict[str, str]) -> list[str]:
    """Which non-probands this record also discloses a carrier genotype for.

    Deliberately **not** a gate. Every record in a one-carrier-parent file
    discloses at least one parent by construction, so as a veto it would fire on
    100% of records and discriminate nothing — which is how it was first written
    here, and it made the review list trivially empty. The disclosure problem is
    real but it is a property of the dataset, so it lives in
    ``NON_PROBAND_DISCLOSURE`` among the cohort abstentions. Per-record, it is an
    annotation a reader may want, not a reason to withhold.
    """
    return [
        role.lower()
        for role in ("FATHER", "MOTHER", "SISTER")
        if genotypes.get(role, "0/0") not in ("0/0", "./.")
    ]


def gate_secondary_finding(genes: list[str]) -> GateHit | None:
    """Withhold findings in ACMG secondary-findings genes.

    This is the gate that costs us something, and it is the reason the skill has
    the name it does. A variant in one of these genes is the *most* clinically
    actionable category in the whole file — which is exactly why returning it
    here would be wrong.

    ACMG's framework is a policy for clinical sequencing. It presupposes a
    clinical context, a pre-test conversation, and a documented opportunity to
    opt out of secondary findings. This dataset has none of the three: no
    phenotype, no clinical relationship, and no evidence that any family member
    was asked. A CC BY licence permits redistributing the data; it does not
    manufacture consent to be told what the data implies about a named person.

    So actionability does not override provenance. It is the argument *for*
    withholding, not against it.
    """
    if ACMG_SF_GENES is None:
        return GateHit(
            code="SF_LIST_UNAVAILABLE",
            evidence=(
                "could not load the ACMG secondary-findings gene list, so this record was "
                "not screened; absence of a hit here is not evidence of absence"
            ),
            source=ACMG_SF_SOURCE,
        )
    hits = sorted({g for g in genes if g.upper() in ACMG_SF_GENES})
    if not hits:
        return None
    return GateHit(
        code="SECONDARY_FINDING_NO_CONSENT",
        evidence=(
            f"{', '.join(hits)} on the ACMG secondary-findings list; returning a finding here "
            "would be opportunistic screening of a named person with no phenotype, no clinical "
            "relationship and no documented opt-out"
        ),
        source=ACMG_SF_SOURCE,
    )


def evaluate_variant(record: dict, summary: EffSummary) -> VariantVerdict:
    """Run every Tier-1 gate over one record."""
    v = VariantVerdict(
        variant_id=record["variant_id"],
        genes=summary.genes,
        supplied_impact=summary.max_impact,
    )
    for hit in (
        gate_transcript_artifact(summary),
        gate_low_complexity(summary.genes, record["chrom"], int(record["pos"])),
        gate_secondary_finding(summary.genes),
    ):
        if hit is not None:
            v.hits.append(hit)
    return v
