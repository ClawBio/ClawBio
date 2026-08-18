"""Parser for the legacy SnpEff ``EFF`` annotation field.

Why this module exists
----------------------
The ClawBio hackathon brief warns, verbatim:

    "Do not feed the historical quartet directly to rare-high-impact-variants:
     it does not parse the legacy EFF field."

That is accurate. ``skills/rare-high-impact-variants/rare_high_impact_variants.py``
reads the consequence from ``MC``, ``Consequence`` or ``ANN`` and never looks at
``EFF``, so every record annotated by classic SnpEff arrives with an empty
consequence. The skill that documents the population-frequency blind spot
therefore cannot read the annotation format the teaching data is written in.

This module closes that gap, and it deliberately does NOT collapse a record to a
single impact. Classic ``EFF`` carries one annotation per transcript, and the
disagreement *between* those annotations is signal, not noise: a variant can be
HIGH impact on a minor transcript and MODERATE on the canonical one. Collapsing
to ``max(impact)`` — which is what most pipelines do — silently converts an
annotation choice into an apparent biological fact.

Format
------
Classic SnpEff (pre-``ANN``) writes a comma-separated list of::

    Effect(Effect_Impact|Functional_Class|Codon_Change|Amino_Acid_Change|
           Amino_Acid_Length|Gene_Name|Transcript_BioType|Gene_Coding|
           Transcript_ID|Exon_Rank|Genotype_Number[|ERRORS|WARNINGS])

Eleven pipe-separated fields inside the parentheses, optionally followed by
error and warning columns. Empty fields are common and meaningful (a MODIFIER
annotation has no codon change).

Reference: SnpEff documentation, "EFF field (classic)".
"""

from __future__ import annotations

import re
from dataclasses import dataclass, asdict, field

# Impact severity, ordered. Index doubles as a comparable rank.
IMPACT_ORDER = ("MODIFIER", "LOW", "MODERATE", "HIGH")

_EFF_RE = re.compile(r"^\s*([A-Za-z0-9_+*/'-]+)\s*\((.*)\)\s*$")

# Field names in positional order inside the parentheses.
_FIELDS = (
    "impact",
    "functional_class",
    "codon_change",
    "aa_change",
    "aa_length",
    "gene",
    "transcript_biotype",
    "gene_coding",
    "transcript_id",
    "exon_rank",
    "genotype_number",
    "errors",
    "warnings",
)


@dataclass(frozen=True)
class EffAnnotation:
    """One transcript-level annotation from an ``EFF`` field."""

    effect: str
    impact: str
    functional_class: str = ""
    codon_change: str = ""
    aa_change: str = ""
    aa_length: str = ""
    gene: str = ""
    transcript_biotype: str = ""
    gene_coding: str = ""
    transcript_id: str = ""
    exon_rank: str = ""
    genotype_number: str = ""
    errors: str = ""
    warnings: str = ""
    raw: str = ""

    @property
    def impact_rank(self) -> int:
        """Comparable severity. Unknown impacts sort below MODIFIER."""
        try:
            return IMPACT_ORDER.index(self.impact)
        except ValueError:
            return -1

    def as_dict(self) -> dict:
        return asdict(self)


class EffParseError(ValueError):
    """Raised when an EFF token cannot be parsed at all."""


def _split_annotations(eff: str) -> list[str]:
    """Split on commas that are not inside parentheses.

    A plain ``eff.split(",")`` is wrong: codon and amino-acid fields can
    themselves contain commas in some SnpEff versions.
    """
    out: list[str] = []
    depth = 0
    buf: list[str] = []
    for ch in eff:
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth = max(0, depth - 1)
        if ch == "," and depth == 0:
            out.append("".join(buf))
            buf = []
        else:
            buf.append(ch)
    if buf:
        out.append("".join(buf))
    return [tok for tok in (t.strip() for t in out) if tok]


def parse_annotation(token: str) -> EffAnnotation:
    """Parse a single ``EFFECT(...)`` token."""
    m = _EFF_RE.match(token)
    if not m:
        raise EffParseError(f"not an EFF annotation: {token!r}")
    effect, body = m.group(1), m.group(2)
    parts = body.split("|")
    kwargs = {name: (parts[i].strip() if i < len(parts) else "") for i, name in enumerate(_FIELDS)}
    return EffAnnotation(effect=effect, raw=token, **kwargs)


def parse_eff(eff: str, *, strict: bool = False) -> list[EffAnnotation]:
    """Parse a full ``EFF`` field into one annotation per transcript.

    Order is preserved. With ``strict=False`` (the default) unparseable tokens
    are skipped rather than raising, because a single malformed annotation in a
    real-world VCF should not discard the rest of the record.
    """
    annotations: list[EffAnnotation] = []
    for token in _split_annotations(eff):
        try:
            annotations.append(parse_annotation(token))
        except EffParseError:
            if strict:
                raise
    return annotations


@dataclass
class EffSummary:
    """What a record's transcript annotations collectively say."""

    annotations: list[EffAnnotation] = field(default_factory=list)

    @property
    def genes(self) -> list[str]:
        seen: dict[str, None] = {}
        for a in self.annotations:
            if a.gene:
                seen.setdefault(a.gene, None)
        return list(seen)

    @property
    def max_impact(self) -> str:
        if not self.annotations:
            return ""
        return max(self.annotations, key=lambda a: a.impact_rank).impact

    def annotations_at(self, impact: str) -> list[EffAnnotation]:
        return [a for a in self.annotations if a.impact == impact]

    @property
    def transcripts_with_max_impact(self) -> list[str]:
        return [a.transcript_id for a in self.annotations_at(self.max_impact) if a.transcript_id]

    def impact_for_gene(self, gene: str) -> str:
        """Highest impact among annotations naming ``gene``."""
        hits = [a for a in self.annotations if a.gene == gene]
        if not hits:
            return ""
        return max(hits, key=lambda a: a.impact_rank).impact

    @property
    def is_impact_transcript_dependent(self) -> bool:
        """True when the record's severity depends on which transcript you pick.

        Restricted to annotations for the *same* gene, so that a HIGH call in
        gene A next to a MODIFIER call in neighbouring gene B — which is merely
        two genes overlapping, not a disagreement — does not trigger.
        """
        for gene in self.genes:
            ranks = {a.impact_rank for a in self.annotations if a.gene == gene}
            if len(ranks) > 1:
                return True
        return False

    def coding_transcripts(self) -> list[EffAnnotation]:
        return [a for a in self.annotations if a.gene_coding.upper() == "CODING"]


def summarise(eff: str, *, strict: bool = False) -> EffSummary:
    return EffSummary(annotations=parse_eff(eff, strict=strict))
