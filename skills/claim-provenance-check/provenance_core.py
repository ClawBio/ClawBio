"""
Provenance core: evidence rows, the citation binder, and the coverage report.

Standard library only, no network, no model. That is deliberate: this is the half of a
report-generation pipeline that must still be trustworthy when the LLM half hallucinates.
Everything here is deterministic and unit tested (see tests/test_claim_provenance_check.py).

TWO IDEAS, IN THE ORDER THEY RUN.

1. THE CITATION BINDER. Every claim carries tags like [trials:2] that must resolve to a real
   row that was actually retrieved. A tag pointing at nothing is an invented citation, which
   is a leading failure mode of LLM-written scientific claims. This check is mechanical, free,
   and runs before any model call, because a model judge costs money and time and cannot be
   trusted to catch a thing that plain arithmetic catches perfectly.

2. THE COVERAGE REPORT. A caller declares which evidence classes it tried to reach for a
   subject. For each class we record whether it was actually reached. A class that could not
   be reached is never silently dropped: it becomes a gap row carrying the named source that
   would close it. Publishing the map including its holes is the whole trust mechanism.

Ported from Sippar's production `reportQuality.ts` and the ClawBio + Nebius Berlin hackathon
demo (`hackathon/clawbio-berlin/provenance/core.py` in Nuru-AI/sippar), trimmed to the citation
binder and coverage report so this skill does one job. The upstream module also carries a
"validity envelope" (bounding a model/simulation's trained region) that is a natural follow-on
skill, not part of this one.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from enum import Enum
from typing import Iterable, Sequence


# ---------------------------------------------------------------------------
# States
# ---------------------------------------------------------------------------

class ClaimState(str, Enum):
    """Evidence-strength language, deliberately not verdict language.

    This skill does not rank, score, or decide whether a claim is TRUE. SUPPORTED and
    CONTESTED describe what the cited evidence does; NO_COVERAGE describes what could not be
    bound to evidence at all, and always ships with the refusal reason attached.
    """

    SUPPORTED = "SUPPORTED"
    CONTESTED = "CONTESTED"
    NO_COVERAGE = "NO COVERAGE"


# ---------------------------------------------------------------------------
# Evidence
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class EvidenceRow:
    """One retrieved fact, from one source, that a claim may cite.

    `tag` is what appears inside a claim, e.g. "trials:2". It must be assigned by the
    retrieval layer, never by the model: a model that can mint its own tags can mint its own
    evidence.
    """

    tag: str
    source_class: str
    text: str
    url: str | None = None
    cohort: str | None = None
    n: int | None = None

    def __post_init__(self) -> None:
        if not self.tag or ":" not in self.tag:
            raise ValueError(f"tag must look like 'class:index', got {self.tag!r}")


@dataclass(frozen=True)
class SourceOutcome:
    """What happened when the caller tried to reach one evidence class.

    A class that returned nothing is recorded with ok=False and a reason, never omitted. When
    ok is False, `route` names what WOULD answer it, because a gap gets a route rather than a
    shrug.
    """

    source_class: str
    ok: bool
    rows: int = 0
    reason: str | None = None
    route: str | None = None
    cost_usd: float | None = None


# ---------------------------------------------------------------------------
# The citation binder
# ---------------------------------------------------------------------------

TAG_RE = re.compile(r"\[([a-z0-9_-]+:\d+)\]", re.IGNORECASE)

#: Shortest run of shared words that counts as copying rather than coincidence. 8 was chosen
#: empirically upstream: shorter fires constantly on ordinary English, longer lets a whole
#: lifted clause through.
VERBATIM_RUN_WORDS = 8


def extract_tags(text: str) -> list[str]:
    """Every citation tag in a claim, in order, lowercased. Duplicates preserved."""
    return [m.group(1).lower() for m in TAG_RE.finditer(text)]


@dataclass
class BinderIssue:
    code: str
    detail: str


def _words(text: str) -> list[str]:
    return re.findall(r"[a-z0-9']+", text.lower())


def check_no_verbatim(claim: str, rows: Sequence[EvidenceRow], run: int = VERBATIM_RUN_WORDS) -> list[BinderIssue]:
    """Flag any run of `run` consecutive words shared between the claim and any source row.

    The claim must be prose ABOUT the source, never a copy of it.
    """
    claim_words = _words(claim)
    if len(claim_words) < run:
        return []
    claim_runs = {
        " ".join(claim_words[i : i + run]) for i in range(len(claim_words) - run + 1)
    }
    issues: list[BinderIssue] = []
    for row in rows:
        row_words = _words(row.text)
        if len(row_words) < run:
            continue
        for i in range(len(row_words) - run + 1):
            candidate = " ".join(row_words[i : i + run])
            if candidate in claim_runs:
                issues.append(
                    BinderIssue("verbatim-run", f"{run}+ words copied from [{row.tag}]: {candidate!r}")
                )
                break
    return issues


def bind_claim(
    claim: str,
    rows: Sequence[EvidenceRow],
    *,
    require_citation: bool = True,
    check_verbatim: bool = True,
) -> list[BinderIssue]:
    """Check one claim against the rows that were actually retrieved.

    Returns an empty list when the claim is bindable. Any issue means the claim does not ship:
    there is no partial credit, because a claim with one invented citation is an invented claim.
    """
    issues: list[BinderIssue] = []
    valid = {r.tag.lower() for r in rows}
    tags = extract_tags(claim)

    if require_citation and not tags:
        issues.append(BinderIssue("citation-floor", "claim carries no citation"))

    for tag in tags:
        if tag not in valid:
            issues.append(BinderIssue("citation-invalid", f"[{tag}] does not match any retrieved row"))

    if check_verbatim:
        cited = [r for r in rows if r.tag.lower() in set(tags)] or list(rows)
        issues.extend(check_no_verbatim(claim, cited))

    return issues


# ---------------------------------------------------------------------------
# The coverage report
# ---------------------------------------------------------------------------

@dataclass
class CoverageReport:
    subject: str
    outcomes: list[SourceOutcome]

    @property
    def reached(self) -> list[str]:
        return [o.source_class for o in self.outcomes if o.ok]

    @property
    def gaps(self) -> list[SourceOutcome]:
        return [o for o in self.outcomes if not o.ok]

    @property
    def total_cost_usd(self) -> float:
        return round(sum(o.cost_usd or 0.0 for o in self.outcomes), 6)

    def as_dict(self) -> dict:
        return {
            "subject": self.subject,
            "reached": self.reached,
            "gaps": [
                {"class": g.source_class, "reason": g.reason, "route": g.route} for g in self.gaps
            ],
            "coverage": f"{len(self.reached)}/{len(self.outcomes)}",
            "costUsd": self.total_cost_usd,
        }

    def render(self) -> str:
        """Plain-text coverage map. This is the headline artifact, holes included."""
        lines = [f"COVERAGE  {self.subject}  ({len(self.reached)}/{len(self.outcomes)} classes)"]
        for o in self.outcomes:
            if o.ok:
                lines.append(f"  [reached] {o.source_class}  ({o.rows} rows)")
            else:
                route = f"  route: {o.route}" if o.route else ""
                lines.append(f"  [NO COVERAGE] {o.source_class}  {o.reason or ''}{route}")
        if self.total_cost_usd:
            lines.append(f"  cost: ${self.total_cost_usd:.4f}")
        return "\n".join(lines)


def state_for_claim(
    claim: str,
    rows: Sequence[EvidenceRow],
    *,
    disagreeing_tags: Iterable[str] = (),
) -> tuple[ClaimState, list[BinderIssue]]:
    """Bind a claim, then state how strong it is.

    A claim that fails the binder is NO_COVERAGE, not a weaker SUPPORTED: an unbindable claim
    has no evidence behind it, which is a coverage fact rather than a confidence one.
    """
    issues = bind_claim(claim, rows)
    if issues:
        return ClaimState.NO_COVERAGE, issues

    tags = set(extract_tags(claim))
    if tags & {t.lower() for t in disagreeing_tags}:
        return ClaimState.CONTESTED, []

    cited_classes = {t.split(":", 1)[0] for t in tags}
    if len(cited_classes) >= 2:
        return ClaimState.SUPPORTED, []
    # One class is not independent corroboration. Say so rather than implying it is.
    return ClaimState.CONTESTED, []
