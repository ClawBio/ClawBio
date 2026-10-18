"""Cross-biobank replication: one verdict per harmonized outcome per dataset pair.

Threshold: a harmonized outcome REPLICATES when both datasets' p are
below 0.05 / n_shared and the signs agree, where n_shared is the number of harmonized
outcomes the two datasets SHARE (both catalogs carry the outcome), whatever happened to
the lookup. Counting only the outcomes that were read would let an outage loosen the
threshold: a full FinnGen outage made it 0.05 / 0, i.e. 0.05. 5e-8 is drawn on figures as
a reference only and decides nothing here.

A missing value is never a null, and the three kinds of missing are kept apart:
`not_tested_in_one` (a permanent boundary: variant absent, not tested in that ancestry,
allele mismatch), `not_read_in_one` (a read that failed or did not happen: outage,
missing file, unparsable record, not queried), and `low_confidence_in_one` (tested, but
the source flags the value). Every view shows each apart from `both_null`.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from itertools import combinations
from typing import Any

import snp_phewas_schema as S

ALPHA = 0.05

V_REPLICATES = "replicates"
V_CONCORDANT_ONE = "concordant_one_significant"
V_DISCORDANT = "discordant"
V_SINGLE = "single_dataset_signal"
V_BOTH_NULL = "both_null"
V_DIRECTION_UNKNOWN = "direction_unknown"
V_NOT_TESTED = "not_tested_in_one"
V_NOT_READ = "not_read_in_one"
V_LOW_CONFIDENCE = "low_confidence_in_one"

VERDICT_ORDER = (V_REPLICATES, V_CONCORDANT_ONE, V_DISCORDANT, V_SINGLE, V_DIRECTION_UNKNOWN,
                 V_BOTH_NULL, V_LOW_CONFIDENCE, V_NOT_READ, V_NOT_TESTED)
# The verdicts that say a comparison could not be made, and why.
NO_COMPARISON = (V_NOT_READ, V_LOW_CONFIDENCE, V_NOT_TESTED)
VERDICT_LABELS = {
    V_REPLICATES: "replicates (significant in both, same direction)",
    V_CONCORDANT_ONE: "same direction, significant in one only",
    V_DISCORDANT: "opposite directions, significant in both",
    V_SINGLE: "significant in one only, direction not shared",
    V_DIRECTION_UNKNOWN: "significant in both, direction unavailable in one",
    V_BOTH_NULL: "tested in both, significant in neither",
    V_NOT_TESTED: "not tested in one of the datasets",
    V_NOT_READ: "not read in one of the datasets (an outage, or outside this lookup; not a null)",
    V_LOW_CONFIDENCE: "flagged low confidence in one of the datasets",
}

# A cell can join a verdict when the source tested the outcome and gave a usable value,
# or when it is known to sit above p = 0.05 (an API-censored FinnGen value).
USABLE = frozenset({S.ST_OK, S.ST_API_CENSORED})


@dataclass(frozen=True)
class CellSignal:
    tested: bool
    significant: bool
    sign: int | None  # +1 / -1, None when withheld or censored
    missing: str | None = None  # the no-comparison verdict this cell forces, when untested


def _missing_kind(row: S.OutcomeRow | None) -> str:
    if row is not None and row.status in S.TRANSIENT_STATUSES:
        return V_NOT_READ
    if row is not None and row.status == S.ST_LOW_CONFIDENCE:
        return V_LOW_CONFIDENCE
    return V_NOT_TESTED


def cell_signal(row: S.OutcomeRow | None, alpha: float) -> CellSignal:
    if row is None or row.status not in USABLE:
        return CellSignal(False, False, None, _missing_kind(row))
    if row.status == S.ST_API_CENSORED or row.neglog10_p is None:
        return CellSignal(True, False, None)
    p = 10.0 ** (-row.neglog10_p)
    sign = None if row.beta is None or row.beta == 0 else (1 if row.beta > 0 else -1)
    return CellSignal(True, p < alpha, sign)


def verdict(a: CellSignal, b: CellSignal) -> str:
    if not (a.tested and b.tested):
        kinds = {c.missing for c in (a, b) if not c.tested}
        return next(v for v in NO_COMPARISON if v in kinds)
    if a.significant and b.significant:
        if a.sign is None or b.sign is None:
            return V_DIRECTION_UNKNOWN
        return V_REPLICATES if a.sign == b.sign else V_DISCORDANT
    if a.significant or b.significant:
        sig, other = (a, b) if a.significant else (b, a)
        if sig.sign is not None and other.sign is not None and sig.sign == other.sign:
            return V_CONCORDANT_ONE
        return V_SINGLE
    return V_BOTH_NULL


@dataclass
class ReplicationRow:
    harmonized_outcome_id: str
    harmonized_label: str
    cells: dict[str, S.OutcomeRow | None]
    verdicts: dict[tuple[str, str], str] = field(default_factory=dict)

    @property
    def best_neglog10_p(self) -> float:
        vals = [c.neglog10_p for c in self.cells.values()
                if c is not None and c.status == S.ST_OK and c.neglog10_p is not None]
        return max(vals) if vals else 0.0


@dataclass
class ReplicationTable:
    dataset_ids: list[str]
    rows: list[ReplicationRow]
    n_shared: dict[tuple[str, str], int]
    alpha: dict[tuple[str, str], float]

    @property
    def pairs(self) -> list[tuple[str, str]]:
        return list(self.n_shared)

    def counts(self, pair: tuple[str, str]) -> dict[str, int]:
        out = {v: 0 for v in VERDICT_ORDER}
        for r in self.rows:
            out[r.verdicts[pair]] += 1
        return out

    def to_records(self) -> list[dict[str, Any]]:
        recs = []
        for r in self.rows:
            d: dict[str, Any] = {"harmonized_outcome_id": r.harmonized_outcome_id,
                                 "harmonized_label": r.harmonized_label}
            for ds, c in r.cells.items():
                d[f"{ds}__outcome_id"] = c.outcome_id if c else None
                d[f"{ds}__status"] = c.status if c else "no_counterpart"
                d[f"{ds}__beta"] = c.beta if c else None
                d[f"{ds}__se"] = c.se if c else None
                d[f"{ds}__neglog10_p"] = c.neglog10_p if c else None
            for (a, b), v in r.verdicts.items():
                d[f"verdict__{a}__vs__{b}"] = v
            recs.append(d)
        return recs


def build_replication(result: Any, dataset_ids: list[str] | None = None) -> ReplicationTable:
    """Harmonized outcomes present in at least two of the datasets, with a verdict per pair."""
    ids = dataset_ids or list(result.dataset_ids)
    reps: dict[str, dict[str, S.OutcomeRow]] = {}
    labels: dict[str, str] = {}
    for r in result.rows:
        if r.dataset_id in ids and r.harmonized_outcome_id and r.harmonized_representative:
            reps.setdefault(r.harmonized_outcome_id, {})[r.dataset_id] = r
            labels.setdefault(r.harmonized_outcome_id, r.harmonized_label or r.outcome_label)
    shared = {h: d for h, d in reps.items() if len(d) >= 2}
    pairs = list(combinations(ids, 2))
    n_shared: dict[tuple[str, str], int] = {}
    alpha: dict[tuple[str, str], float] = {}
    for a, b in pairs:
        n = sum(1 for d in shared.values() if a in d and b in d)
        n_shared[(a, b)] = n
        alpha[(a, b)] = ALPHA / n if n else ALPHA
    rows = []
    for h, d in shared.items():
        row = ReplicationRow(h, labels[h], {ds: d.get(ds) for ds in ids})
        for a, b in pairs:
            row.verdicts[(a, b)] = verdict(cell_signal(d.get(a), alpha[(a, b)]),
                                           cell_signal(d.get(b), alpha[(a, b)]))
        rows.append(row)
    rows.sort(key=lambda r: (-r.best_neglog10_p, r.harmonized_outcome_id))
    return ReplicationTable(ids, rows, n_shared, alpha)
