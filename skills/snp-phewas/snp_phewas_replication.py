"""Cross-biobank replication: one verdict per harmonized outcome per dataset pair.

Threshold: a harmonized outcome REPLICATES when both datasets' p are
below 0.05 / n_shared and the signs agree, where n_shared is the number of harmonized
outcomes the two datasets SHARE (both catalogs carry the outcome), whatever happened to
the lookup. Counting only the outcomes that were read would let an outage loosen the
threshold: a full FinnGen outage made it 0.05 / 0, i.e. 0.05. The one exclusion is known
from the catalog before any read, so no outage can cause it: a study the source holds with
no readable file (an MVP study the GWAS Catalog has not harmonized) can never be compared and
is outside the denominator, as it is outside phenome_wide_n. 5e-8 is drawn on figures as a
reference only and decides nothing here.

A missing value is never a null, and the three kinds of missing are kept apart:
`not_tested_in_one` (a permanent boundary: variant absent, not tested in that ancestry,
allele mismatch), `not_read_in_one` (a read that failed or did not happen: outage,
missing file, unparsable record, not queried), and `low_confidence_in_one` (tested, but
the source flags the value). Every view shows each apart from `both_null`.

A fourth kind is not a missing value at all: `no_counterpart` (with three or more datasets,
one of a pair has no outcome with this harmonized id in its catalog). It is excluded from
that pair's counts and drawn blank, never hatched (hatched means "not tested").

Three-way (exactly three datasets): for an outcome all three catalogs share, `three_way`
summarizes the three cells at `alpha_3 = 0.05 / n_shared_3`, n_shared_3 = the harmonized
outcomes all three catalogs carry, readable files only (same reasoning as the pair
threshold: an outage cannot loosen it). Fewer outcomes are shared by all three than by any
pair, so alpha_3 is looser than a pair's alpha: a row can replicate in all three while a
pair strip, judged at that pair's stricter alpha, does not. The figure's legend says so.
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
V_NO_COUNTERPART = "no_counterpart"

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
    V_NOT_READ: ("not read in one of the datasets (an outage, outside this lookup, or a study "
                 "the source holds with no readable file; not a null)"),
    V_LOW_CONFIDENCE: "flagged low confidence in one of the datasets",
    V_NO_COUNTERPART: "one of the datasets has no matching outcome",
}

# Three-way summary classes (one column, `three_way`; the pair that replicates when two of
# three do is in `three_way_pair`). The no-comparison classes reuse the pair verdict names.
T_ALL_THREE = "replicates_in_all_three"
T_TWO = "replicates_in_two"
T_DISCORDANT = "discordant"
T_DIRECTION_UNKNOWN = "direction_unknown"
T_ONE = "significant_in_one_only"
T_NULL = "null_in_all_three"
THREE_WAY_ORDER = (T_ALL_THREE, T_TWO, T_DISCORDANT, T_DIRECTION_UNKNOWN, T_ONE, T_NULL,
                   V_LOW_CONFIDENCE, V_NOT_READ, V_NOT_TESTED)
THREE_WAY_LABELS = {
    T_ALL_THREE: "replicates in all three (significant in each, same direction)",
    T_TWO: "replicates in two, tested and not significant in the third",
    T_DISCORDANT: "opposite directions, significant in two",
    T_DIRECTION_UNKNOWN: "significant in two or three, a direction unavailable",
    T_ONE: "significant in one only",
    T_NULL: "tested in all three, significant in none",
    V_LOW_CONFIDENCE: "flagged low confidence in one of the three",
    V_NOT_READ: ("not read in one of the three (an outage, outside this lookup, or a study "
                 "the source holds with no readable file; not a null)"),
    V_NOT_TESTED: "not tested in one of the three",
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


def three_way(a: CellSignal, b: CellSignal, c: CellSignal) -> tuple[str, tuple[int, int] | None]:
    """(class, the indices of the replicating pair when two of three replicate)."""
    cells = (a, b, c)
    if not all(x.tested for x in cells):
        kinds = {x.missing for x in cells if not x.tested}
        return next(v for v in NO_COMPARISON if v in kinds), None
    sig = [i for i, x in enumerate(cells) if x.significant]
    signs = {cells[i].sign for i in sig}
    if len(sig) >= 2 and {1, -1} <= signs:
        return T_DISCORDANT, None
    if len(sig) >= 2 and None in signs:
        return T_DIRECTION_UNKNOWN, None
    if len(sig) == 3:
        return T_ALL_THREE, None
    if len(sig) == 2:
        return T_TWO, (sig[0], sig[1])
    if len(sig) == 1:
        return T_ONE, None
    return T_NULL, None


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
    three_way: str | None = None  # set only with exactly three datasets, all sharing the id
    three_way_pair: tuple[str, str] | None = None

    @property
    def n_present(self) -> int:
        return sum(1 for c in self.cells.values() if c is not None)

    def significant_anywhere(self, alpha: dict[tuple[str, str], float]) -> bool:
        """Significant in at least one dataset at the threshold of any pair it is in."""
        for (a, b), v in self.verdicts.items():
            if v == V_NO_COUNTERPART:
                continue
            for ds in (a, b):
                if cell_signal(self.cells.get(ds), alpha[(a, b)]).significant:
                    return True
        return False

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
    n_shared_3: int | None = None  # outcomes all three datasets share (three datasets only)
    alpha_3: float | None = None

    @property
    def pairs(self) -> list[tuple[str, str]]:
        return list(self.n_shared)

    def counts(self, pair: tuple[str, str]) -> dict[str, int]:
        """Verdict counts over the outcomes BOTH datasets of the pair carry (an outcome with
        no counterpart in one of them is not part of the pair)."""
        out = {v: 0 for v in VERDICT_ORDER}
        for r in self.rows:
            if r.verdicts[pair] != V_NO_COUNTERPART:
                out[r.verdicts[pair]] += 1
        return out

    def counts_three_way(self) -> dict[str, int]:
        out = {v: 0 for v in THREE_WAY_ORDER}
        for r in self.rows:
            if r.three_way is not None:
                out[r.three_way] += 1
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
            if self.n_shared_3 is not None:
                d["three_way"] = r.three_way or V_NO_COUNTERPART
                d["three_way_pair"] = "__vs__".join(r.three_way_pair) if r.three_way_pair else None
                d["three_way_n_shared"] = self.n_shared_3
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

    def counted(d: dict[str, S.OutcomeRow], among) -> bool:
        # A study the source holds with no readable file can never be compared, so it is
        # outside every denominator, as it is outside phenome_wide_n.
        return all(ds in d and d[ds].status != S.ST_NOT_READABLE for ds in among)

    for a, b in pairs:
        n = sum(1 for d in shared.values() if counted(d, (a, b)))
        n_shared[(a, b)] = n
        alpha[(a, b)] = ALPHA / n if n else ALPHA
    n3 = a3 = None
    if len(ids) == 3:
        n3 = sum(1 for d in shared.values() if counted(d, ids))
        a3 = ALPHA / n3 if n3 else ALPHA
    rows = []
    for h, d in shared.items():
        row = ReplicationRow(h, labels[h], {ds: d.get(ds) for ds in ids})
        for a, b in pairs:
            if a not in d or b not in d:
                row.verdicts[(a, b)] = V_NO_COUNTERPART
                continue
            row.verdicts[(a, b)] = verdict(cell_signal(d.get(a), alpha[(a, b)]),
                                           cell_signal(d.get(b), alpha[(a, b)]))
        if a3 is not None and all(ds in d for ds in ids):
            cls, pair_ix = three_way(*(cell_signal(d[ds], a3) for ds in ids))
            row.three_way = cls
            row.three_way_pair = (ids[pair_ix[0]], ids[pair_ix[1]]) if pair_ix else None
        rows.append(row)
    rows.sort(key=lambda r: (-r.best_neglog10_p, r.harmonized_outcome_id))
    return ReplicationTable(ids, rows, n_shared, alpha, n3, a3)


def _verdict_lines(labels: dict[str, str], counts: dict[str, int], named: dict[str, list[str]]
                   ) -> list[str]:
    out = []
    for vkey, k in sorted(counts.items(), key=lambda kv: -kv[1]):
        if not k:
            continue
        # Name the signals; a null bucket of hundreds is a count, not a list.
        names = named.get(vkey, [])
        tail = (f" ({', '.join(names)})" if vkey not in NO_COMPARISON
                and vkey not in (V_BOTH_NULL, T_NULL) and len(names) <= 12 else "")
        out.append(f"- {labels[vkey]}: {k:,}{tail}")
    return out


def replication_lines(table: ReplicationTable, *, heading: str = "###") -> list[str]:
    """The replication section in plain language: one block per dataset pair, and the
    three-way summary when three datasets were compared. Every count comes from the table.
    The report and the command-line output both print this one text, so a third dataset
    cannot be reported by one and silently cut to the first pair by the other."""
    from snp_phewas_catalog import get_dataset, labels_listed

    two = len(table.dataset_ids) == 2
    lines = ["", f"{heading} Replication across " + ("the two biobanks" if two else "biobanks"),
             ""]
    lines.append(
        f"{len(table.rows):,} outcomes are harmonized across "
        f"{'the two datasets' if two else 'at least two of the datasets'} (the same ICD-10 "
        f"block or phecode, FinnGen's own published match, or a reviewed match of the same "
        f"lab measurement). An outcome replicates in a pair when both p-values are below "
        f"0.05 / n, n = the outcomes both datasets of the pair share, with the same "
        f"direction.")
    for pair in table.pairs:
        both = labels_listed(list(pair))
        n_shared = table.n_shared[pair]
        c = table.counts(pair)
        if not two:
            lines += ["", f"**{both}**: {n_shared:,} shared outcomes, threshold "
                          f"0.05 / {n_shared:,}." if n_shared else
                      f"**{both}**: no shared outcomes."]
        named: dict[str, list[str]] = {}
        for r in table.rows:
            named.setdefault(r.verdicts[pair], []).append(r.harmonized_label)
        lines += _verdict_lines(VERDICT_LABELS, c, named)
        no_cmp = sum(k for vkey, k in c.items() if vkey in NO_COMPARISON)
        if no_cmp:
            lines.append(f"{no_cmp:,} harmonized outcomes could not be compared; they are "
                         f"listed apart from the nulls, never counted as nulls.")
    if table.n_shared_3 is not None:
        n3 = table.n_shared_3
        short = [get_dataset(d).short_label for d in table.dataset_ids]
        names = ", ".join(short[:-1]) + " and " + short[-1]
        three_families = len({get_dataset(d).family for d in table.dataset_ids}) == 3
        head = "All three biobanks" if three_families else "All three datasets"
        lines += ["", f"**{head}** ({names}): {n3:,} outcomes with a readable file are shared "
                      f"by all three; each is compared at 0.05 / {n3:,}." if n3 else
                  f"**{head}**: no outcome is shared by all three."]
        if n3:
            named3: dict[str, list[str]] = {}
            for r in table.rows:
                if r.three_way is None:
                    continue
                lab = r.harmonized_label
                if r.three_way == T_TWO and r.three_way_pair:
                    lab += " (" + " and ".join(get_dataset(d).short_label
                                               for d in r.three_way_pair) + ")"
                named3.setdefault(r.three_way, []).append(lab)
            c3 = table.counts_three_way()
            lines += _verdict_lines(THREE_WAY_LABELS, {k: c3[k] for k in THREE_WAY_ORDER},
                                    named3)
    return lines
