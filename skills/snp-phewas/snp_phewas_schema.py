"""Row schema, the status enum, and the two write-time invariants.

Grain: one row per (dataset_id, dataset_release, outcome_id, variant_id_grch38). EVERY
outcome in a dataset's catalog yields exactly one row whatever happened to it, so
"which outcomes were not tested, and why?" is a query on `status`, never a log line.

Invariants (raised, not `assert`ed, so they survive `python -O`):
  * `assert_partition`: the rows of one dataset cover the catalog's outcome ids EXACTLY,
    compared as SETS (a cardinality match can hide one id swapped for another), with no
    duplicate key.
  * `assert_source_ids_present`: the source's own row ids (`SOURCE_ID_COLUMNS`) are
    populated wherever the status says the source returned a record, and never all-null.
"""

from __future__ import annotations

import csv
from dataclasses import asdict, dataclass, fields
from typing import Iterable, Sequence

# --- status: the excluded-item enum -------------------------------------------------------

ST_OK = "ok"
ST_LOW_CONFIDENCE = "low_confidence"
ST_NOT_TESTED = "not_tested_in_outcome"
ST_ABSENT = "absent_in_dataset"
ST_ALLELE_MISMATCH = "allele_mismatch"
ST_BUILD_UNMAPPABLE = "build_unmappable"
ST_FILE_MISSING = "source_file_missing"
ST_UNREACHABLE = "dataset_unreachable"
ST_NOT_QUERIED = "not_queried"
ST_NOT_IN_CATALOG = "outcome_not_in_catalog"
ST_API_CENSORED = "above_api_reporting_threshold"
ST_RECORD_UNPARSABLE = "record_unparsable"
ST_NOT_READABLE = "source_file_not_readable"

STATUS_LABELS: dict[str, str] = {
    ST_OK: "tested",
    ST_LOW_CONFIDENCE: "tested, flagged low confidence by the source",
    ST_NOT_TESTED: "variant present in the file but not tested for this outcome in this ancestry group",
    ST_ABSENT: "variant not in this dataset's results for this outcome",
    ST_ALLELE_MISMATCH: "a different variant sits at this position (alleles do not match)",
    ST_BUILD_UNMAPPABLE: "variant has no unique position on this dataset's genome build",
    ST_FILE_MISSING: "the source did not serve the file its catalog lists for this outcome",
    ST_RECORD_UNPARSABLE: (
        "the source returned a record at this position that could not be read (wrong number "
        "of fields); this is a read failure, not an absence"
    ),
    ST_UNREACHABLE: "dataset could not be reached; this is an outage, not an absence",
    ST_NOT_READABLE: (
        "the GWAS Catalog holds this study without a position index (or without summary "
        "statistics), so it cannot be read; not a null"
    ),
    ST_NOT_QUERIED: "outcome was not part of this lookup",
    ST_NOT_IN_CATALOG: "the source returned an outcome its published catalog does not list",
    ST_API_CENSORED: (
        "tested, but the source's web API only reports values below p = 0.05, so only "
        "'p at least 0.05' is known; the full value needs a file read"
    ),
}

# Which statuses carry numbers, which are a permanent boundary for this release, and which
# are "not yet" (a re-run can change them). Kept apart so a real boundary is never confused
# with a transient failure.
TESTED_STATUSES = frozenset({ST_OK, ST_LOW_CONFIDENCE})
# A file the catalog lists that the source does not serve is stale, not a boundary: the
# next release (or a re-upload) can serve it.
PERMANENT_STATUSES = frozenset({ST_NOT_TESTED, ST_ABSENT, ST_ALLELE_MISMATCH,
                                ST_BUILD_UNMAPPABLE})
# `source_file_not_readable` (a GWAS Catalog study with no position-indexed file) is "not
# read", like an outage: the Catalog can harmonize the study later. Which of its two reasons
# applies (not yet harmonized / no summary statistics deposited) is in `status_detail`.
TRANSIENT_STATUSES = frozenset({ST_UNREACHABLE, ST_NOT_QUERIED, ST_API_CENSORED,
                                ST_FILE_MISSING, ST_RECORD_UNPARSABLE, ST_NOT_READABLE})
CATALOG_STATUSES = TESTED_STATUSES | PERMANENT_STATUSES | TRANSIENT_STATUSES
ALL_STATUSES = CATALOG_STATUSES | {ST_NOT_IN_CATALOG}
#: Rows outside the phenome-wide denominator besides the permanent statuses.
PHENOME_WIDE_EXCLUDED = frozenset({ST_NOT_QUERIED, ST_NOT_IN_CATALOG, ST_NOT_READABLE})

def display_label(label: str) -> str:
    """An outcome label as shown to a reader. A single all-caps word longer than four letters
    is the source shouting, not an acronym, and is shown in sentence case; the raw label
    stays in every table. Measured over both catalogs (9,629 labels): 5 are all caps, and
    this changes only FinnGen's STROKE (IGF-1, SHBG, GERD, COPD are acronyms and stay)."""
    if label.isalpha() and label.isupper() and len(label) > 4:
        return label.capitalize()
    return label


def phenome_wide_n(rows) -> int:
    """The phenome-wide Bonferroni denominator for one dataset's rows: every outcome this
    lookup asked about, minus those the source can never test for this variant (absent,
    not tested in this ancestry, allele mismatch, no position on the build).

    A value the FinnGen API withholds above p = 0.05 was still tested, and an outcome lost to
    an outage was asked about: counting only rows that carry a value let either shrink the
    denominator and loosen the threshold (the API-only route reported 0.05 / 195 where the
    lookup covered 2,469 endpoints). Outcomes outside the lookup (`not_queried`), source
    results outside the catalog, and studies the source holds with no readable file
    (`source_file_not_readable`: outside what this release can look up by construction, not
    by an outage) are not part of it.
    """
    return sum(1 for r in rows if r.status not in PERMANENT_STATUSES
               and r.status not in PHENOME_WIDE_EXCLUDED)


# Statuses for which the source returned a record, so its native variant id must be set.
SOURCE_RECORD_STATUSES = TESTED_STATUSES | {ST_NOT_TESTED, ST_ALLELE_MISMATCH, ST_API_CENSORED}

# Direction-withheld reasons (the plain-language label rides alongside).
WITHHELD_PALINDROMIC_ORDER = "palindromic_allele_order_mismatch"
# The GWAS Catalog harmonizer could not orient the record to the forward strand (its
# `hm_code` 9 or 14-18); the record keeps its p, never a guessed direction.
WITHHELD_SOURCE_NOT_ORIENTATED = "source_not_orientated"
DIRECTION_WITHHELD_LABELS: dict[str, str] = {
    WITHHELD_PALINDROMIC_ORDER: (
        "no direction shown: the variant is strand-ambiguous (A/T or C/G) and the source "
        "lists its alleles in the opposite order, so the effect allele cannot be confirmed"
    ),
    WITHHELD_SOURCE_NOT_ORIENTATED: (
        "no direction shown: the source could not orient this record to the forward strand"
    ),
}

#: A dataset that released autosomal results only (`contigs: autosomes`, MVP) answers an X /
#: Y / MT variant `absent_in_dataset` with this detail, without reading a file.
AUTOSOMES = frozenset(str(i) for i in range(1, 23))
AUTOSOMES_ONLY_DETAIL = "released autosomal results only"

SOURCE_ID_COLUMNS: tuple[str, ...] = ("outcome_id", "variant_id_native")
GRAIN_COLUMNS: tuple[str, ...] = ("dataset_id", "dataset_release", "outcome_id", "variant_id_grch38")


class PartitionError(AssertionError):
    """A dataset's rows do not account for its outcome catalog exactly."""


class SourceIdMissingError(ValueError):
    """A write is missing the source's own row ids."""


@dataclass
class OutcomeRow:
    dataset_id: str
    dataset_release: str
    outcome_id: str
    outcome_label: str
    outcome_category_native: str
    outcome_kind: str
    n_cases: int | None
    n_controls: int | None
    n_total: int | None
    variant_id_grch38: str
    variant_id_native: str | None
    rs_id: str | None
    effect_allele: str
    other_allele: str
    status: str
    allele_flip_applied: bool = False
    palindromic: bool = False
    direction_withheld_reason: str | None = None
    beta: float | None = None
    se: float | None = None
    neglog10_p: float | None = None
    eaf: float | None = None
    low_confidence: bool = False
    source_route: str | None = None
    status_detail: str | None = None
    harmonized_outcome_id: str | None = None
    harmonized_label: str | None = None
    harmonization_basis: str | None = None
    harmonized_representative: bool | None = None

    @property
    def status_label(self) -> str:
        return STATUS_LABELS.get(self.status, self.status)

    @property
    def p_value(self) -> float | None:
        return None if self.neglog10_p is None else 10.0 ** (-self.neglog10_p)

    @property
    def z(self) -> float | None:
        if self.beta is None or not self.se:
            return None
        return self.beta / self.se

    def key(self) -> tuple[str, str, str, str]:
        return (self.dataset_id, self.dataset_release, self.outcome_id, self.variant_id_grch38)


ROW_COLUMNS: tuple[str, ...] = tuple(f.name for f in fields(OutcomeRow))


def assert_partition(dataset_id: str, rows: Sequence[OutcomeRow],
                     catalog_ids: Iterable[str]) -> None:
    """Rows of one dataset == the catalog, as sets, one row per outcome, known statuses."""
    expected = set(catalog_ids)
    got = [r.outcome_id for r in rows if r.dataset_id == dataset_id]
    stray = [r for r in rows if r.dataset_id != dataset_id]
    if stray:
        raise PartitionError(f"{len(stray)} rows of another dataset passed as {dataset_id}")
    dupes = len(got) - len(set(got))
    missing = expected - set(got)
    extra = set(got) - expected
    bad_status = {r.status for r in rows} - CATALOG_STATUSES
    if dupes or missing or extra or bad_status:
        raise PartitionError(
            f"{dataset_id}: rows do not partition the {len(expected)}-outcome catalog: "
            f"{len(missing)} catalog outcomes have no row (e.g. {sorted(missing)[:3]}), "
            f"{len(extra)} rows name outcomes outside the catalog (e.g. {sorted(extra)[:3]}), "
            f"{dupes} duplicate rows, unknown statuses {sorted(bad_status)}"
        )


def assert_source_ids_present(rows: Sequence[OutcomeRow]) -> None:
    """Write-time guard: source ids are present wherever the source returned a record."""
    if not rows:
        return
    if all(not r.outcome_id for r in rows):
        raise SourceIdMissingError(f"outcome_id is null on all {len(rows)} rows")
    for r in rows:
        if not r.outcome_id:
            raise SourceIdMissingError(f"a {r.dataset_id} row has no outcome_id")
        if r.status in SOURCE_RECORD_STATUSES and not r.variant_id_native:
            raise SourceIdMissingError(
                f"{r.dataset_id}/{r.outcome_id} has status {r.status} but no "
                f"variant_id_native: the source's record id was read and dropped"
            )
    has_record = [r for r in rows if r.status in SOURCE_RECORD_STATUSES]
    if has_record and all(r.variant_id_native is None for r in has_record):
        raise SourceIdMissingError("variant_id_native is null on every row with a source record")


def _cell(v: object) -> str:
    if v is None:
        return ""
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, float):
        return repr(v)
    return str(v)


def write_rows_tsv(rows: Sequence[OutcomeRow], path: str) -> str:
    """Write rows (all statuses) with a plain-language status label column."""
    assert_source_ids_present(rows)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow([*ROW_COLUMNS, "status_label"])
        for r in rows:
            d = asdict(r)
            w.writerow([*(_cell(d[c]) for c in ROW_COLUMNS), r.status_label])
    return path
