"""Biobank outcome crosswalk (read side): FinnGen R12 endpoints and Pan-UKBB outcomes on
one join key.

The committed table `catalogs/biobank_outcome_crosswalk.tsv.gz` has one row per outcome of
each biobank (2,469 FinnGen endpoints + 7,223 Pan-UKBB outcomes), mapped or not: an
outcome that cannot be harmonized is a row with an `unmapped_reason`, never a dropped id.
It was built offline from pinned source files (the Pan-UKBB phenotype manifest, the
FinnGen R12 manifest, and FinnGen's two published R12 phenotype-mapping files; md5s in
`SOURCES`). This module only reads it.

The join key (`harmonized_outcome_id`) is anchored on a Pan-UKBB outcome:

* a Pan-UKBB `icd10-<X>-<sex>` outcome -> `ICD10:<X>`;
* any other Pan-UKBB outcome that FinnGen published a match for -> `PANUKBB:<outcome_id>`;
* a FinnGen endpoint takes the id of the Pan-UKBB outcome it matches.

FinnGen endpoint crossings, in precedence order (the first that applies wins):

1. `finngen_published_match`: FinnGen's own R12 mapping files pair the endpoint with a
   Pan-UKBB outcome.
2. `icd10_exact`: the endpoint's inclusion definition is ICD-10 codes only and they fall
   in ONE 3-character block that Pan-UKBB serves. The case definitions still differ
   (FinnGen adds control exclusions and register sources).
3. `icd10_broad_only`: only a broader or narrower ICD-10 counterpart exists. Listed, but
   kept OUT of the cross-dataset views, so a broader case definition is never equated with
   a narrower one.

The FinnGen-derived content (endpoint ids and labels, and the matches taken from FinnGen's
mapping files) is used under the FinnGen public-results terms: free download, citation and
acknowledgment required (see `catalogs/README.md`).
"""

from __future__ import annotations

import csv
import gzip
import os
from dataclasses import dataclass
from functools import lru_cache

_HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_CROSSWALK_PATH = os.path.join(_HERE, "catalogs", "biobank_outcome_crosswalk.tsv.gz")

# The exact source files the committed crosswalk was built from (url, md5).
SOURCES = {
    "pan_ukbb_manifest": (
        "https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_release/phenotype_manifest.tsv.bgz",
        "0735065d06ea8ffb098c119ef419c821"),
    "finngen_manifest": (
        "https://storage.googleapis.com/finngen-public-data-r12/summary_stats/finngen_R12_manifest.tsv",
        "3034fcfbf61bb1ba33a0efb3a36ac482"),
    "finngen_ukbb_mapping": (
        "https://storage.googleapis.com/finngen-public-data-r12/meta_analysis/ukbb/"
        "finngen_R12_meta_analysis_mapping_with_definitions.tsv",
        "a92386421e2f177ea3ad5bc5d7bb0e12"),
    "finngen_mvp_ukbb_mapping": (
        "https://storage.googleapis.com/finngen-public-data-r12/meta_analysis/mvp_ukbb/"
        "FinnGen_R12_MVP_UKBB_mapping.tsv",
        "c3928d62dcf834714722e87e40007bc0"),
}

FAMILY_PAN_UKBB = "pan_ukbb"
FAMILY_FINNGEN = "finngen"

BASIS_PUBLISHED = "finngen_published_match"
BASIS_ICD10_EXACT = "icd10_exact"
BASIS_ICD10_CODE = "icd10_code"  # a Pan-UKBB icd10 outcome's own code
BASIS_BROAD_ONLY = "icd10_broad_only"
BASIS_UNMAPPED = "unmapped"

HARMONIZATION_BASIS_LABELS = {
    BASIS_PUBLISHED: "matched by FinnGen's own published phenotype mapping",
    BASIS_ICD10_EXACT: "same single ICD-10 code block (case definitions may still differ)",
    BASIS_ICD10_CODE: "the outcome is an ICD-10 code block",
    BASIS_BROAD_ONLY: "only a broader or narrower ICD-10 counterpart exists; not compared",
    BASIS_UNMAPPED: "no counterpart in the other biobank",
}

# Unmapped reasons (an enum: a new failure mode cannot merge silently into an old bucket).
UNMAPPED_REASONS = {
    "no_published_definition": "FinnGen publishes no mapping or ICD-10 definition for this endpoint",
    "definition_not_icd10_only": "the endpoint is defined by medication, procedure or other criteria, not ICD-10 codes alone",
    "icd10_block_not_served": "its ICD-10 code block has no Pan-UKBB outcome",
    "definition_unparsable": "its ICD-10 definition uses a code prefix that is not one code block",
    "published_match_not_served": "FinnGen's published match is not a Pan-UKBB outcome",
    "no_icd10_code": "a Pan-UKBB disease code (phecode) with no FinnGen published match",
    "measurement_no_counterpart": "a measurement with no FinnGen counterpart",
    "questionnaire_or_medication": "a questionnaire answer or medication-use outcome, not coded to a disease vocabulary",
}

# Which basis wins when several outcomes of one dataset share a harmonized id.
BASIS_PRECEDENCE = {BASIS_PUBLISHED: 0, BASIS_ICD10_EXACT: 1, BASIS_ICD10_CODE: 2}

CROSSWALK_COLUMNS = (
    "dataset_family", "outcome_id", "outcome_label", "harmonized_outcome_id",
    "harmonized_label", "harmonization_basis", "predicate_id",
    "anchor_outcome_id", "broad_candidates", "unmapped_reason", "mapping_source",
    "source_record_id",
)


@dataclass(frozen=True)
class CrosswalkRow:
    dataset_family: str
    outcome_id: str
    outcome_label: str
    harmonized_outcome_id: str | None
    harmonized_label: str | None
    harmonization_basis: str
    predicate_id: str | None
    anchor_outcome_id: str | None
    broad_candidates: tuple[str, ...]
    unmapped_reason: str | None
    mapping_source: str
    source_record_id: str


@lru_cache(maxsize=4)
def load_crosswalk(path: str = DEFAULT_CROSSWALK_PATH) -> dict[tuple[str, str], CrosswalkRow]:
    """{(dataset_family, outcome_id): row} from the committed crosswalk."""
    out: dict[tuple[str, str], CrosswalkRow] = {}
    with gzip.open(path, "rt", encoding="utf-8", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            row = CrosswalkRow(
                dataset_family=r["dataset_family"], outcome_id=r["outcome_id"],
                outcome_label=r["outcome_label"],
                harmonized_outcome_id=r["harmonized_outcome_id"] or None,
                harmonized_label=r["harmonized_label"] or None,
                harmonization_basis=r["harmonization_basis"],
                predicate_id=r["predicate_id"] or None,
                anchor_outcome_id=r["anchor_outcome_id"] or None,
                broad_candidates=tuple(x for x in r["broad_candidates"].split(";") if x),
                unmapped_reason=r["unmapped_reason"] or None,
                mapping_source=r["mapping_source"], source_record_id=r["source_record_id"])
            out[(row.dataset_family, row.outcome_id)] = row
    return out


@lru_cache(maxsize=4)
def harmonized_pan_ukbb_outcomes(path: str = DEFAULT_CROSSWALK_PATH) -> frozenset[str]:
    """Pan-UKBB outcome ids whose harmonized id a FinnGen endpoint shares: the outcomes
    the cross-biobank replication view can compare (the default lookup scope)."""
    cw = load_crosswalk(path)
    fg = {r.harmonized_outcome_id for r in cw.values()
          if r.dataset_family == FAMILY_FINNGEN and r.harmonized_outcome_id}
    return frozenset(r.outcome_id for r in cw.values()
                     if r.dataset_family == FAMILY_PAN_UKBB and r.harmonized_outcome_id in fg)
