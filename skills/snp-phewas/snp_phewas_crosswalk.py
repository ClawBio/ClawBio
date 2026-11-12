"""Biobank outcome crosswalk (read side): FinnGen R12 endpoints, Pan-UKBB outcomes and
Million Veteran Program studies on one join key.

The committed table `catalogs/biobank_outcome_crosswalk.tsv.gz` has one row per outcome of
each biobank (2,469 FinnGen endpoints + 7,223 Pan-UKBB outcomes + 7,537 MVP GWAS Catalog
studies), mapped or not: an outcome that cannot be harmonized is a row with an
`unmapped_reason`, never a dropped id. It was built offline from pinned source files (the
Pan-UKBB phenotype manifest, the FinnGen R12 manifest, FinnGen's two published R12
phenotype-mapping files, the MVP outcome catalog and the reviewed MVP measurement table;
md5s in `SOURCES`). This module only reads it.

The join key (`harmonized_outcome_id`) is anchored on a Pan-UKBB outcome, and a row carries
exactly ONE key, compared with every dataset sharing it:

* a Pan-UKBB `icd10-<X>-<sex>` outcome -> `ICD10:<X>`;
* a Pan-UKBB `phecode-<code>-both_sexes` outcome -> `PHECODE:<code>` (normalized code);
* any other Pan-UKBB outcome that FinnGen published a match for, or that an MVP
  measurement matches by the reviewed name table -> `PANUKBB:<outcome_id>`;
* a FinnGen endpoint takes the id of the Pan-UKBB outcome it matches;
* an MVP phecode study takes `PHECODE:<code>` when Pan-UKBB serves that code for both
  sexes (the two phecode codings are the same 1.2 family: a natural join on the code); an
  MVP measurement takes the reviewed match's id. The mean (or single) value is the one
  compared; the maximum / minimum keep the id with basis `summary_statistic_variant` and are
  never compared.

FinnGen endpoint crossings, in precedence order (the first that applies wins):

1. `finngen_published_match`: FinnGen's own R12 mapping files pair the endpoint with a
   Pan-UKBB outcome.
2. `icd10_exact`: the endpoint's inclusion definition is ICD-10 codes only and they fall
   in ONE 3-character block that Pan-UKBB serves. The case definitions still differ
   (FinnGen adds control exclusions and register sources).
3. `icd10_broad_only`: only a broader or narrower ICD-10 counterpart exists. Listed, but
   kept OUT of the cross-dataset views, so a broader case definition is never equated with
   a narrower one.

A FinnGen endpoint's MVP link from FinnGen's own FinnGen/MVP/UKBB file that the one-key
rule cannot compare is kept in `secondary_match` with a `secondary_reason`, never dropped.

The FinnGen-derived content (endpoint ids and labels, and the matches taken from FinnGen's
mapping files) is used under the FinnGen public-results terms: free download, citation and
acknowledgment required (see `catalogs/README.md`).
"""

from __future__ import annotations

import csv
import gzip
import os
import re
from dataclasses import dataclass
from functools import lru_cache

_HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_CROSSWALK_PATH = os.path.join(_HERE, "catalogs", "biobank_outcome_crosswalk.tsv.gz")

# The exact source files the committed crosswalk was built from (url or vendored path, md5).
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
    # Vendored: the MVP outcome catalog of the 7,537 GWAS Catalog studies (rebuilt by
    # snp_phewas_build_catalogs.py from its own pinned snapshots), and the reviewed MVP
    # lab / vital -> Pan-UKBB field table, one row per MVP measurement (76), matched or not.
    # The crosswalk build reads the table's data rows and skips its `#` comment lines.
    "mvp_outcome_catalog": ("catalogs/mvp_outcomes.tsv.gz", "e9f131d6461b87e79a3d857f7caad568"),
    "mvp_measurement_matches": ("catalogs/mvp_measurement_matches.tsv",
                                "c328ed9398bb1ec65978ae40162ddb44"),
}

FAMILY_PAN_UKBB = "pan_ukbb"
FAMILY_FINNGEN = "finngen"
FAMILY_MVP = "mvp"
FAMILIES = (FAMILY_FINNGEN, FAMILY_PAN_UKBB, FAMILY_MVP)

BASIS_PUBLISHED = "finngen_published_match"
BASIS_ICD10_EXACT = "icd10_exact"
BASIS_ICD10_CODE = "icd10_code"  # a Pan-UKBB icd10 outcome's own code
BASIS_PHECODE_CODE = "phecode_code"  # a Pan-UKBB both-sexes phecode outcome's own code
BASIS_PHECODE_SAME = "phecode_same_code"  # an MVP phecode study on the same code
BASIS_CURATED = "curated_name_match"  # an MVP measurement in the reviewed name table
BASIS_SUMMARY_VARIANT = "summary_statistic_variant"  # MVP max / min of a matched measurement
BASIS_BROAD_ONLY = "icd10_broad_only"
BASIS_NEAR_EQUIVALENT = "measurement_near_equivalent"
BASIS_SEX_STRATUM = "sex_stratum_only"
BASIS_UNMAPPED = "unmapped"

HARMONIZATION_BASIS_LABELS = {
    BASIS_PUBLISHED: "matched by FinnGen's own published phenotype mapping",
    BASIS_ICD10_EXACT: "same single ICD-10 code block (case definitions may still differ)",
    BASIS_ICD10_CODE: "the outcome is an ICD-10 code block",
    BASIS_PHECODE_CODE: "the outcome is a phecode",
    BASIS_PHECODE_SAME: "same phecode (case definitions may still differ by cohort)",
    BASIS_CURATED: "same measurement by a reviewed name match (units and summaries may differ)",
    BASIS_SUMMARY_VARIANT: ("the maximum or minimum of a measurement whose mean is the one "
                            "compared; listed, not compared"),
    BASIS_BROAD_ONLY: "only a broader or narrower ICD-10 counterpart exists; not compared",
    BASIS_NEAR_EQUIVALENT: "only a near-equivalent measurement exists; not compared",
    BASIS_SEX_STRATUM: "a one-sex version of an outcome; not compared",
    BASIS_UNMAPPED: "no counterpart in another biobank",
}

# Unmapped reasons (an enum: a new failure mode cannot merge silently into an old bucket).
UNMAPPED_REASONS = {
    "no_published_definition": "FinnGen publishes no mapping or ICD-10 definition for this endpoint",
    "definition_not_icd10_only": "the endpoint is defined by medication, procedure or other criteria, not ICD-10 codes alone",
    "icd10_block_not_served": "its ICD-10 code block has no Pan-UKBB outcome",
    "definition_unparsable": "its ICD-10 definition uses a code prefix that is not one code block",
    "published_match_not_served": "FinnGen's published match is not a Pan-UKBB outcome",
    "measurement_no_counterpart": "a measurement with no counterpart in another biobank",
    "questionnaire_or_medication": "a questionnaire answer or medication-use outcome, not coded to a disease vocabulary",
    "no_counterpart_code": "an MVP phecode that Pan-UK Biobank does not serve for both sexes",
    "survey_item_no_counterpart": "an MVP survey answer with no counterpart in another biobank",
}

# A FinnGen endpoint's link to an MVP phecode (FinnGen's FinnGen/MVP/UKBB file) that the
# one-key rule does not compare. Kept in the row, never dropped.
SECONDARY_REASONS = {
    "mvp_phecode_not_released": "FinnGen names an MVP phecode that MVP did not release",
    "anchor_code_differs": "the endpoint's Pan-UKBB anchor is a different phecode than its MVP code",
    "anchored_elsewhere": "the endpoint is compared on another key (an ICD-10 block or a Pan-UKBB outcome)",
    "mvp_link_only": "FinnGen links the endpoint to MVP but to no Pan-UKBB outcome",
}

# The precedence orders the representative when several outcomes of ONE dataset share an
# id; a basis in NON_REPRESENTATIVE_BASES keeps its id and is never the representative.
BASIS_PRECEDENCE = {BASIS_PUBLISHED: 0, BASIS_ICD10_EXACT: 1, BASIS_ICD10_CODE: 2,
                    BASIS_PHECODE_CODE: 2, BASIS_PHECODE_SAME: 2, BASIS_CURATED: 3,
                    BASIS_SUMMARY_VARIANT: 9}
NON_REPRESENTATIVE_BASES = frozenset({BASIS_SUMMARY_VARIANT})

CROSSWALK_COLUMNS = (
    "dataset_family", "outcome_id", "outcome_label", "harmonized_outcome_id",
    "harmonized_label", "harmonization_basis", "predicate_id",
    "anchor_outcome_id", "broad_candidates", "unmapped_reason", "mapping_source",
    "source_record_id", "secondary_match", "secondary_reason",
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
    secondary_match: str | None = None
    secondary_reason: str | None = None


def normalize_phecode(code: str) -> str:
    """One spelling of a phecode for a natural join: Pan-UKBB writes `008.5`, the GWAS
    Catalog's MVP records write `8.5`. Leading zeros of the integer part and trailing
    fractional zeros are dropped (`008` -> `8`, `250.10` -> `250.1`, `250.0` -> `250`)."""
    c = (code or "").strip()
    if not re.fullmatch(r"\d+(\.\d+)?", c):
        raise ValueError(f"not a numeric phecode: {code!r}")
    whole, _, frac = c.partition(".")
    whole = whole.lstrip("0") or "0"
    frac = frac.rstrip("0")
    return f"{whole}.{frac}" if frac else whole


def quant_key(base: str, qualifier: str) -> str:
    """The reviewed-table key of an MVP quantitative trait: base + qualifier, case-folded
    (the GWAS Catalog writes `Height` and `height` for one measurement)."""
    b = " ".join(base.split()).casefold()
    q = " ".join(qualifier.split()).casefold()
    return f"{b} [{q}]" if q else b


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
                mapping_source=r["mapping_source"], source_record_id=r["source_record_id"],
                secondary_match=r.get("secondary_match") or None,
                secondary_reason=r.get("secondary_reason") or None)
            out[(row.dataset_family, row.outcome_id)] = row
    return out


def is_comparable(r: CrosswalkRow) -> bool:
    """A row the replication view can compare: it carries a harmonized id and is not a
    summary variant (an MVP measurement's maximum / minimum)."""
    return bool(r.harmonized_outcome_id) and r.harmonization_basis not in NON_REPRESENTATIVE_BASES


@lru_cache(maxsize=16)
def harmonized_outcomes(family: str, with_families: tuple[str, ...],
                        path: str = DEFAULT_CROSSWALK_PATH) -> frozenset[str]:
    """Outcome ids of `family` whose harmonized id an outcome of any of `with_families`
    shares, by the crosswalk alone (whether or not the files can be read)."""
    cw = load_crosswalk(path)
    others = {r.harmonized_outcome_id for r in cw.values()
              if r.dataset_family in with_families and is_comparable(r)}
    return frozenset(r.outcome_id for r in cw.values()
                     if r.dataset_family == family and is_comparable(r)
                     and r.harmonized_outcome_id in others)


def harmonized_pan_ukbb_outcomes(path: str = DEFAULT_CROSSWALK_PATH) -> frozenset[str]:
    """Pan-UKBB outcome ids whose harmonized id a FinnGen endpoint shares (589 on the
    committed crosswalk)."""
    return harmonized_outcomes(FAMILY_PAN_UKBB, (FAMILY_FINNGEN,), path)


def harmonized_finngen_endpoints(path: str = DEFAULT_CROSSWALK_PATH) -> frozenset[str]:
    """FinnGen endpoints with a Pan-UKBB counterpart (622 on the committed crosswalk)."""
    return harmonized_outcomes(FAMILY_FINNGEN, (FAMILY_PAN_UKBB,), path)
