"""The curated dataset catalog and the per-dataset outcome catalogs.

`datasets.yaml` names every dataset the skill serves (license verdict, release pin,
attribution, access route). Each dataset family ships an in-repo OUTCOME catalog built
from the source's own published manifest by `snp_phewas_build_catalogs.py`; the catalog is the
universe a lookup is accounted against (every catalog outcome yields exactly one row,
whatever happened to it).

Outcome ids are the source's own unique keys (persisted, never re-minted):
  Pan-UKBB  `trait_type-phenocode-pheno_sex[-coding][-modifier]`, composed from the
            phenotype manifest (`snp_phewas_build_catalogs.pan_ukbb_outcome_id`)
            (e.g. `biomarkers-30780-both_sexes-irnt`).
  FinnGen   the endpoint phenocode (e.g. `I9_CHD`).
"""

from __future__ import annotations

import csv
import gzip
import os
from dataclasses import dataclass, field, replace
from functools import lru_cache
from typing import Any

import yaml

_SKILL_DIR = os.path.dirname(os.path.abspath(__file__))
DATASETS_YAML = os.path.join(_SKILL_DIR, "datasets.yaml")

PAN_UKBB_POPS: tuple[str, ...] = ("AFR", "AMR", "CSA", "EAS", "EUR", "MID")

# Plain-language labels for the Pan-UKBB trait_type server enum (user-friendly output rule).
TRAIT_TYPE_LABELS: dict[str, str] = {
    "phecode": "disease (phecode)",
    "icd10": "disease (ICD-10 code)",
    "categorical": "questionnaire answer",
    "continuous": "continuous measurement",
    "prescriptions": "medication use",
    "biomarkers": "blood biomarker",
}


class UnknownDatasetError(ValueError):
    """A dataset id that the curated catalog does not list."""


@dataclass(frozen=True)
class DatasetSpec:
    """One served dataset (a family + an ancestry stratum)."""

    dataset_id: str
    family: str
    ancestry: str
    label: str
    default: bool
    release: str
    build: str
    reference_anchored: bool
    license: str
    license_verdict: str
    license_verified: str
    attribution: str
    family_spec: dict[str, Any] = field(repr=False, default_factory=dict)

    @property
    def outcome_catalog_path(self) -> str:
        return os.path.join(_SKILL_DIR, self.family_spec["outcome_catalog"])


@dataclass(frozen=True)
class OutcomeSpec:
    """One outcome in a dataset's catalog."""

    dataset_id: str
    outcome_id: str
    outcome_label: str
    outcome_category_native: str
    outcome_kind: str  # binary | continuous
    trait_type: str | None
    n_cases: int | None
    n_controls: int | None
    filename: str | None
    filename_tabix: str | None = None

    @property
    def n_total(self) -> int | None:
        if self.n_cases is None:
            return None
        return self.n_cases + (self.n_controls or 0)


@lru_cache(maxsize=1)
def _load_yaml() -> dict[str, Any]:
    with open(DATASETS_YAML, encoding="utf-8") as fh:
        return yaml.safe_load(fh)


def load_datasets() -> dict[str, DatasetSpec]:
    """Every served dataset, keyed by dataset_id, in catalog order."""
    raw = _load_yaml()
    out: dict[str, DatasetSpec] = {}
    for dsid, d in raw["datasets"].items():
        fam = raw["families"][d["family"]]
        out[dsid] = DatasetSpec(
            dataset_id=dsid,
            family=d["family"],
            ancestry=d["ancestry"],
            label=d["label"],
            default=bool(d.get("default", False)),
            release=fam["release"],
            build=fam["build"],
            reference_anchored=bool(fam["reference_anchored"]),
            license=fam["license"],
            license_verdict=fam["license_verdict"],
            license_verified=str(fam["license_verified"]),
            attribution=" ".join(str(fam["attribution"]).split()),
            family_spec=fam,
        )
    return out


def not_served() -> dict[str, str]:
    """Datasets considered and deliberately not served, with the plain-language reason."""
    return dict(_load_yaml().get("not_served", {}))


def get_dataset(dataset_id: str) -> DatasetSpec:
    ds = load_datasets()
    if dataset_id not in ds:
        raise UnknownDatasetError(
            f"{dataset_id!r} is not a curated dataset; choose from {', '.join(ds)}"
        )
    return ds[dataset_id]


def default_dataset_ids() -> list[str]:
    return [k for k, v in load_datasets().items() if v.default]


def resolve_dataset_ids(spec: str | list[str] | None) -> list[str]:
    """`None`/`"default"` -> the defaults; `"all"` -> every dataset; else a list/CSV."""
    if spec is None or spec == "default":
        return default_dataset_ids()
    if spec == "all":
        return list(load_datasets())
    ids = [s.strip() for s in spec.split(",")] if isinstance(spec, str) else list(spec)
    for i in ids:
        get_dataset(i)
    return ids


def _int_or_none(v: str | None) -> int | None:
    if v is None or v in ("", "NA", "nan"):
        return None
    return int(float(v))


def _read_tsv_gz(path: str) -> list[dict[str, str]]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


@lru_cache(maxsize=16)
def _load_outcome_catalog_cached(dataset_id: str) -> tuple[OutcomeSpec, ...]:
    ds = get_dataset(dataset_id)
    rows = _read_tsv_gz(ds.outcome_catalog_path)
    out: list[OutcomeSpec] = []
    if ds.family == "pan_ukbb":
        pop = ds.ancestry
        for r in rows:
            if pop not in r["pops"].split(","):
                continue
            n_cases = _int_or_none(r[f"n_cases_{pop}"])
            n_controls = _int_or_none(r[f"n_controls_{pop}"])
            out.append(OutcomeSpec(
                dataset_id=dataset_id,
                outcome_id=r["outcome_id"],
                outcome_label=_pan_label(r),
                outcome_category_native=r["category"],
                outcome_kind="binary" if n_controls is not None else "continuous",
                trait_type=r["trait_type"],
                n_cases=n_cases,
                n_controls=n_controls,
                filename=r["filename"],
                filename_tabix=r["filename_tabix"],
            ))
        out = _disambiguate(out, {r["outcome_id"]: (r.get("modifier") or "").strip() for r in rows})
        out = _disambiguate(out, {r["outcome_id"]: (f"UK Biobank field {r['phenocode']}"
                                                    if r["phenocode"].isdigit() else r["trait_type"])
                                  for r in rows})
        # Last: two codings of one field that the source itself describes identically.
        out = _disambiguate(out, {r["outcome_id"]: f"code {r['coding']}" for r in rows
                                  if (r.get("coding") or "").strip()})
    elif ds.family == "finngen":
        for r in rows:
            n_controls = _int_or_none(r["num_controls"])
            out.append(OutcomeSpec(
                dataset_id=dataset_id,
                outcome_id=r["outcome_id"],
                outcome_label=r["phenotype"] or f"FinnGen endpoint {r['outcome_id']}",
                outcome_category_native=r["category"],
                outcome_kind="binary" if n_controls else "continuous",
                trait_type=None,
                n_cases=_int_or_none(r["num_cases"]),
                n_controls=n_controls,
                filename=r["filename"],
                filename_tabix=r["filename"] + ".tbi",
            ))
    else:  # pragma: no cover - the yaml lists only these families
        raise UnknownDatasetError(f"no outcome-catalog reader for family {ds.family!r}")
    ids = [o.outcome_id for o in out]
    if len(set(ids)) != len(ids):
        raise ValueError(f"outcome catalog for {dataset_id} has duplicate outcome ids")
    return tuple(out)


_SEX_LABELS = {"females": "women", "males": "men"}
_MODIFIER_WORDS = {"irnt": "rank-normalized", "raw": "raw units", "medadj": "medication-adjusted"}


def _pan_label(r: dict) -> str:
    """A label that tells two Pan-UKBB outcomes apart. The manifest's `description` alone
    is shared: every self-reported illness is "Non-cancer illness code, self-reported" (the
    illness is in `coding_description`), and the rank-normalized and raw-unit versions of
    one measurement carry the same text. The modifier is added by `_disambiguate` only where
    a label would otherwise repeat."""
    desc = (r["description"] or "").strip()
    label = desc if desc and desc.upper() != "NA" else r["phenocode"]
    cd = (r.get("coding_description") or "").strip()
    if cd and cd.upper() != "NA":
        label = f"{label}: {cd}"
    if r.get("pheno_sex") in _SEX_LABELS:
        label = f"{label} ({_SEX_LABELS[r['pheno_sex']]})"
    return label


def _modifier_words(mod: str) -> str:
    if " " in mod:  # already words (the field-number fallback)
        return mod
    return ", ".join(_MODIFIER_WORDS.get(t, t) for t in mod.split("_") if t)


def _disambiguate(specs: list, modifiers: dict[str, str]) -> list:
    """Append the plain-language modifier (the manifest's `modifier` column) to a Pan-UKBB
    label that repeats within a dataset."""
    from collections import Counter

    counts = Counter(o.outcome_label for o in specs)
    out = []
    for o in specs:
        mod = modifiers.get(o.outcome_id, "")
        if counts[o.outcome_label] > 1 and mod:
            o = replace(o, outcome_label=f"{o.outcome_label} ({_modifier_words(mod)})")
        out.append(o)
    return out


def load_outcome_catalog(dataset_id: str) -> list[OutcomeSpec]:
    """The outcome universe of one dataset (one entry per outcome, unique outcome_id)."""
    return list(_load_outcome_catalog_cached(dataset_id))
