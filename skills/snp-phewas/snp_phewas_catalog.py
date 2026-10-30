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
  MVP       the GWAS Catalog study accession (e.g. `GCST90475929`); one study per (trait,
            population), so each MVP dataset is one population's studies.
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
    short_label: str = ""  # the column / strip header a figure uses ("MVP EUR")
    #: The participant groups the results are computed from (`samples` in datasets.yaml);
    #: two datasets whose groups intersect share participants (`sample_overlap`).
    samples: tuple[str, ...] = ()
    #: The genetic-ancestry group effect sizes are compared within (`ancestry_group`, else
    #: `ancestry`): FinnGen's FIN is EUR.
    ancestry_group: str = ""

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
    # MVP (GWAS Catalog) only: whether the study has a position-indexed file
    # (`snp_phewas_mvp_catalog.ACCESS_STATUSES`); None for sources whose every catalog file
    # is readable.
    access_status: str | None = None
    md5_hex: str | None = None  # the source's published md5 of the file, where it has one

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
    groups = raw["sample_groups"]
    for dsid, d in raw["datasets"].items():
        fam = raw["families"][d["family"]]
        samples = tuple(d.get("samples") or ())
        # Fail loud: a dataset with no declared participants would look independent of all.
        if not samples or any(a not in groups for a in samples):
            raise ValueError(f"datasets.yaml: {dsid} must list its `samples`, each one of "
                             f"`sample_groups`; got {list(samples)}")
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
            short_label=str(d.get("short_label") or d["label"]),
            samples=samples,
            ancestry_group=str(d.get("ancestry_group") or d["ancestry"]),
        )
    return out


def family_label(family: str) -> str:
    """The plain-language name of a dataset family (`Pan-UK Biobank`, `FinnGen`, ...)."""
    return str(_load_yaml()["families"][family]["label"])


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


ANCESTRY_WORDS: dict[str, str] = {
    "EUR": "European", "AFR": "African", "AMR": "admixed American", "EAS": "East Asian",
    "CSA": "Central/South Asian", "MID": "Middle Eastern",
    "META": "multi-population meta-analysis",
}


#: The `ancestry` value of a multi-population meta-analysis (MVP's `mvp_meta`).
META = "META"


def sample_group_label(atom: str) -> str:
    """The plain-language name of a participant group (`sample_groups` in datasets.yaml)."""
    return str(_load_yaml()["sample_groups"][atom])


def _listed_label(spec: DatasetSpec) -> str:
    head, sep, tail = spec.label.partition(", ")
    return f"{head} ({tail})" if sep else spec.label


def overlap_between(da: DatasetSpec, db: DatasetSpec) -> str | None:
    """Why two datasets share participants, as a plain-language sentence; None when their
    declared `samples` do not intersect. The one home of the rule: overlap is
    DECLARED in datasets.yaml, never derived from ids or ancestry strings, so a meta-analysis
    pooling two biobanks is caught against either one."""
    shared = [a for a in da.samples if a in db.samples]
    if not shared:
        return None
    who = " and the ".join(sample_group_label(a) for a in shared)
    if set(db.samples) < set(da.samples):
        return f"{_listed_label(da)} includes the {who}."
    if set(da.samples) < set(db.samples):
        return f"{_listed_label(db)} includes the {who}."
    return f"{_listed_label(da)} and {_listed_label(db)} both include the {who}."


def sample_overlap(a: str, b: str) -> str | None:
    """`overlap_between` for two catalog dataset ids."""
    return overlap_between(get_dataset(a), get_dataset(b))


def cross_ancestry(a: str, b: str) -> bool:
    """True when effect sizes of `a` and `b` are not compared: different `ancestry_group`s, or
    either a multi-population meta-analysis (its pooled effect mixes groups)."""
    ga, gb = get_dataset(a).ancestry_group, get_dataset(b).ancestry_group
    return ga != gb or META in (ga, gb)


def label_in_list(dataset_id: str) -> str:
    """A dataset's label as it stands in a list: "Pan-UK Biobank, European ancestry" becomes
    "Pan-UK Biobank (European ancestry)", so the label's own comma cannot read as a list
    separator ("A, B and C" of three labels read as five sources)."""
    label = get_dataset(dataset_id).label
    head, sep, tail = label.partition(", ")
    return f"{head} ({tail})" if sep else label


def labels_listed(dataset_ids: list[str]) -> str:
    """Several datasets named in a sentence: "A", "A and B", "A, B and C". The one join every
    answer, heading and progress message uses."""
    names = [label_in_list(d) for d in dataset_ids]
    if len(names) <= 2:
        return " and ".join(names)
    return ", ".join(names[:-1]) + " and " + names[-1]


def datasets_for_ancestry(ancestry: str) -> tuple[list[str], list[str]]:
    """(dataset ids, plain sentences) for one genetic-ancestry group: per biobank family,
    the dataset of that group; a family with no such group is left out and a sentence says
    so (the answer must say what it substituted). EUR is the default lookup. META is MVP's
    multi-population meta-analysis, shown beside the other biobanks' default groups."""
    anc = (ancestry or "EUR").upper()
    if anc not in ANCESTRY_WORDS:
        raise UnknownDatasetError(f"no dataset for the ancestry group {ancestry!r}")
    if anc == "EUR":
        return default_dataset_ids(), []
    ds = load_datasets()
    families = list(dict.fromkeys(s.family for s in ds.values()))
    ids: list[str] = []
    said: list[str] = []
    for fam in families:
        group = [d for d, s in ds.items() if s.family == fam and s.ancestry == anc]
        if group:
            ids.append(group[0])
            continue
        if anc == "META":
            fallback = [d for d, s in ds.items() if s.family == fam and s.default]
            ids += fallback[:1]
            # true for FinnGen (one population) and not for Pan-UKBB, whose files carry a
            # meta-analysis this lookup does not serve, so the sentence says what the
            # LOOKUP serves
            said.append(f"This lookup serves no multi-population meta-analysis for "
                        f"{family_label(fam)}, so its default group is shown instead: "
                        f"{ds[fallback[0]].label}.")
            continue
        said.append(f"{family_label(fam)} has no {ANCESTRY_WORDS[anc]} genetic-ancestry "
                    f"results, so it is not part of this lookup.")
    return ids, said


def resolve_dataset_ids(spec: str | list[str] | None) -> list[str]:
    """`None`/`"default"` -> the defaults; `"all"` -> every dataset; else a list/CSV."""
    if spec is None or spec == "default":
        return default_dataset_ids()
    if spec == "all":
        return list(load_datasets())
    ids = [s.strip() for s in spec.split(",")] if isinstance(spec, str) else list(spec)
    ids = list(dict.fromkeys(ids))  # a repeated id would repeat every row
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
    elif ds.family == "mvp":
        for r in rows:
            if r["population"] != ds.ancestry:
                continue
            n_controls = _int_or_none(r["n_controls"])
            binary = n_controls is not None
            out.append(OutcomeSpec(
                dataset_id=dataset_id,
                outcome_id=r["outcome_id"],
                outcome_label=r["outcome_label"],
                outcome_category_native=r["category_native"],
                outcome_kind="binary" if binary else "continuous",
                trait_type=r["trait_kind"],
                n_cases=_int_or_none(r["n_cases"] if binary else r["n"]),
                n_controls=n_controls,
                filename=r["filename"] or None,
                filename_tabix=r["filename_tabix"] or None,
                access_status=r["access_status"],
                md5_hex=r["md5_hex"] or None,
            ))
        # Two studies of one population can share a trait text only when the Catalog lists
        # the trait twice; the accession tells them apart.
        out = _disambiguate(out, {r["outcome_id"]: f"study {r['outcome_id']}" for r in rows})
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
