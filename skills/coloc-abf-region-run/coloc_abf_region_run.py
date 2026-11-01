#!/usr/bin/env python3
"""coloc-abf-region-run: colocalise two traits over one genomic region.

Reads a small config naming two regional summary-statistics files, pairs their
variants, runs coloc.abf, sweeps the p12 prior, and writes a report, a machine
readable result, the per-p12 table, a sensitivity figure and a reproducibility
bundle.

    python skills/coloc-abf-region-run/coloc_abf_region_run.py \
        --input config.yaml --output results/

    python skills/coloc-abf-region-run/coloc_abf_region_run.py --demo --output /tmp/coloc

The arithmetic lives in `coloc_abf_math.py`; this module is the I/O, the variant
pairing and the report.

Config schema (YAML or JSON; paths are relative to the config file):

    trait1:
      sumstats: sumstats/eqtl_region.tsv
      trait_type: quant        # "quant" or "cc"
      sdY: 1.0                 # quant: give sdY, or give n and an eaf/maf column
      n: 50000                 # quant: used to estimate sdY when sdY is absent
      label: SORT1 eQTL (liver)
    trait2:
      sumstats: sumstats/gwas_region.tsv
      trait_type: quant
      sdY: 1.0
      label: LDL cholesterol
    priors:                    # optional
      p1: 1.0e-4
      p2: 1.0e-4
      p12: 1.0e-5
    rule_h4_min: 0.5           # optional: the rule the p12 sweep is judged by
"""

from __future__ import annotations

import argparse
import csv
import json
import shlex
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Sequence

import numpy as np
import yaml

_SKILL_DIR = Path(__file__).resolve().parent
_PROJECT_ROOT = _SKILL_DIR.parents[1]
if str(_SKILL_DIR) not in sys.path:
    sys.path.insert(0, str(_SKILL_DIR))
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from clawbio.common.report import skill_version  # noqa: E402
from clawbio.common.reproducibility import (  # noqa: E402
    write_checksums,
    write_commands_sh,
    write_environment_yml,
)
from coloc_abf_math import (  # noqa: E402
    DEFAULT_P1,
    DEFAULT_P2,
    DEFAULT_P12,
    DEFAULT_RULE_H4_MIN,
    ColocDataset,
    ColocError,
    aligned_log_abf,
    coloc_sensitivity,
    run_coloc_abf,
)

SKILL_NAME = "coloc-abf-region-run"
SKILL_VERSION = skill_version(__file__)  # SKILL.md is the only place it is written

DISCLAIMER = (
    "ClawBio is a research and educational tool. It is not a medical device "
    "and does not provide clinical diagnoses. Consult a healthcare "
    "professional before making any medical decisions."
)

# Why a variant read from a file did not reach the colocalisation. Every excluded
# variant is written to tables/excluded_variants.tsv with one of these, and the
# kept and excluded sets are asserted to partition what was read.
EXCLUDED_MISSING_STATS = "missing_beta_or_se"
EXCLUDED_INVALID_SE = "non_positive_or_non_finite_se"
EXCLUDED_MISSING_ALLELES = "missing_chromosome_position_or_alleles"
EXCLUDED_DUPLICATE_ROW = "duplicate_row_same_statistics"
EXCLUDED_NOT_SHARED = "variant_not_present_in_the_other_trait"

_REQUIRED_COLUMNS = ("variant_id", "chromosome", "position", "ref", "alt", "beta", "se")

# Header aliases, matched case-insensitively, so a file straight out of a region
# fetch or a harmoniser loads without being rewritten.
_ALIASES = {
    "snp": "variant_id",
    "rsid": "variant_id",
    "variant": "variant_id",
    "chr": "chromosome",
    "chrom": "chromosome",
    "pos": "position",
    "bp": "position",
    "base_pair_location": "position",
    "effect_allele": "alt",
    "other_allele": "ref",
    "ea": "alt",
    "nea": "ref",
    "a1": "alt",
    "a2": "ref",
    "p": "p_value",
    "pval": "p_value",
    "pvalue": "p_value",
    "effect_allele_frequency": "eaf",
    "freq": "eaf",
    "af": "eaf",
    "sample_size": "n",
    "n_samples": "n",
}


class ColocConfigError(ValueError):
    """Raised when the config is missing a field or names a file that is absent."""


class ColocInputError(ValueError):
    """Raised when a summary-statistics file cannot be read."""


@dataclass(frozen=True)
class SumstatRow:
    """One variant's summary statistics, as read from a file."""

    variant_id: str
    chromosome: str
    position: int
    ref: str
    alt: str
    beta: float | None
    se: float | None
    p_value: float | None
    eaf: float | None
    n: float | None

    @property
    def identity(self) -> str:
        """A key naming the variant rather than the row.

        Chromosome and position plus the allele pair, sorted, so one variant
        written in either allele order gets one key and the two alleles of a
        multi-allelic site get two.
        """
        a, b = sorted((self.alt.upper(), self.ref.upper()))
        chrom = self.chromosome.lower().removeprefix("chr").upper()
        return f"{chrom}_{self.position}_{a}_{b}"


# ---------------------------------------------------------------------------
# Reading summary statistics
# ---------------------------------------------------------------------------

def _canonical(column: str) -> str:
    c = column.strip().lower().lstrip("#")
    return _ALIASES.get(c, c)


def _to_float(value: str | None) -> float | None:
    if value is None:
        return None
    v = value.strip()
    if v == "" or v.upper() in {"NA", "NAN", "NULL", "."}:
        return None
    try:
        return float(v)
    except ValueError:
        return None


def load_regional_sumstats(path: str | Path) -> list[SumstatRow]:
    """Read a regional summary-statistics table into rows.

    Required columns, after alias mapping: variant_id, chromosome, position, ref,
    alt, beta, se. Optional: p_value, eaf (or maf), n. Tab-separated by default;
    a comma-separated file with no tabs is also accepted.
    """
    path = Path(path)
    if not path.exists():
        raise ColocInputError(f"summary statistics file not found: {path}")

    with path.open(newline="", encoding="utf-8") as fh:
        first = fh.readline()
        fh.seek(0)
        delimiter = "," if ("," in first and "\t" not in first) else "\t"
        reader = csv.DictReader(fh, delimiter=delimiter)
        if reader.fieldnames is None:
            raise ColocInputError(f"empty summary statistics file: {path}")
        field_map = {name: _canonical(name) for name in reader.fieldnames}
        present = set(field_map.values())
        missing = [c for c in _REQUIRED_COLUMNS if c not in present]
        if missing:
            raise ColocInputError(
                f"{path}: missing required column(s) {missing}. Columns read: "
                f"{sorted(present)}."
            )

        rows: list[SumstatRow] = []
        for line_no, raw in enumerate(reader, start=2):
            row = {field_map[k]: v for k, v in raw.items() if k in field_map}
            eaf = _to_float(row.get("eaf"))
            if eaf is None:
                eaf = _to_float(row.get("maf"))
            try:
                position = int(float(row["position"]))
            except (TypeError, ValueError) as exc:
                raise ColocInputError(f"{path} line {line_no}: bad position: {exc}") from exc
            rows.append(SumstatRow(
                variant_id=(row.get("variant_id") or "").strip(),
                chromosome=str(row.get("chromosome") or "").strip(),
                position=position,
                ref=(row.get("ref") or "").strip(),
                alt=(row.get("alt") or "").strip(),
                beta=_to_float(row.get("beta")),
                se=_to_float(row.get("se")),
                p_value=_to_float(row.get("p_value")),
                eaf=eaf,
                n=_to_float(row.get("n")),
            ))

    if not rows:
        raise ColocInputError(f"{path}: no data rows")
    return rows


# ---------------------------------------------------------------------------
# Pairing the two traits' variants
# ---------------------------------------------------------------------------

@dataclass
class PairedRegion:
    """The variants both traits measured, plus every variant that was excluded."""

    keys: list[str]
    display_ids: list[str]
    rows1: list[SumstatRow]
    rows2: list[SumstatRow]
    excluded: list[dict[str, Any]]
    n_read1: int
    n_read2: int


def _usable(row: SumstatRow) -> str | None:
    """The exclusion reason for a row, or None when it can be used."""
    if not (row.chromosome and row.ref and row.alt) or row.position <= 0:
        return EXCLUDED_MISSING_ALLELES
    if row.beta is None or row.se is None:
        return EXCLUDED_MISSING_STATS
    if not np.isfinite(row.beta) or not np.isfinite(row.se) or row.se <= 0:
        return EXCLUDED_INVALID_SE
    return None


def _excluded_record(trait: str, row: SumstatRow, reason: str) -> dict[str, Any]:
    return {
        "trait": trait,
        "variant_id": row.variant_id,
        "chromosome": row.chromosome,
        "position": row.position,
        "ref": row.ref,
        "alt": row.alt,
        "reason": reason,
    }


def _index_one_trait(
    trait: str, rows: Sequence[SumstatRow], excluded: list[dict[str, Any]]
) -> dict[str, SumstatRow]:
    """Index one trait's rows by variant, recording every row that is dropped.

    Two rows for one variant carrying the same statistics are a duplicated row,
    so the second is dropped and recorded. Two rows for one variant whose
    statistics differ cannot be reconciled here, so the run stops: picking one
    would silently choose an answer.
    """
    by_variant: dict[str, SumstatRow] = {}
    for row in rows:
        reason = _usable(row)
        if reason is not None:
            excluded.append(_excluded_record(trait, row, reason))
            continue
        key = row.identity
        seen = by_variant.get(key)
        if seen is None:
            by_variant[key] = row
            continue
        if (seen.beta, seen.se) == (row.beta, row.se):
            excluded.append(_excluded_record(trait, row, EXCLUDED_DUPLICATE_ROW))
            continue
        raise ColocError(
            f"{trait}: variant {key} appears on two rows whose beta or standard error "
            f"differ (identifiers {seen.variant_id!r} and {row.variant_id!r}). "
            "coloc.abf needs one row per variant and cannot choose between them, so "
            "keep one row, or split the rows by molecular trait if they belong to "
            "different transcripts or splicing events."
        )
    return by_variant


def pair_region(
    rows1: Sequence[SumstatRow],
    rows2: Sequence[SumstatRow],
    *,
    label1: str = "trait1",
    label2: str = "trait2",
) -> PairedRegion:
    """Pair two traits' variants on chromosome, position and allele pair.

    Returns the shared variants in the order the first trait lists them, plus one
    record per excluded variant. The kept and excluded sets are asserted to
    partition what was read from each file.
    """
    excluded: list[dict[str, Any]] = []
    by_variant1 = _index_one_trait(label1, rows1, excluded)
    by_variant2 = _index_one_trait(label2, rows2, excluded)

    shared = [k for k in by_variant1 if k in by_variant2]
    for key, row in by_variant1.items():
        if key not in by_variant2:
            excluded.append(_excluded_record(label1, row, EXCLUDED_NOT_SHARED))
    for key, row in by_variant2.items():
        if key not in by_variant1:
            excluded.append(_excluded_record(label2, row, EXCLUDED_NOT_SHARED))

    # A variant keeps the identifier the first trait gave it, so the reported lead
    # comes back in the user's own terms; an identifier naming more than one
    # variant falls back to the key, which names exactly one.
    id_counts: dict[str, int] = {}
    for key in by_variant1:
        vid = by_variant1[key].variant_id
        id_counts[vid] = id_counts.get(vid, 0) + 1
    display_ids = []
    for key in shared:
        vid = by_variant1[key].variant_id
        display_ids.append(vid if vid and id_counts.get(vid) == 1 else key)

    for trait, n_read in ((label1, len(rows1)), (label2, len(rows2))):
        n_excluded = sum(1 for e in excluded if e["trait"] == trait)
        if len(shared) + n_excluded != n_read:
            raise ColocError(
                f"{trait}: {len(shared)} shared plus {n_excluded} excluded does not "
                f"account for the {n_read} rows read. This is a defect in the pairing, "
                "not in your data; please report it."
            )

    return PairedRegion(
        keys=shared,
        display_ids=display_ids,
        rows1=[by_variant1[k] for k in shared],
        rows2=[by_variant2[k] for k in shared],
        excluded=excluded,
        n_read1=len(rows1),
        n_read2=len(rows2),
    )


def dataset_from_rows(
    keys: Sequence[str],
    rows: Sequence[SumstatRow],
    trait_type: str,
    *,
    label: str,
    sdY: float | None,
    n: float | None,
) -> ColocDataset:
    """Build a dataset for the computation, folding EAF to MAF when sdY is absent."""
    maf = None
    if trait_type == "quant" and sdY is None:
        eafs = [r.eaf for r in rows]
        if any(e is None for e in eafs):
            raise ColocError(
                f"{label}: a quantitative trait without sdY needs an allele frequency "
                "on every shared variant so that sdY can be estimated; supply sdY in "
                "the config instead, or add an eaf (or maf) column."
            )
        maf = np.minimum(np.asarray(eafs, dtype=float), 1.0 - np.asarray(eafs, dtype=float))
        if n is None:
            per_row = [r.n for r in rows]
            if any(v is None for v in per_row):
                raise ColocError(
                    f"{label}: a quantitative trait without sdY needs a sample size, "
                    "either as n in the config or as an n column in the file."
                )
            n = np.asarray(per_row, dtype=float)
    return ColocDataset(
        snp=list(keys),
        beta=np.asarray([r.beta for r in rows], dtype=float),
        se=np.asarray([r.se for r in rows], dtype=float),
        trait_type=trait_type,  # type: ignore[arg-type]
        label=label,
        sdY=sdY,
        maf=maf,
        n=n,
    )


# ---------------------------------------------------------------------------
# Caveats
# ---------------------------------------------------------------------------

def build_caveats(result, sensitivity, paired: PairedRegion) -> list[str]:
    """What the numbers do not say on their own."""
    caveats: list[str] = list(result.notes)
    caveats.extend(sensitivity.notes)

    if result.h3 > result.h4 and result.h3 > 0.5:
        caveats.append(
            "PP.H3 is the largest posterior: both traits are associated here but at "
            "different causal variants, which is not a colocalisation."
        )
    if result.h3 > 0.1 and result.h4 > 0.1:
        caveats.append(
            "PP.H3 and PP.H4 both hold non-trivial mass, which is what a region with "
            "more than one independent signal looks like. coloc.abf assumes at most "
            "one causal variant per trait, so fine-map the region and colocalise the "
            "credible sets before concluding."
        )
    if result.n_colocalising_variants <= 2:
        caveats.append(
            f"95% of the shared-variant evidence sits on {result.n_colocalising_variants} "
            "variant(s), so the posteriors hinge on that row's effect and standard error. "
            "One mis-imputed or mis-aligned row at the peak would move the answer."
        )
    if result.n_snps < 50:
        caveats.append(
            f"Only {result.n_snps} variants were shared by both traits. A window this "
            "sparse can miss the causal variant altogether, which pushes mass onto H3."
        )
    if any(result.sdy_estimated.values()):
        caveats.append(
            "sdY was estimated from allele frequency, standard error and sample size "
            "for at least one quantitative trait. Supplying the measured trait standard "
            "deviation gives a tighter prior."
        )
    n_excluded = len(paired.excluded)
    if n_excluded:
        caveats.append(
            f"{n_excluded} variant record(s) did not reach the computation; each one is "
            "listed with its reason in tables/excluded_variants.tsv."
        )
    return caveats


# ---------------------------------------------------------------------------
# Running one pair
# ---------------------------------------------------------------------------

def run_pair(
    rows1: Sequence[SumstatRow],
    rows2: Sequence[SumstatRow],
    *,
    trait1_type: str,
    trait2_type: str,
    label1: str,
    label2: str,
    sdY1: float | None = None,
    sdY2: float | None = None,
    n1: float | None = None,
    n2: float | None = None,
    p1: float = DEFAULT_P1,
    p2: float = DEFAULT_P2,
    p12: float = DEFAULT_P12,
    rule_h4_min: float = DEFAULT_RULE_H4_MIN,
) -> dict[str, Any]:
    """Pair, colocalise, sweep p12, and return everything the writers need."""
    paired = pair_region(rows1, rows2, label1=label1, label2=label2)
    d1 = dataset_from_rows(paired.keys, paired.rows1, trait1_type,
                           label=label1, sdY=sdY1, n=n1)
    d2 = dataset_from_rows(paired.keys, paired.rows2, trait2_type,
                           label=label2, sdY=sdY2, n=n2)

    aligned = aligned_log_abf(d1, d2)
    result = run_coloc_abf(d1, d2, p1=p1, p2=p2, p12=p12, aligned=aligned)
    sensitivity = coloc_sensitivity(
        d1, d2, p1=p1, p2=p2, p12_default=p12,
        rule_h4_min=rule_h4_min, aligned=aligned,
    )

    lead_display = paired.display_ids[paired.keys.index(result.best_snp)]
    return {
        "paired": paired,
        "result": result,
        "sensitivity": sensitivity,
        "lead_display_id": lead_display,
        "labels": (label1, label2),
        "trait_types": (trait1_type, trait2_type),
        "caveats": build_caveats(result, sensitivity, paired),
    }


# ---------------------------------------------------------------------------
# Writing the output
# ---------------------------------------------------------------------------

def hypothesis_labels(label1: str = "the first trait",
                      label2: str = "the second trait") -> dict[str, str]:
    """What each of the five hypotheses means, in the traits' own names."""
    return {
        "PP.H0": "H0: neither trait is associated in this region",
        "PP.H1": f"H1: only {label1} is associated",
        "PP.H2": f"H2: only {label2} is associated",
        "PP.H3": "H3: both are associated, at different causal variants",
        "PP.H4": "H4: both are associated, at one shared causal variant",
    }


def write_outputs(run: dict[str, Any], output_dir: Path, invocation: str) -> dict[str, Path]:
    """Write the report, the result, the tables, the figure and the bundle."""
    output_dir = Path(output_dir)
    (output_dir / "tables").mkdir(parents=True, exist_ok=True)
    (output_dir / "figures").mkdir(parents=True, exist_ok=True)

    result = run["result"]
    sens = run["sensitivity"]
    paired: PairedRegion = run["paired"]
    label1, label2 = run["labels"]

    labels = hypothesis_labels(label1, label2)
    posteriors_tsv = output_dir / "tables" / "posteriors.tsv"
    with posteriors_tsv.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["hypothesis", "posterior", "meaning"])
        for key, value in result.pp.items():
            w.writerow([key, f"{value:.6g}", labels[key]])

    sensitivity_tsv = output_dir / "tables" / "sensitivity_p12.tsv"
    with sensitivity_tsv.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["p12", "PP.H3", "PP.H4", "rule_passes"])
        for pt in sens.points:
            w.writerow([f"{pt.p12:.6g}", f"{pt.h3:.6g}", f"{pt.h4:.6g}",
                        "yes" if pt.rule_passes else "no"])

    excluded_tsv = output_dir / "tables" / "excluded_variants.tsv"
    with excluded_tsv.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["trait", "variant_id", "chromosome", "position", "ref", "alt", "reason"])
        for rec in paired.excluded:
            w.writerow([rec["trait"], rec["variant_id"], rec["chromosome"],
                        rec["position"], rec["ref"], rec["alt"], rec["reason"]])

    figure = _write_sensitivity_figure(sens, output_dir / "figures" / "sensitivity_p12.png",
                                       label1, label2)

    reason_counts: dict[str, int] = {}
    for rec in paired.excluded:
        reason_counts[rec["reason"]] = reason_counts.get(rec["reason"], 0) + 1

    payload = {
        "skill": SKILL_NAME,
        "skill_version": SKILL_VERSION,
        "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "method": "coloc.abf: single causal variant per trait, Wakefield approximate "
                  "Bayes factors combined as in Giambartolomei 2014",
        "traits": {
            "trait1": {"label": label1, "trait_type": run["trait_types"][0],
                       "variants_read": paired.n_read1,
                       "sdY": result.sdy_used["trait1"],
                       "sdY_estimated": result.sdy_estimated["trait1"]},
            "trait2": {"label": label2, "trait_type": run["trait_types"][1],
                       "variants_read": paired.n_read2,
                       "sdY": result.sdy_used["trait2"],
                       "sdY_estimated": result.sdy_estimated["trait2"]},
        },
        "priors": result.priors,
        "prior_sd": result.prior_sd,
        "posteriors": result.pp,
        "n_shared_variants": result.n_snps,
        "n_colocalising_variants": result.n_colocalising_variants,
        "lead_variant": run["lead_display_id"],
        "lead_variant_key": result.best_snp,
        "lead_variant_h4_share": result.best_snp_pp_h4,
        "sensitivity": {
            "rule": sens.rule,
            "verdict": sens.verdict,
            "shared_signal_at_default_p12": sens.passes_at_default,
            "prior_driven": sens.prior_driven,
            "robust": sens.robust,
            "passing_p12_range": list(sens.passing_p12_range) if sens.passing_p12_range else None,
            "points": [{"p12": pt.p12, "PP.H3": pt.h3, "PP.H4": pt.h4,
                        "rule_passes": pt.rule_passes} for pt in sens.points],
        },
        "excluded_variants": {"total": len(paired.excluded), "by_reason": reason_counts},
        "caveats": run["caveats"],
        "disclaimer": DISCLAIMER,
    }
    result_json = output_dir / "result.json"
    result_json.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n",
                           encoding="utf-8")

    report_md = output_dir / "report.md"
    report_md.write_text(_render_report(run, payload), encoding="utf-8")

    write_commands_sh(output_dir, invocation)
    write_environment_yml(
        output_dir,
        env_name="clawbio-coloc-abf-region-run",
        pip_deps=[],
        conda_deps=["numpy>=2.0", "matplotlib>=3.7", "pyyaml>=6.0"],
    )
    write_checksums(
        [report_md, result_json, posteriors_tsv, sensitivity_tsv, excluded_tsv, figure],
        output_dir,
        anchor=output_dir,
    )
    return {
        "report": report_md,
        "result": result_json,
        "posteriors": posteriors_tsv,
        "sensitivity": sensitivity_tsv,
        "excluded": excluded_tsv,
        "figure": figure,
    }


def _write_sensitivity_figure(sens, path: Path, label1: str, label2: str) -> Path:
    """Plot PP.H4 and PP.H3 against p12, with the default p12 and the rule marked."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    p12 = [pt.p12 for pt in sens.points]
    h4 = [pt.h4 for pt in sens.points]
    h3 = [pt.h3 for pt in sens.points]

    fig, ax = plt.subplots(figsize=(6.0, 4.0))
    ax.plot(p12, h4, marker="o", label="PP.H4 (shared variant)")
    ax.plot(p12, h3, marker="s", label="PP.H3 (distinct variants)")
    ax.axhline(sens.rule_h4_min, linestyle=":", linewidth=1,
               label=f"rule: PP.H4 = {sens.rule_h4_min:g}")
    ax.axvline(sens.p12_default, linestyle="--", linewidth=1,
               label=f"p12 = {sens.p12_default:g}")
    ax.set_xscale("log")
    ax.set_xlabel("p12 (prior that one variant is causal for both traits)")
    ax.set_ylabel("posterior probability")
    ax.set_ylim(-0.05, 1.05)
    ax.set_title(f"Prior sensitivity: {label1} vs {label2}", fontsize=10)
    ax.legend(fontsize=7, loc="best")
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def _render_report(run: dict[str, Any], payload: dict[str, Any]) -> str:
    result = run["result"]
    sens = run["sensitivity"]
    paired: PairedRegion = run["paired"]
    label1, label2 = run["labels"]

    lines: list[str] = []
    lines.append(f"# Colocalisation: {label1} vs {label2}")
    lines.append("")
    lines.append(f"**Generated**: {payload['generated_utc']}")
    lines.append("**Method**: coloc.abf, one causal variant per trait "
                 "(Giambartolomei 2014; Wakefield 2009)")
    lines.append(f"**Shared variants**: {result.n_snps} "
                 f"({paired.n_read1} read for {label1}, {paired.n_read2} for {label2})")
    lines.append("")
    lines.append("## Posteriors")
    lines.append("")
    lines.append("| Hypothesis | Posterior |")
    lines.append("|---|---|")
    labels = hypothesis_labels(label1, label2)
    for key, value in result.pp.items():
        lines.append(f"| {labels[key]} | {value:.4f} |")
    lines.append("")
    lines.append(f"- Candidate shared variant: `{run['lead_display_id']}` "
                 f"(holds {result.best_snp_pp_h4:.3f} of the shared-variant evidence)")
    lines.append(f"- 95% of that evidence sits on {result.n_colocalising_variants} variant(s)")
    lines.append("")
    lines.append("## Prior sensitivity")
    lines.append("")
    lines.append(f"A single PP.H4 depends on `p12`, the prior that one variant is causal "
                 f"for both traits, so the posteriors are recomputed across a range of it "
                 f"(Wallace 2020).")
    lines.append("")
    lines.append(f"- Rule swept: {sens.rule}, over p12 from "
                 f"{sens.points[0].p12:.0e} to {sens.points[-1].p12:.0e}")
    lines.append(f"- Verdict: **{sens.verdict}**")
    if sens.passing_p12_range:
        lo, hi = sens.passing_p12_range
        lines.append(f"- The rule holds for p12 between {lo:.0e} and {hi:.0e}")
    lines.append(f"- See `figures/sensitivity_p12.png` and `tables/sensitivity_p12.tsv`")
    lines.append("")
    lines.append("## Priors used")
    lines.append("")
    p = result.priors
    lines.append(f"- p1 = {p['p1']:g}, p2 = {p['p2']:g}, p12 = {p['p12']:g}")
    lines.append(f"- Prior effect SD: {result.prior_sd['trait1']:g} for {label1}, "
                 f"{result.prior_sd['trait2']:g} for {label2}")
    for side, label in (("trait1", label1), ("trait2", label2)):
        sdy = result.sdy_used[side]
        if sdy is not None:
            how = "estimated from allele frequency, standard error and sample size" \
                if result.sdy_estimated[side] else "supplied"
            lines.append(f"- sdY for {label}: {sdy:.4g} ({how})")
    lines.append("")
    if run["caveats"]:
        lines.append("## What this does not say")
        lines.append("")
        for caveat in run["caveats"]:
            lines.append(f"- {caveat}")
        lines.append("")
    lines.append("---")
    lines.append("")
    lines.append(f"*{DISCLAIMER}*")
    lines.append("")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Config and demos
# ---------------------------------------------------------------------------

def _load_config(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise ColocConfigError(f"config not found: {path}")
    with path.open(encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)  # reads YAML, and JSON as a subset of it
    if not isinstance(cfg, dict):
        raise ColocConfigError(f"{path}: the config must be a mapping")
    for side in ("trait1", "trait2"):
        block = cfg.get(side)
        if not isinstance(block, dict):
            raise ColocConfigError(f"{path}: missing the '{side}' block")
        if "sumstats" not in block:
            raise ColocConfigError(f"{path}: {side}.sumstats is required")
        if block.get("trait_type") not in ("quant", "cc"):
            raise ColocConfigError(
                f"{path}: {side}.trait_type must be 'quant' (a quantitative trait) "
                "or 'cc' (a case-control trait)"
            )
    return cfg


def run_from_config(config_path: Path, output_dir: Path, invocation: str) -> dict[str, Any]:
    """Read a config, run the pair, write everything, and return the run."""
    cfg = _load_config(config_path)
    base = config_path.parent
    t1, t2 = cfg["trait1"], cfg["trait2"]
    priors = cfg.get("priors") or {}

    run = run_pair(
        load_regional_sumstats(base / t1["sumstats"]),
        load_regional_sumstats(base / t2["sumstats"]),
        trait1_type=t1["trait_type"],
        trait2_type=t2["trait_type"],
        label1=t1.get("label", "trait 1"),
        label2=t2.get("label", "trait 2"),
        sdY1=t1.get("sdY"),
        sdY2=t2.get("sdY"),
        n1=t1.get("n"),
        n2=t2.get("n"),
        p1=float(priors.get("p1", DEFAULT_P1)),
        p2=float(priors.get("p2", DEFAULT_P2)),
        p12=float(priors.get("p12", DEFAULT_P12)),
        rule_h4_min=float(cfg.get("rule_h4_min", DEFAULT_RULE_H4_MIN)),
    )
    run["paths"] = write_outputs(run, output_dir, invocation)
    return run


def available_demos() -> dict[str, Path]:
    """The bundled demo configs, by name."""
    return {p.stem: p for p in sorted((_SKILL_DIR / "examples").glob("*.yaml"))}


DEFAULT_DEMO = "sort1_eqtl_x_ldl_gwas"


def headline(run: dict[str, Any]) -> str:
    """One line for a log or a chat surface."""
    result, sens = run["result"], run["sensitivity"]
    return (
        f"PP.H4={result.h4:.3f} PP.H3={result.h3:.3f} "
        f"shared_variants={result.n_snps} lead={run['lead_display_id']} "
        f"sensitivity={sens.verdict}"
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog=SKILL_NAME,
        description=(
            "Colocalise two traits over one genomic region with coloc.abf. Computes "
            "the five posteriors from each trait's regional beta and standard error, "
            "and reports whether the conclusion survives the p12 prior."
        ),
    )
    parser.add_argument("--input", type=Path,
                        help="Config file (YAML or JSON) naming the two traits.")
    parser.add_argument("--output", type=Path,
                        help="Output directory; created if missing.")
    parser.add_argument("--demo", nargs="?", const=DEFAULT_DEMO, default=None,
                        metavar="NAME",
                        help="Run a bundled demo on synthetic data. Bare --demo runs "
                             f"{DEFAULT_DEMO}; see --list-demos for the others.")
    parser.add_argument("--list-demos", action="store_true",
                        help="List the bundled demos and exit.")
    args = parser.parse_args(argv)

    demos = available_demos()
    if args.list_demos:
        print("Bundled demos (synthetic data):")
        for name, path in demos.items():
            print(f"  {name}  ({path.relative_to(_SKILL_DIR)})")
        return 0

    if args.demo is not None:
        if args.demo not in demos:
            print(f"[{SKILL_NAME}] unknown demo {args.demo!r}; choose one of "
                  f"{sorted(demos)}", file=sys.stderr)
            return 2
        config_path = demos[args.demo]
    elif args.input is not None:
        config_path = args.input
    else:
        parser.error("one of --input or --demo is required")

    if args.output is None:
        parser.error("--output is required")

    invocation = shlex.join([
        "python", f"skills/{_SKILL_DIR.name}/{Path(__file__).name}",
        *(["--demo", args.demo] if args.demo is not None else ["--input", str(config_path)]),
        "--output", str(args.output),
    ])

    try:
        run = run_from_config(config_path, args.output, invocation)
    except (ColocConfigError, ColocInputError, ColocError) as exc:
        print(f"[{SKILL_NAME}] cannot run: {exc}", file=sys.stderr)
        return 2

    print(f"[{SKILL_NAME}] {headline(run)}")
    for name, path in run["paths"].items():
        print(f"[{SKILL_NAME}] {name}: {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
