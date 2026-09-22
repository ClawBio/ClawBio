# nfcore-sarek-wrapper / outputs_parser.py
"""Discover and catalog outputs produced by an nf-core/sarek 3.8.1 run.

This module ONLY reads from <outdir>/.  It performs no IO writes, no
downstream handoff, and no provenance bookkeeping.  Callers should treat
the resulting :class:`OutputsReport` as immutable.

Implements the file layout documented in AGENTS.md / instructions, derived
from the official nf-core/sarek 3.8.1 output spec.
"""

from __future__ import annotations

import csv
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

_SKILL_DIR = Path(__file__).resolve().parent
if str(_SKILL_DIR) in sys.path:
    sys.path.remove(str(_SKILL_DIR))
sys.path.insert(0, str(_SKILL_DIR))


sys.modules.pop("_isolated_imports", None)
from _isolated_imports import purge_foreign_bare_modules

purge_foreign_bare_modules("errors", "schemas")

from errors import ErrorCode, SkillError  # noqa: E402

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_PAIR_RE = re.compile(r"^(?P<tumor>[^/]+)_vs_(?P<normal>[^/]+)$")

# Subset of SUPPORTED_TOOLS that produce per-sample/per-pair VCFs we expect to
# detect under <outdir>/variant_calling/<tool>/.  Order matters only for
# deterministic iteration in tests.
_VARIANT_CALLERS: tuple[str, ...] = (
    "strelka",
    "mutect2",
    "freebayes",
    "deepvariant",
    "haplotypecaller",
    "sentieon_haplotyper",
    "sentieon_dnascope",
    "sentieon_tnscope",
    "lofreq",
    "muse",
    "mpileup",
    "manta",
    "tiddit",
    "ascat",
    "cnvkit",
    "controlfreec",
    "indexcov",
    "msisensor2",
    "msisensorpro",
    "varlociraptor",
)

# Some callers publish under a directory whose name differs from the --tools
# token. The rendered output docs call the MSIsensorPro directory `msisensor/`,
# while conf/modules/msisensorpro.config in tag 3.8.1 publishes `msisensorpro/`;
# accept both, preferring the executable configuration's path.
_TOOL_DIRS = {
    "mpileup": ("bcftools",),
    "msisensorpro": ("msisensorpro", "msisensor"),
}

# Some callers write VCFs under a suffix that differs from the directory name.
# (e.g. dir `sentieon_haplotyper/` → file `<sample>.haplotyper.*.vcf.gz`;
#  mpileup → `<sample>.bcftools.vcf.gz`).
_VCF_SUFFIX = {
    "sentieon_haplotyper": "haplotyper",
    "sentieon_dnascope": "dnascope",
    "sentieon_tnscope": "tnscope",
    "mpileup": "bcftools",
}

# Annotation tokens recognised by --tools that emit `*.ann.vcf.gz`.
# `bcfann`, `merge` and `snpsift` are official --tools values too.
_ANNOTATORS: tuple[str, ...] = ("snpeff", "vep", "bcfann", "merge", "snpsift")

# Suffix for the primary aligned BAM/CRAM file -- depends on --save_output_as_bam.
def _ext_and_idx(save_output_as_bam: bool) -> tuple[str, str]:
    return ("bam", "bai") if save_output_as_bam else ("cram", "crai")


# ---------------------------------------------------------------------------
# Public dataclass
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class OutputsReport:
    step_completed: str
    preprocessing: dict
    variant_calling: dict
    annotation: dict
    qc: dict
    pipeline_info: Path | None
    reference_outputs: dict | None
    samples_detected: list[str]
    pairs_detected: list[str]
    csv_handoff: dict[str, Path]
    handoff_available: bool
    warnings: list[str] = field(default_factory=list)
    missing_outputs: list[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Helpers (file/path)
# ---------------------------------------------------------------------------


def _nonempty(path: Path, warnings: list[str]) -> Path | None:
    """Return ``path`` only if it exists and is non-empty.

    Logs a warning for zero-byte files so we never silently return junk.
    """
    try:
        st = path.stat()
    except OSError:
        return None
    if not path.is_file():
        return None
    if st.st_size <= 0:
        warnings.append(f"Zero-byte output skipped: {path}")
        return None
    return path


def _filter_nonempty(paths: Iterable[Path], warnings: list[str]) -> list[Path]:
    keep: list[Path] = []
    for p in paths:
        kept = _nonempty(p, warnings)
        if kept is not None:
            keep.append(kept)
    return keep


def _first_nonempty(candidates: Iterable[Path], warnings: list[str]) -> Path | None:
    for c in candidates:
        ok = _nonempty(c, warnings)
        if ok is not None:
            return ok
    return None


def _vcf_tbi(vcf: Path, warnings: list[str]) -> Path | None:
    return _nonempty(vcf.with_suffix(vcf.suffix + ".tbi"), warnings)


def _sorted_subdirs(parent: Path) -> list[Path]:
    if not parent.is_dir():
        return []
    return sorted([p for p in parent.iterdir() if p.is_dir()])


# ---------------------------------------------------------------------------
# §1 Preprocessing
# ---------------------------------------------------------------------------


def _parse_preprocessing(
    output_dir: Path,
    *,
    save_mapped: bool,
    save_output_as_bam: bool,
    warnings: list[str],
) -> tuple[dict, dict[str, Path], list[str]]:
    pre_dir = output_dir / "preprocessing"
    ext, idx = _ext_and_idx(save_output_as_bam)

    # Sarek 3.8.1 publishes restart CSVs at top-level <outdir>/csv/ in the
    # usage docs and official rendered results tree.  The output docs also
    # mention <outdir>/preprocessing/csv/, so keep that as a compatibility
    # fallback for older/generated fixtures.
    csv_dirs = (output_dir / "csv", pre_dir / "csv")
    csv_handoff: dict[str, Path] = {}
    for key, fname in (
        ("mapped", "mapped.csv"),
        ("markduplicates_no_table", "markduplicates_no_table.csv"),
        ("markduplicates", "markduplicates.csv"),
        ("recalibrated", "recalibrated.csv"),
        ("variantcalled", "variantcalled.csv"),
    ):
        p = next((candidate for d in csv_dirs if (candidate := d / fname).is_file()), csv_dirs[0] / fname)
        kept = _nonempty(p, warnings)
        if kept is not None:
            csv_handoff[key] = kept

    # Per-sample CRAM/BAM under markduplicates/recalibrated/mapped
    md_dir = pre_dir / "markduplicates"
    rc_dir = pre_dir / "recalibrated"
    mp_dir = pre_dir / "mapped"

    markduplicates: dict[str, dict[str, Path]] = {}
    for sub in _sorted_subdirs(md_dir):
        sample = sub.name
        cram = _first_nonempty([sub / f"{sample}.md.{ext}"], warnings)
        crai = _first_nonempty([sub / f"{sample}.md.{ext}.{idx}"], warnings)
        entry: dict[str, Path] = {}
        if cram is not None:
            entry["cram" if ext == "cram" else "bam"] = cram
        if crai is not None:
            entry["crai" if idx == "crai" else "bai"] = crai
        if entry:
            markduplicates[sample] = entry

    recalibrated: dict[str, dict[str, Path]] = {}
    for sub in _sorted_subdirs(rc_dir):
        sample = sub.name
        cram = _first_nonempty([sub / f"{sample}.recal.{ext}"], warnings)
        crai = _first_nonempty([sub / f"{sample}.recal.{ext}.{idx}"], warnings)
        entry = {}
        if cram is not None:
            entry["cram" if ext == "cram" else "bam"] = cram
        if crai is not None:
            entry["crai" if idx == "crai" else "bai"] = crai
        if entry:
            recalibrated[sample] = entry

    mapped: dict[str, dict[str, Path]] = {}
    if save_mapped:
        for sub in _sorted_subdirs(mp_dir):
            sample = sub.name
            cram = _first_nonempty([sub / f"{sample}.sorted.{ext}"], warnings)
            crai = _first_nonempty([sub / f"{sample}.sorted.{ext}.{idx}"], warnings)
            entry = {}
            if cram is not None:
                entry["cram" if ext == "cram" else "bam"] = cram
            if crai is not None:
                entry["crai" if idx == "crai" else "bai"] = crai
            if entry:
                mapped[sample] = entry

    # BQSR recalibration tables: preprocessing/recal_table/<sample>/<sample>.recal.table
    recal_table: dict[str, Path] = {}
    for sub in _sorted_subdirs(pre_dir / "recal_table"):
        sample = sub.name
        table = _first_nonempty([sub / f"{sample}.recal.table"], warnings)
        if table is not None:
            recal_table[sample] = table

    # Optional preprocessing artefacts (conditional on the relevant save/UMI/
    # BBSplit/Sentieon/aligner flags) - catalogued as file lists when present.
    extras = {
        name: _glob_files(pre_dir / name, "*", warnings)
        for name in (
            "fastp",
            "umi",
            "bbsplit",
            "sentieon_dedup",
            # conf/modules/sentieon_dedup.config switches to this directory
            # when --sentieon_consensus is enabled.
            "sentieon_consensus",
            # Documented output directory for --aligner parabricks.
            "parabricks",
            # Official module configs publish BAM/CRAM conversions here,
            # including variant-calling CRAM_TO_BAM restart conversion.
            "converted",
        )
    }

    # Sample IDs from CSVs (preferred), falling back to path scan
    samples_from_csvs = _samples_from_csvs(csv_handoff, warnings)

    preprocessing = {
        "csv_handoff": dict(csv_handoff),
        "markduplicates": markduplicates,
        "recalibrated": recalibrated,
        "recal_table": recal_table,
        "mapped": mapped if save_mapped else {},
        "fastp": extras["fastp"],
        "umi": extras["umi"],
        "bbsplit": extras["bbsplit"],
        "sentieon_dedup": extras["sentieon_dedup"],
        "sentieon_consensus": extras["sentieon_consensus"],
        "parabricks": extras["parabricks"],
        "converted": extras["converted"],
        "samples_detected_from_csvs": samples_from_csvs,
    }
    return preprocessing, csv_handoff, samples_from_csvs


def _samples_from_csvs(csv_handoff: dict[str, Path], warnings: list[str]) -> list[str]:
    """Return unique sample IDs gathered from the most-recent CSV available.

    Preference order follows pipeline progress, newest first: variantcalled >
    recalibrated > markduplicates > markduplicates_no_table > mapped.
    """
    for key in ("variantcalled", "recalibrated", "markduplicates", "markduplicates_no_table", "mapped"):
        csv_path = csv_handoff.get(key)
        if csv_path is None:
            continue
        try:
            with csv_path.open(newline="", encoding="utf-8") as fh:
                reader = csv.DictReader(fh)
                if reader.fieldnames is None or "sample" not in reader.fieldnames:
                    warnings.append(f"CSV without 'sample' column: {csv_path}")
                    continue
                seen: list[str] = []
                seen_set: set[str] = set()
                for row in reader:
                    s = (row.get("sample") or "").strip()
                    if s and s not in seen_set:
                        seen.append(s)
                        seen_set.add(s)
                return seen
        except OSError as exc:
            warnings.append(f"Could not read CSV {csv_path}: {exc}")
            continue
    return []


# ---------------------------------------------------------------------------
# §2 Variant calling
# ---------------------------------------------------------------------------


def _is_pair_dir(name: str) -> bool:
    return bool(_PAIR_RE.match(name))


def _parse_strelka(tool_dir: Path, warnings: list[str]) -> dict[str, dict]:
    out: dict[str, dict] = {}
    for sub in _sorted_subdirs(tool_dir):
        key = sub.name
        entry: dict[str, object] = {}
        # Germline
        genome_vcf = _nonempty(sub / f"{key}.strelka.genome.vcf.gz", warnings)
        germline_vcf = _nonempty(sub / f"{key}.strelka.variants.vcf.gz", warnings)
        if genome_vcf is not None:
            entry["genome_vcf"] = genome_vcf
            if tbi := _vcf_tbi(genome_vcf, warnings):
                entry["genome_tbi"] = tbi
        if germline_vcf is not None:
            entry["vcf"] = [germline_vcf]
            if tbi := _vcf_tbi(germline_vcf, warnings):
                entry["tbi"] = [tbi]
        # Somatic
        indels_vcf = _nonempty(sub / f"{key}.strelka.somatic_indels.vcf.gz", warnings)
        snvs_vcf = _nonempty(sub / f"{key}.strelka.somatic_snvs.vcf.gz", warnings)
        if indels_vcf is not None:
            entry["indels_vcf"] = indels_vcf
            if tbi := _vcf_tbi(indels_vcf, warnings):
                entry["indels_tbi"] = tbi
        if snvs_vcf is not None:
            entry["snvs_vcf"] = snvs_vcf
            if tbi := _vcf_tbi(snvs_vcf, warnings):
                entry["snvs_tbi"] = tbi
        if entry:
            out[key] = entry
    return out


def _parse_mutect2(tool_dir: Path, warnings: list[str]) -> dict[str, dict]:
    """Mutect2 emits an unfiltered VCF, a filtered VCF, and several auxiliary
    tables/archives used by FilterMutectCalls — capture them all.
    See https://nf-co.re/sarek/3.8.1/docs/output/#mutect2 .
    """
    out: dict[str, dict] = {}
    for sub in _sorted_subdirs(tool_dir):
        key = sub.name
        entry: dict[str, object] = {}
        candidates = (
            ("vcf", sub / f"{key}.mutect2.vcf.gz"),
            ("vcf_stats", sub / f"{key}.mutect2.vcf.gz.stats"),
            ("filtered_vcf", sub / f"{key}.mutect2.filtered.vcf.gz"),
            ("filtering_stats", sub / f"{key}.mutect2.filtered.vcf.gz.filteringStats.tsv"),
            ("contamination_table", sub / f"{key}.mutect2.contamination.table"),
            ("segmentation_table", sub / f"{key}.mutect2.segmentation.table"),
            ("artifact_prior", sub / f"{key}.mutect2.artifactprior.tar.gz"),
        )
        for label, path in candidates:
            kept = _nonempty(path, warnings)
            if kept is not None:
                entry[label] = kept
                if label.endswith("vcf") and (tbi := _vcf_tbi(kept, warnings)):
                    entry[f"{label}_tbi"] = tbi
        pileups = _filter_nonempty(sorted(sub.glob("*.mutect2.pileups.table")), warnings)
        if pileups:
            entry["pileups_table"] = pileups[0]
            entry["pileups_tables"] = pileups
        if entry:
            out[key] = entry
    return out


def _parse_simple_single(tool: str, tool_dir: Path, warnings: list[str]) -> dict[str, dict]:
    """Generic ``<sample>/<sample>.<suffix>.vcf.gz`` layout.

    ``suffix`` defaults to the tool name (freebayes, deepvariant, …) but the
    Sentieon callers publish under a shortened suffix (haplotyper/dnascope/
    tnscope). Preserve raw, filtered and unfiltered products when more than
    one is published, while exposing the most processed output as ``vcf``.
    """
    suffix = _VCF_SUFFIX.get(tool, tool)
    out: dict[str, dict] = {}
    for sub in _sorted_subdirs(tool_dir):
        key = sub.name
        filtered = _nonempty(sub / f"{key}.{suffix}.filtered.vcf.gz", warnings)
        unfiltered = _nonempty(sub / f"{key}.{suffix}.unfiltered.vcf.gz", warnings)
        raw = _nonempty(sub / f"{key}.{suffix}.vcf.gz", warnings)
        if filtered is not None:
            out.setdefault(key, {})["filtered_vcf"] = filtered
            if tbi := _vcf_tbi(filtered, warnings):
                out[key]["filtered_tbi"] = tbi
        if unfiltered is not None:
            out.setdefault(key, {})["unfiltered_vcf"] = unfiltered
            if tbi := _vcf_tbi(unfiltered, warnings):
                out[key]["unfiltered_tbi"] = tbi
        if raw is not None:
            out.setdefault(key, {})["raw_vcf"] = raw
            if tbi := _vcf_tbi(raw, warnings):
                out[key]["raw_tbi"] = tbi
        # Sentieon Haplotyper publishes the pre-filter merge as
        # `<id>.haplotyper.unfiltered.vcf.gz` and the FilterVariantTranches
        # product as `<id>.haplotyper.vcf.gz`; prefer the base/final callset.
        vcf = filtered or raw or unfiltered
        if vcf is None:
            # The rendered output guide documents LoFreq as `<sample>.vcf.gz`,
            # while conf/modules/lofreq.config prefixes `<sample>.lofreq`.
            # Accept the documented form as well as the executable config form
            # handled by `raw` above.
            fallback = [
                p for p in sorted(sub.glob("*.vcf.gz"))
                if not p.name.endswith(".g.vcf.gz")
            ]
            vcf = _first_nonempty(fallback, warnings)
        if vcf is not None:
            out.setdefault(key, {})["vcf"] = vcf
            if tbi := _vcf_tbi(vcf, warnings):
                out[key]["tbi"] = tbi
        gvcf = _nonempty(sub / f"{key}.{suffix}.g.vcf.gz", warnings)
        if gvcf is not None:
            out.setdefault(key, {})["gvcf"] = gvcf
            if tbi := _vcf_tbi(gvcf, warnings):
                out.setdefault(key, {})["gvcf_tbi"] = tbi
        if tool == "muse":
            # MuSE publishes the score/call table alongside its compressed
            # VCF (docs/output.md: <tumorsample_vs_normalsample>.MuSE.txt).
            call_table = _nonempty(sub / f"{key}.MuSE.txt", warnings)
            if call_table is not None:
                out.setdefault(key, {})["call_table"] = call_table
    return out


def _parse_haplotypecaller(tool_dir: Path, warnings: list[str]) -> dict[str, dict]:
    out: dict[str, dict] = {}
    for sub in _sorted_subdirs(tool_dir):
        key = sub.name
        entry: dict[str, object] = {}
        vcf = _nonempty(sub / f"{key}.haplotypecaller.vcf.gz", warnings)
        filtered = _nonempty(sub / f"{key}.haplotypecaller.filtered.vcf.gz", warnings)
        gvcf = _nonempty(sub / f"{key}.haplotypecaller.g.vcf.gz", warnings)
        if vcf is not None:
            entry["vcf"] = vcf
            if tbi := _vcf_tbi(vcf, warnings):
                entry["tbi"] = tbi
        if filtered is not None:
            entry["filtered_vcf"] = filtered
            if tbi := _vcf_tbi(filtered, warnings):
                entry["filtered_tbi"] = tbi
        if gvcf is not None:
            entry["gvcf"] = gvcf
            if tbi := _vcf_tbi(gvcf, warnings):
                entry["gvcf_tbi"] = tbi
        if entry:
            out[key] = entry
    return out


def _parse_manta(tool_dir: Path, warnings: list[str]) -> dict[str, dict]:
    out: dict[str, dict] = {}
    for sub in _sorted_subdirs(tool_dir):
        key = sub.name
        entry: dict[str, object] = {}
        diploid = _nonempty(sub / f"{key}.manta.diploid_sv.vcf.gz", warnings)
        small_indels = _nonempty(sub / f"{key}.manta.candidate_small_indels.vcf.gz", warnings)
        candidate_sv = _nonempty(sub / f"{key}.manta.candidate_sv.vcf.gz", warnings)
        somatic = _nonempty(sub / f"{key}.manta.somatic_sv.vcf.gz", warnings)
        tumor = _nonempty(sub / f"{key}.manta.tumor_sv.vcf.gz", warnings)
        if diploid is not None:
            entry["sv_vcf"] = diploid
            if tbi := _vcf_tbi(diploid, warnings):
                entry["sv_tbi"] = tbi
        if somatic is not None:
            # Somatic key takes precedence for the canonical ``sv_vcf`` slot when
            # only a pair was called; we still keep diploid (germline) separately
            # if both happen to be present.
            entry.setdefault("sv_vcf", somatic)
            entry["somatic_sv_vcf"] = somatic
            if tbi := _vcf_tbi(somatic, warnings):
                entry["somatic_sv_tbi"] = tbi
        if tumor is not None:
            entry.setdefault("sv_vcf", tumor)
            entry["tumor_sv_vcf"] = tumor
            if tbi := _vcf_tbi(tumor, warnings):
                entry["tumor_sv_tbi"] = tbi
        if small_indels is not None:
            entry["small_indels_vcf"] = small_indels
        if candidate_sv is not None:
            entry["candidate_sv_vcf"] = candidate_sv
        if entry:
            out[key] = entry
    return out


def _parse_tiddit(tool_dir: Path, warnings: list[str]) -> dict[str, dict]:
    out: dict[str, dict] = {}
    for sub in _sorted_subdirs(tool_dir):
        key = sub.name
        entry: dict[str, object] = {}
        single = _nonempty(sub / f"{key}.tiddit.vcf.gz", warnings)
        normal = _nonempty(sub / f"{key}.tiddit.normal.vcf.gz", warnings)
        tumor = _nonempty(sub / f"{key}.tiddit.tumor.vcf.gz", warnings)
        # docs/output.md and the executable 3.8.1 config disagree on the
        # merged paired-call prefix; support both official representations.
        merged = _first_nonempty(
            (
                sub / f"{key}.tiddit_sv_merge.vcf.gz",
                sub / f"{key}_sv_merge.tiddit.vcf.gz",
            ),
            warnings,
        )
        ploidies = _nonempty(sub / f"{key}.tiddit.ploidies.tab", warnings)
        if single is not None:
            entry["vcf"] = single
            if tbi := _vcf_tbi(single, warnings):
                entry["tbi"] = tbi
        if normal is not None:
            entry["normal_vcf"] = normal
            if tbi := _vcf_tbi(normal, warnings):
                entry["normal_tbi"] = tbi
        if tumor is not None:
            entry["tumor_vcf"] = tumor
            if tbi := _vcf_tbi(tumor, warnings):
                entry["tumor_tbi"] = tbi
        if merged is not None:
            entry["merged_vcf"] = merged
            if tbi := _vcf_tbi(merged, warnings):
                entry["merged_tbi"] = tbi
        if ploidies is not None:
            entry["ploidies"] = ploidies
        if entry:
            out[key] = entry
    return out


def _parse_post_variant(output_dir: Path, warnings: list[str]) -> dict[str, dict]:
    out: dict[str, dict] = {}
    root = output_dir / "variant_calling"
    for name in ("filtered", "normalized", "consensus", "concat"):
        parsed = _parse_generic_dir(name, root / name, warnings)
        if parsed:
            out[name] = parsed
    return out


def _parse_joint_variant_calling(output_dir: Path, warnings: list[str]) -> dict[str, dict]:
    out: dict[str, dict] = {}
    vc_root = output_dir / "variant_calling"
    for caller in ("haplotypecaller", "sentieon_dnascope", "sentieon_haplotyper"):
        parsed = _parse_joint_dir(vc_root / caller / "joint_variant_calling", warnings)
        if parsed:
            out[f"joint_{caller}"] = parsed
    legacy_root = output_dir / "joint_variant_calling"
    for caller_dir in _sorted_subdirs(legacy_root):
        parsed = _parse_joint_dir(caller_dir, warnings)
        if parsed:
            out[f"joint_{caller_dir.name}"] = parsed
    return out


def _parse_joint_dir(joint_dir: Path, warnings: list[str]) -> dict[str, dict]:
    out: dict[str, dict] = {}
    if not joint_dir.is_dir():
        return out
    for vcf in _filter_nonempty(sorted(joint_dir.glob("*.vcf.gz")), warnings):
        name = vcf.name
        key = name.removesuffix(".vcf.gz")
        entry: dict[str, object] = {"vcf": vcf}
        tbi = _nonempty(vcf.with_suffix(vcf.suffix + ".tbi"), warnings)
        if tbi is not None:
            entry["tbi"] = tbi
        out[key] = entry
    for sub in _sorted_subdirs(joint_dir):
        parsed = _parse_generic_dir(sub.name, sub, warnings)
        if parsed:
            out[sub.name] = parsed.get(sub.name, parsed)
    return out


def _parse_ascat(tool_dir: Path, warnings: list[str]) -> dict[str, dict]:
    out: dict[str, dict] = {}
    for sub in _sorted_subdirs(tool_dir):
        key = sub.name
        entry: dict[str, object] = {}
        # ASCAT writes <name>.cnvs.txt etc. where <name> typically matches the dir.
        cnvs_candidates = sorted(sub.glob("*.cnvs.txt"))
        purity_candidates = sorted(sub.glob("*.purityploidy.txt"))
        segments_candidates = sorted(sub.glob("*.segments.txt"))
        metrics_candidates = sorted(sub.glob("*.metrics.txt"))
        # ASCAT also publishes raw tumour BAF / LogR txt traces.
        # See https://nf-co.re/sarek/3.8.1/docs/output/#ascat
        tumour_baf_candidates = sorted(sub.glob("*_tumourBAF.txt"))
        tumour_logr_candidates = sorted(sub.glob("*_tumourLogR.txt"))
        normal_baf_candidates = sorted(sub.glob("*_normalBAF.txt"))
        normal_logr_candidates = sorted(sub.glob("*_normalLogR.txt"))
        plot_candidates = sorted(sub.glob("*.png"))

        cnvs = _first_nonempty(cnvs_candidates, warnings)
        purity = _first_nonempty(purity_candidates, warnings)
        segments = _first_nonempty(segments_candidates, warnings)
        metrics = _first_nonempty(metrics_candidates, warnings)
        tumour_baf = _first_nonempty(tumour_baf_candidates, warnings)
        tumour_logr = _first_nonempty(tumour_logr_candidates, warnings)
        normal_baf = _first_nonempty(normal_baf_candidates, warnings)
        normal_logr = _first_nonempty(normal_logr_candidates, warnings)
        plots = _filter_nonempty(plot_candidates, warnings)

        if cnvs is not None:
            entry["cnvs"] = cnvs
        if purity is not None:
            entry["purity_ploidy"] = purity
        if segments is not None:
            entry["segments"] = segments
        if metrics is not None:
            entry["metrics"] = metrics
        if tumour_baf is not None:
            entry["tumour_baf"] = tumour_baf
        if tumour_logr is not None:
            entry["tumour_logr"] = tumour_logr
        if normal_baf is not None:
            entry["normal_baf"] = normal_baf
        if normal_logr is not None:
            entry["normal_logr"] = normal_logr
        if plots:
            entry["plots"] = plots
        if entry:
            out[key] = entry
    return out


def _parse_cnvkit(tool_dir: Path, warnings: list[str]) -> dict[str, dict]:
    out: dict[str, dict] = {}
    for sub in _sorted_subdirs(tool_dir):
        key = sub.name
        entry: dict[str, object] = {}
        # Files are prefixed by the sample name; for tumor/normal pairs the dir
        # is <tumor_vs_normal> but the files use the individual sample names,
        # so we glob rather than assume the directory name is the prefix.
        cnr = _first_nonempty(sorted(sub.glob("*.cnr")), warnings)
        call_cns = _first_nonempty(sorted(sub.glob("*.call.cns")), warnings)
        bintest = _first_nonempty(sorted(sub.glob("*.bintest.cns")), warnings)
        coverage = _filter_nonempty(sorted(sub.glob("*.cnn")), warnings)
        gene_metrics = _filter_nonempty(sorted(sub.glob("*.genemetrics.tsv")), warnings)
        diagrams = _filter_nonempty(sorted(sub.glob("*-diagram.pdf")), warnings)
        scatters = _filter_nonempty(sorted(sub.glob("*-scatter.png")), warnings)
        exported_vcf = _first_nonempty(sorted(sub.glob("*.vcf")), warnings)
        beds = _filter_nonempty(sorted(sub.glob("*.bed")), warnings)
        # Plain .cns must exclude the .call.cns / .bintest.cns variants.
        plain_cns = _first_nonempty(
            [
                p for p in sorted(sub.glob("*.cns"))
                if not p.name.endswith((".call.cns", ".bintest.cns"))
            ],
            warnings,
        )
        if cnr is not None:
            entry["cnr"] = cnr
        if plain_cns is not None:
            entry["cns"] = plain_cns
        if call_cns is not None:
            entry["call_cns"] = call_cns
        if bintest is not None:
            entry["bintest_cns"] = bintest
        if coverage:
            entry["coverage"] = coverage
        if gene_metrics:
            entry["gene_metrics"] = gene_metrics
        if diagrams:
            entry["diagrams"] = diagrams
        if scatters:
            entry["scatters"] = scatters
        if exported_vcf is not None:
            entry["vcf"] = exported_vcf
        if beds:
            entry["bed"] = beds
        if entry:
            out[key] = entry
    return out


def _parse_msisensor(tool_dir: Path, warnings: list[str]) -> dict[str, dict]:
    out: dict[str, dict] = {}
    for sub in _sorted_subdirs(tool_dir):
        key = sub.name
        entry: dict[str, object] = {}
        score = _nonempty(sub / key, warnings)
        distribution = _nonempty(sub / f"{key}_dis", warnings)
        germline = _nonempty(sub / f"{key}_germline", warnings)
        somatic = _nonempty(sub / f"{key}_somatic", warnings)
        if score is not None:
            entry["score"] = score
        if distribution is not None:
            entry["distribution"] = distribution
        if germline is not None:
            entry["germline_sites"] = germline
        if somatic is not None:
            entry["somatic_sites"] = somatic
        if entry:
            out[key] = entry
    return out


def _parse_generic_dir(tool: str, tool_dir: Path, warnings: list[str]) -> dict[str, dict]:
    """Fallback: catalog files per child dir, plus files written directly under
    the tool dir (e.g. indexcov emits ``<sample>-indexcov.bed.gz`` at the root)."""
    out: dict[str, dict] = {}

    def _collect(scope: Path) -> dict[str, object]:
        entry: dict[str, object] = {}
        vcfs = _filter_nonempty(sorted(scope.glob("*.vcf.gz")), warnings)
        tbis = _filter_nonempty(sorted(scope.glob("*.vcf.gz.tbi")), warnings)
        txts = _filter_nonempty(sorted(scope.glob("*.txt")), warnings)
        tsvs = _filter_nonempty(sorted(scope.glob("*.tsv")), warnings)
        tabs = _filter_nonempty(sorted(scope.glob("*.tab")), warnings)
        jsons = _filter_nonempty(sorted(scope.glob("*.json")), warnings)
        yamls = _filter_nonempty(sorted(scope.glob("*.yml")) + sorted(scope.glob("*.yaml")), warnings)
        htmls = _filter_nonempty(sorted(scope.glob("*.html")), warnings)
        logs = _filter_nonempty(sorted(scope.glob("*.log")), warnings)
        metrics = _filter_nonempty(sorted(scope.glob("*.metrics")), warnings)
        pngs = _filter_nonempty(sorted(scope.glob("*.png")), warnings)
        beds = _filter_nonempty(sorted(scope.glob("*.bed.gz")), warnings)
        files = _filter_nonempty([p for p in sorted(scope.glob("*")) if p.is_file()], warnings)
        if vcfs:
            entry["vcf"] = vcfs
        if tbis:
            entry["tbi"] = tbis
        if txts:
            entry["txt"] = txts
        if tsvs:
            entry["tsv"] = tsvs
        if tabs:
            entry["tab"] = tabs
        if jsons:
            entry["json"] = jsons
        if yamls:
            entry["yaml"] = yamls
        if htmls:
            entry["html"] = htmls
        if logs:
            entry["log"] = logs
        if metrics:
            entry["metrics"] = metrics
        if pngs:
            entry["plots"] = pngs
        if beds:
            entry["bed"] = beds
        if files:
            entry["files"] = files
        return entry

    for sub in _sorted_subdirs(tool_dir):
        entry = _collect(sub)
        if entry:
            out[sub.name] = entry

    # Files placed directly in the tool dir (no per-sample subdir), e.g. indexcov.
    root_files = [p for p in tool_dir.glob("*") if p.is_file()] if tool_dir.is_dir() else []
    if root_files:
        root_entry = _collect(tool_dir)
        if root_entry:
            # These are cohort-level files, not a biological sample named
            # after the tool. Keep a synthetic key that sample inference can
            # explicitly ignore.
            out.setdefault("_cohort", root_entry)
    return out


_PARSERS = {
    "strelka": _parse_strelka,
    "mutect2": _parse_mutect2,
    "haplotypecaller": _parse_haplotypecaller,
    "manta": _parse_manta,
    "tiddit": lambda tool_dir, warnings: _parse_tiddit(tool_dir, warnings),
    "ascat": _parse_ascat,
    "cnvkit": _parse_cnvkit,
    "msisensor2": _parse_msisensor,
    "msisensorpro": _parse_msisensor,
}


def _parse_variant_calling(
    output_dir: Path,
    *,
    tools_lower: set[str],
    warnings: list[str],
) -> tuple[dict, list[str]]:
    vc_root = output_dir / "variant_calling"
    variant_calling: dict[str, dict] = {}
    pairs: list[str] = []
    pair_set: set[str] = set()

    if not vc_root.is_dir():
        joint = _parse_joint_variant_calling(output_dir, warnings)
        if joint:
            variant_calling.update(joint)
        return variant_calling, pairs

    for tool in _VARIANT_CALLERS:
        parsed_for_tool: dict[str, dict] = {}
        for dir_name in _TOOL_DIRS.get(tool, (tool,)):
            tool_dir = vc_root / dir_name
            if not tool_dir.is_dir():
                continue
            parser = _PARSERS.get(tool)
            if parser is not None:
                parsed = parser(tool_dir, warnings)
            elif tool in {
                "freebayes", "deepvariant", "lofreq", "muse", "mpileup",
                "sentieon_haplotyper", "sentieon_dnascope", "sentieon_tnscope",
            }:
                parsed = _parse_simple_single(tool, tool_dir, warnings)
            else:
                # controlfreec, indexcov, msisensor2, msisensorpro, varlociraptor
                parsed = _parse_generic_dir(tool, tool_dir, warnings)

            for key, entry in parsed.items():
                parsed_for_tool.setdefault(key, entry)
            # Pairs come from any tool's child dir names.
            for sub in _sorted_subdirs(tool_dir):
                if _is_pair_dir(sub.name) and sub.name not in pair_set:
                    pair_set.add(sub.name)
                    pairs.append(sub.name)
        if parsed_for_tool:
            variant_calling[tool] = parsed_for_tool

    pairs.sort()
    variant_calling.update(_parse_joint_variant_calling(output_dir, warnings))
    variant_calling.update(_parse_post_variant(output_dir, warnings))
    return variant_calling, pairs


# ---------------------------------------------------------------------------
# §3 Annotation
# ---------------------------------------------------------------------------


# Map annotation token (official --tools value) -> output filename suffix, as
# emitted by conf/modules/annotate.config in nf-core/sarek 3.8.1. `merge` runs
# snpEff THEN VEP and is published as `_snpEff_VEP`; bcfann is published `_BCF`.
# Suffixes are matched longest-first so `_snpEff_VEP` is not mis-detected as
# `_VEP` (or `_snpEff`).
_ANN_SUFFIXES = {
    "merge": "_snpEff_VEP",
    "snpeff": "_snpEff",
    "vep": "_VEP",
    "bcfann": "_BCF",
    "snpsift": "_snpSift",
}
_ANN_SUFFIX_ALIASES = {
    # The rendered output docs show `_bcf.ann.vcf.gz`, while
    # conf/modules/annotate.config in 3.8.1 uses `_BCF.ann`; accept both.
    "bcfann": ("_BCF", "_bcf"),
}

# Annotated-output extensions. VEP/merge honour --vep_out_format (vcf|json|tab);
# all published variants are gzipped (annotate.config pattern `*{gz,tbi}`).
_ANN_DATA_EXTS = (".ann.vcf.gz", ".ann.json.gz", ".ann.tab.gz")
_ANN_DATA_KEY = {
    ".ann.vcf.gz": "vcf",
    ".ann.json.gz": "json",
    ".ann.tab.gz": "tab",
}


def _annotation_report_dirs(output_dir: Path, ann_root: Path, data_file: Path, report_kind: str) -> list[Path]:
    """Return candidate report directories for an annotated output.

    The executable 3.8.1 config publishes ``reports/<tool>/<variantcaller>/<id>/``.
    The rendered output page describes ``reports/<tool>/<id>/<variantcaller>/``;
    both are official sources, so resolve both plus the legacy co-located layout.
    """
    dirs = [data_file.parent]
    try:
        relative = data_file.relative_to(ann_root)
    except ValueError:
        return dirs
    if len(relative.parts) < 3:
        return dirs
    variantcaller, sample_id = relative.parts[-3], relative.parts[-2]
    report_roots = ("snpeff", "SnpEff") if report_kind == "snpeff" else ("EnsemblVEP",)
    for root_name in report_roots:
        root = output_dir / "reports" / root_name
        dirs.extend((root / variantcaller / sample_id, root / sample_id / variantcaller))
    return dirs


def _first_report_html(dirs: Iterable[Path], patterns: Iterable[str], warnings: list[str]) -> Path | None:
    candidates: list[Path] = []
    for directory in dirs:
        for pattern in patterns:
            candidates.extend(sorted(directory.glob(pattern)))
    return _first_nonempty(candidates, warnings)


def _parse_annotation(output_dir: Path, warnings: list[str]) -> dict:
    ann_root = output_dir / "annotation"
    out: dict[str, dict] = {k: {} for k in _ANNOTATORS}
    if not ann_root.is_dir():
        return out

    # Sarek 3.8.1 publishes annotation as annotation/<variantcaller>/<sample_or_pair>/
    # (conf/modules/annotate.config), i.e. two levels deep. Recurse so both that
    # layout and any flat one are discovered.
    candidates: list[tuple[Path, str]] = []
    seen: set = set()
    for data_ext in _ANN_DATA_EXTS:
        for p in sorted(ann_root.rglob(f"*{data_ext}")):
            if p not in seen:
                seen.add(p)
                candidates.append((p, data_ext))
    for data_file, data_ext in candidates:
        kept = _nonempty(data_file, warnings)
        if kept is None:
            continue
        stem = data_file.name[: -len(data_ext)]
        kind = None
        base_stem = None
        # Longest suffix first so `_snpEff_VEP` is not matched as `_VEP`/`_snpEff`.
        suffix_candidates = []
        for candidate_kind, primary_suffix in _ANN_SUFFIXES.items():
            for candidate_suffix in _ANN_SUFFIX_ALIASES.get(candidate_kind, (primary_suffix,)):
                suffix_candidates.append((candidate_kind, candidate_suffix))
        for k, suf in sorted(suffix_candidates, key=lambda kv: -len(kv[1])):
            if stem.endswith(suf):
                kind = k
                base_stem = stem[: -len(suf)]
                break
        if kind is None:
            continue
        key = base_stem
        data_key = _ANN_DATA_KEY[data_ext]
        entry: dict[str, object] = {data_key: kept}
        if data_key == "vcf" and (tbi := _nonempty(data_file.with_suffix(data_file.suffix + ".tbi"), warnings)) is not None:
            entry["tbi"] = tbi
        if kind == "snpeff":
            entry["html"] = _first_report_html(
                _annotation_report_dirs(output_dir, ann_root, data_file, "snpeff"),
                (
                    # sarek 3.8.1's snpEff module passes `-csvStats <prefix>.csv`
                    # but no `-stats <prefix>`, so the stats HTML keeps snpEff's
                    # default constant name. The per-sample report dir resolved
                    # above scopes this so it never leaks across samples.
                    "snpEff_summary.html",
                    f"{base_stem}{_ANN_SUFFIXES['snpeff']}.html",
                    f"{base_stem}{_ANN_SUFFIXES['snpeff']}.ann.summary.html",
                ),
                warnings,
            )
        elif kind in {"vep", "merge"}:
            suffix = _ANN_SUFFIXES[kind]
            entry["html"] = _first_report_html(
                _annotation_report_dirs(output_dir, ann_root, data_file, "vep"),
                (
                    f"{base_stem}{suffix}.ann.summary.html",
                    f"{base_stem}{suffix}.summary.html",
                ),
                warnings,
            )
        out[kind][key] = entry
    return out


# ---------------------------------------------------------------------------
# §4 QC reports
# ---------------------------------------------------------------------------


def _glob_files(root: Path, pattern: str, warnings: list[str]) -> list[Path]:
    if not root.is_dir():
        return []
    matches = sorted(root.rglob(pattern))
    return _filter_nonempty(matches, warnings)


def _parse_qc(output_dir: Path, warnings: list[str]) -> dict:
    reports = output_dir / "reports"
    # MultiQC is published at the top level of outdir (multiqc/multiqc_report.html),
    # NOT under reports/. See the official nf-core/sarek 3.8.1 output spec.
    multiqc_report: Path | None = None
    mq_dir = output_dir / "multiqc"
    candidate = mq_dir / "multiqc_report.html"
    multiqc_report = _nonempty(candidate, warnings)
    if multiqc_report is None and mq_dir.is_dir():
        matches = sorted(mq_dir.rglob("multiqc_report.html"))
        multiqc_report = _first_nonempty(matches, warnings)
    if multiqc_report is None:
        # Fallback: a stray top-level multiqc_report.html.
        multiqc_report = _nonempty(output_dir / "multiqc_report.html", warnings)

    ngscheckmate: Path | None = None
    nc_dir = reports / "ngscheckmate"
    if nc_dir.is_dir():
        # NGSCheckMate writes a couple of summary files; we don't enforce a name.
        matches = sorted(nc_dir.rglob("*"))
        # Return the directory itself as a marker, but pick first file as canonical.
        files = [p for p in matches if p.is_file() and _nonempty(p, warnings) is not None]
        ngscheckmate = files[0] if files else None

    return {
        "multiqc_report": multiqc_report,
        "fastqc": _glob_files(reports / "fastqc", "*", warnings),
        "fastp": _glob_files(reports / "fastp", "*", warnings),
        "mosdepth": _glob_files(reports / "mosdepth", "*", warnings),
        "markduplicates_metrics": _glob_files(reports / "markduplicates", "*", warnings),
        "sentieon_dedup_metrics": _glob_files(reports / "sentieon_dedup", "*", warnings),
        "samtools_stats": _glob_files(reports / "samtools", "*", warnings),
        "bcftools_stats": _glob_files(reports / "bcftools", "*", warnings),
        "vcftools": _glob_files(reports / "vcftools", "*", warnings),
        # Per conf/modules/annotate.config: snpEff reports go to reports/snpeff/
        # (lowercase) and VEP to reports/EnsemblVEP/ (mixed case). Accept the
        # capitalised snpeff variant too for forward/backward robustness.
        "snpeff_reports": (
            _glob_files(reports / "snpeff", "*", warnings)
            or _glob_files(reports / "SnpEff", "*", warnings)
        ),
        "vep_reports": _glob_files(reports / "EnsemblVEP", "*", warnings),
        "ngscheckmate": ngscheckmate,
        "umi": _glob_files(reports / "umi", "*", warnings),
    }


# ---------------------------------------------------------------------------
# §6 Reference outputs (only when --save_reference)
# ---------------------------------------------------------------------------


def _parse_reference(output_dir: Path, warnings: list[str]) -> dict:
    ref_root = output_dir / "reference"
    if not ref_root.is_dir():
        return {}
    out: dict[str, list[Path]] = {}
    root_files = _filter_nonempty(sorted(p for p in ref_root.iterdir() if p.is_file()), warnings)
    if root_files:
        out["root"] = root_files
    for sub in _sorted_subdirs(ref_root):
        files = _filter_nonempty(sorted(p for p in sub.rglob("*") if p.is_file()), warnings)
        if files:
            out[sub.name] = files
    return out


# ---------------------------------------------------------------------------
# Sample/pair detection (fallbacks)
# ---------------------------------------------------------------------------


def _detect_samples(
    csv_samples: list[str],
    preprocessing: dict,
    variant_calling: dict,
) -> list[str]:
    if csv_samples:
        return list(csv_samples)
    seen: list[str] = []
    seen_set: set[str] = set()
    # markduplicates path scan
    for sample in preprocessing.get("markduplicates", {}):
        if sample not in seen_set:
            seen_set.add(sample)
            seen.append(sample)
    for sample in preprocessing.get("recalibrated", {}):
        if sample not in seen_set:
            seen_set.add(sample)
            seen.append(sample)
    if seen:
        return seen
    # Fallback: variant_calling parent dirs that are NOT pair dirs
    for tool, by_key in variant_calling.items():
        for key in by_key:
            if key.startswith("_"):
                continue
            if not _is_pair_dir(key) and key not in seen_set:
                seen_set.add(key)
                seen.append(key)
    return seen


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------


def _validate_outputs(
    *,
    step: str,
    tools_lower: set[str],
    variant_calling: dict,
    annotation: dict,
    output_dir: Path,
) -> list[str]:
    """Return ``missing_outputs`` list and raise on hard validation failures."""
    missing: list[str] = []
    requested_callers = tools_lower & set(_VARIANT_CALLERS)

    step_rank = {name: idx for idx, name in enumerate(_STEP_ORDER)}
    starts_before_or_at_variant_calling = step_rank.get(step, 999) <= step_rank["variant_calling"]
    starts_before_or_at_annotate = step_rank.get(step, 999) <= step_rank["annotate"]

    if starts_before_or_at_variant_calling and requested_callers:
        any_found = False
        for tool in requested_callers:
            if (tool in variant_calling and variant_calling[tool]) or variant_calling.get(f"joint_{tool}"):
                any_found = True
            else:
                missing.append(tool)
        if not any_found:
            vc_root = output_dir / "variant_calling"
            found_glob = []
            if vc_root.is_dir():
                found_glob = sorted(str(p.relative_to(output_dir)) for p in vc_root.rglob("*") if p.is_file())[:20]
            raise SkillError(
                stage="outputs",
                error_code=ErrorCode.EXPECTED_OUTPUTS_NOT_FOUND,
                message="No variant-calling outputs were detected for any requested caller.",
                fix="Inspect Nextflow logs to see whether the variant-calling step ran. Confirm the chosen callers actually produced output files.",
                details={
                    "requested_callers": sorted(requested_callers),
                    "variant_calling_dir": str(vc_root),
                    "files_found_sample": found_glob,
                },
            )

    if starts_before_or_at_annotate:
        requested_annotators = tools_lower & set(_ANNOTATORS)
        if requested_annotators:
            any_found = False
            for ann in requested_annotators:
                if annotation.get(ann):
                    any_found = True
                else:
                    missing.append(ann)
            if not any_found:
                ann_root = output_dir / "annotation"
                files = []
                if ann_root.is_dir():
                    files = sorted(
                        str(p.relative_to(output_dir))
                        for ext in _ANN_DATA_EXTS
                        for p in ann_root.rglob(f"*{ext}")
                    )[:20]
                raise SkillError(
                    stage="outputs",
                    error_code=ErrorCode.EXPECTED_OUTPUTS_NOT_FOUND,
                    message="No annotation outputs were detected for any requested annotator.",
                    fix="Check that snpeff/vep ran successfully and that the annotation step was reached.",
                    details={
                        "requested_annotators": sorted(requested_annotators),
                        "annotation_dir": str(ann_root),
                        "files_found_sample": files,
                    },
                )

    return missing


def _missing_for_report(
    *,
    step: str,
    tools_lower: set[str],
    variant_calling: dict,
    annotation: dict,
) -> list[str]:
    missing: list[str] = []
    requested_callers = tools_lower & set(_VARIANT_CALLERS)
    step_rank = {name: idx for idx, name in enumerate(_STEP_ORDER)}
    starts_before_or_at_variant_calling = step_rank.get(step, 999) <= step_rank["variant_calling"]
    starts_before_or_at_annotate = step_rank.get(step, 999) <= step_rank["annotate"]
    if starts_before_or_at_variant_calling:
        for tool in sorted(requested_callers):
            if not variant_calling.get(tool) and not variant_calling.get(f"joint_{tool}"):
                missing.append(tool)
    if starts_before_or_at_annotate:
        for ann in sorted(tools_lower & set(_ANNOTATORS)):
            if not annotation.get(ann):
                missing.append(ann)
    return missing


# ---------------------------------------------------------------------------
# step_completed inference
# ---------------------------------------------------------------------------


_STEP_ORDER = ("mapping", "markduplicates", "prepare_recalibration", "recalibrate", "variant_calling", "annotate")


def _infer_step_completed(
    *,
    requested_step: str,
    csv_handoff: dict[str, Path],
    variant_calling: dict,
    annotation: dict,
) -> str:
    achieved = requested_step
    if any(annotation.get(a) for a in _ANNOTATORS):
        achieved = "annotate"
    elif variant_calling:
        achieved = "variant_calling"
    elif "recalibrated" in csv_handoff:
        achieved = "recalibrate"
    elif "markduplicates" in csv_handoff or "markduplicates_no_table" in csv_handoff:
        achieved = "markduplicates"
    elif "mapped" in csv_handoff:
        achieved = "mapping"
    return achieved


# ---------------------------------------------------------------------------
# Public entry
# ---------------------------------------------------------------------------


def parse_outputs(
    output_dir: Path,
    *,
    step: str,
    tools: list[str],
    skip_tools: list[str] | None = None,
    analysis_mode: str,
    aligner: str,
    save_mapped: bool = False,
    save_output_as_bam: bool = False,
    save_reference: bool = False,
    validate: bool = True,
) -> OutputsReport:
    output_dir = Path(output_dir)
    warnings: list[str] = []

    tools_lower: set[str] = {t.lower() for t in (tools or [])}

    # Sarek publishes mapping results either when users explicitly request
    # them, or when duplicate marking is skipped so that mapping is itself the
    # deliverable. The skipped-MarkDuplicates route publishes BAM when
    # save_output_as_bam is enabled (aligner.config) and CRAM otherwise
    # (markduplicates.config:BAM_TO_CRAM_MAPPING). Sentieon dedup takes over
    # that latter output path when enabled.
    skip_tools_lower = {str(t).lower() for t in (skip_tools or [])}
    mapped_published = save_mapped or (
        "markduplicates" in skip_tools_lower
        and "sentieon_dedup" not in tools_lower
    )

    preprocessing, csv_handoff, csv_samples = _parse_preprocessing(
        output_dir,
        save_mapped=mapped_published,
        save_output_as_bam=save_output_as_bam,
        warnings=warnings,
    )

    variant_calling, pairs = _parse_variant_calling(
        output_dir,
        tools_lower=tools_lower,
        warnings=warnings,
    )

    annotation = _parse_annotation(output_dir, warnings)
    qc = _parse_qc(output_dir, warnings)

    pipeline_info_path = output_dir / "pipeline_info"
    pipeline_info = pipeline_info_path if pipeline_info_path.is_dir() else None

    reference_outputs = _parse_reference(output_dir, warnings) if save_reference else None

    samples_detected = _detect_samples(csv_samples, preprocessing, variant_calling)

    if validate:
        # Raises on hard failure; otherwise populates missing list below.
        _validate_outputs(
            step=step,
            tools_lower=tools_lower,
            variant_calling=variant_calling,
            annotation=annotation,
            output_dir=output_dir,
        )

    missing_outputs = _missing_for_report(
        step=step,
        tools_lower=tools_lower,
        variant_calling=variant_calling,
        annotation=annotation,
    )

    step_completed = _infer_step_completed(
        requested_step=step,
        csv_handoff=csv_handoff,
        variant_calling=variant_calling,
        annotation=annotation,
    )

    handoff_available = bool(csv_handoff)

    return OutputsReport(
        step_completed=step_completed,
        preprocessing=preprocessing,
        variant_calling=variant_calling,
        annotation=annotation,
        qc=qc,
        pipeline_info=pipeline_info,
        reference_outputs=reference_outputs,
        samples_detected=samples_detected,
        pairs_detected=pairs,
        csv_handoff=dict(csv_handoff),
        handoff_available=handoff_available,
        warnings=warnings,
        missing_outputs=missing_outputs,
    )


__all__ = ["OutputsReport", "parse_outputs"]
