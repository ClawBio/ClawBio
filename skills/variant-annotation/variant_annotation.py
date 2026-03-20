#!/usr/bin/env python3
"""variant-annotation: annotate VCF variants and summarize prioritized findings."""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib
import json
import shlex
import sys
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import requests

try:
    pysam = importlib.import_module("pysam")
except ImportError:  # pragma: no cover - handled at runtime in main()
    pysam = None

_PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from clawbio.common.checksums import sha256_file
from clawbio.common.report import (
    DISCLAIMER,
    generate_report_footer,
    generate_report_header,
    write_result_json,
)


SKILL_NAME = "variant-annotation"
SKILL_VERSION = "0.1.0"
SKILL_DIR = Path(__file__).resolve().parent
DEFAULT_DEMO_INPUT = SKILL_DIR / "example_data" / "synthetic_clinvar_panel.vcf"

VEP_BASE_URL = "https://rest.ensembl.org"
VEP_ENDPOINT = "vep/homo_sapiens/region"
DEFAULT_ASSEMBLY = "GRCh38"
DEFAULT_BATCH_SIZE = 200
DEFAULT_TIMEOUT = 30
RATE_INTERVAL = 1.0 / 15.0
DEFAULT_CACHE_TTL = 86400
DEFAULT_CACHE_DIR = Path.home() / ".clawbio" / "variant_annotation_cache"
USER_AGENT = f"ClawBio-VariantAnnotation/{SKILL_VERSION}"

CLINVAR_PATHOGENIC_TERMS = {
    "pathogenic",
    "likely_pathogenic",
    "pathogenic/likely_pathogenic",
    "likely_pathogenic/pathogenic",
}

SEVERITY_RANKING = [
    "transcript_ablation",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
    "stop_lost",
    "start_lost",
    "transcript_amplification",
    "feature_elongation",
    "feature_truncation",
    "inframe_insertion",
    "inframe_deletion",
    "missense_variant",
    "protein_altering_variant",
    "splice_donor_5th_base_variant",
    "splice_region_variant",
    "splice_donor_region_variant",
    "splice_polypyrimidine_tract_variant",
    "incomplete_terminal_codon_variant",
    "start_retained_variant",
    "stop_retained_variant",
    "synonymous_variant",
    "coding_sequence_variant",
    "mature_miRNA_variant",
    "5_prime_UTR_variant",
    "3_prime_UTR_variant",
    "non_coding_transcript_exon_variant",
    "intron_variant",
    "NMD_transcript_variant",
    "non_coding_transcript_variant",
    "coding_transcript_variant",
    "upstream_gene_variant",
    "downstream_gene_variant",
    "TFBS_ablation",
    "TFBS_amplification",
    "TF_binding_site_variant",
    "regulatory_region_ablation",
    "regulatory_region_amplification",
    "regulatory_region_variant",
    "intergenic_variant",
]
SEVERITY_INDEX = {term: idx for idx, term in enumerate(SEVERITY_RANKING)}

ANNOTATED_COLUMNS = [
    "input_id",
    "chrom",
    "pos",
    "ref",
    "alt",
    "genotype",
    "gene",
    "gene_id",
    "feature_type",
    "feature_id",
    "consequence",
    "all_consequences",
    "impact",
    "impact_tier",
    "consequence_severity_rank",
    "hgvsc",
    "hgvsp",
    "clinvar_significance",
    "clinvar_bucket",
    "clinvar_accessions",
    "is_rare",
    "gnomad_af",
    "global_af",
    "max_af",
    "min_af",
    "highest_frequency_population",
    "population_frequency_spread",
    "population_freq_summary",
    "population_outlier_flag",
    "variant_class",
    "existing_variation",
    "priority_score",
    "priority_bucket",
    "review_reasons",
    "clinically_relevant",
    "vep_status",
    "warning",
]

IMPACT_SCORES = {
    "HIGH": 4,
    "MODERATE": 3,
    "LOW": 2,
    "MODIFIER": 1,
}

CLINVAR_BUCKET_SCORES = {
    "pathogenic_or_likely_pathogenic": 4,
    "conflicting": 2,
    "risk_factor": 2,
    "drug_response": 2,
    "other_or_uncertain": 1,
    "benign_or_likely_benign": 0,
}


class VEPClient:
    """Rate-limited, caching client for Ensembl VEP POST requests."""

    def __init__(
        self,
        cache_dir: Path = DEFAULT_CACHE_DIR,
        use_cache: bool = True,
        timeout: int = DEFAULT_TIMEOUT,
        cache_ttl: int = DEFAULT_CACHE_TTL,
        rate_interval: float = RATE_INTERVAL,
    ):
        self.cache_dir = cache_dir
        self.use_cache = use_cache
        self.timeout = timeout
        self.cache_ttl = cache_ttl
        self.rate_interval = rate_interval
        self._last_request_time = 0.0

        self.session = requests.Session()
        self.session.headers.update(
            {
                "Accept": "application/json",
                "Content-Type": "application/json",
                "User-Agent": USER_AGENT,
            }
        )

        if self.use_cache:
            self.cache_dir.mkdir(parents=True, exist_ok=True)

    def _throttle(self) -> None:
        elapsed = time.time() - self._last_request_time
        if elapsed < self.rate_interval:
            time.sleep(self.rate_interval - elapsed)
        self._last_request_time = time.time()

    def _cache_key(
        self, endpoint: str, params: dict[str, Any], payload: dict[str, Any]
    ) -> str:
        raw = f"POST|{endpoint}|{json.dumps(params, sort_keys=True)}|{json.dumps(payload, sort_keys=True)}"
        return hashlib.sha256(raw.encode("utf-8")).hexdigest()[:20]

    def _get_cached(self, key: str) -> Any | None:
        if not self.use_cache:
            return None
        path = self.cache_dir / f"{key}.json"
        if not path.exists():
            return None
        try:
            payload = json.loads(path.read_text())
        except json.JSONDecodeError:
            return None
        cached_at = payload.get("_cached_at", 0)
        if time.time() - cached_at >= self.cache_ttl:
            return None
        return payload.get("response")

    def _set_cached(self, key: str, response_data: Any) -> None:
        if not self.use_cache:
            return
        path = self.cache_dir / f"{key}.json"
        path.write_text(
            json.dumps(
                {"_cached_at": time.time(), "response": response_data},
                indent=2,
                default=str,
            )
        )

    def annotate_batch(
        self, variant_strings: list[str], assembly: str = DEFAULT_ASSEMBLY
    ) -> list[dict[str, Any]]:
        params: dict[str, Any] = {
            "canonical": 1,
            "variant_class": 1,
            "hgvs": 1,
            "vcf_string": 1,
            "af_gnomad": 1,
            "af": 1,
            "mane": 1,
        }
        if assembly:
            params["assembly"] = assembly

        payload = {"variants": variant_strings}
        cache_key = self._cache_key(VEP_ENDPOINT, params, payload)
        cached = self._get_cached(cache_key)
        if cached is not None:
            return cached

        self._throttle()
        url = f"{VEP_BASE_URL}/{VEP_ENDPOINT}"
        response = self.session.post(
            url, params=params, json=payload, timeout=self.timeout
        )
        if response.status_code == 429:
            retry_after = float(response.headers.get("Retry-After", "2"))
            time.sleep(retry_after)
            response = self.session.post(
                url, params=params, json=payload, timeout=self.timeout
            )

        response.raise_for_status()
        data = response.json()
        self._set_cached(cache_key, data)
        return data


def parse_vcf(input_path: Path) -> list[dict[str, Any]]:
    """Parse a single-sample or multi-sample VCF via pysam."""
    records: list[dict[str, Any]] = []

    if pysam is None:
        raise ImportError(
            "pysam is required for VCF parsing. Install it with `pip install pysam`."
        )

    with pysam.VariantFile(str(input_path)) as vcf:
        sample_names = list(vcf.header.samples)
        sample_name = sample_names[0] if sample_names else None

        for record in vcf:
            genotype = extract_genotype(record, sample_name)
            variant_id = (
                record.id
                or f"{record.chrom}:{record.pos}:{record.ref}:{','.join(record.alts or [])}"
            )

            for alt in record.alts or []:
                records.append(
                    {
                        "input_id": variant_id,
                        "chrom": str(record.chrom),
                        "pos": int(record.pos),
                        "ref": record.ref,
                        "alt": alt,
                        "qual": record.qual,
                        "filter": list(record.filter.keys()),
                        "info": {
                            key: normalise_info_value(value)
                            for key, value in record.info.items()
                        },
                        "sample_name": sample_name,
                        "genotype": genotype,
                    }
                )

    return records


def extract_genotype(record: Any, sample_name: str | None) -> str:
    if not sample_name:
        return "."

    sample = record.samples.get(sample_name)
    if sample is None:
        return "."

    alleles = sample.get("GT")
    if not alleles:
        return "."

    rendered: list[str] = []
    all_alleles = (record.ref,) + tuple(record.alts or ())
    for allele_index in alleles:
        if allele_index is None:
            rendered.append(".")
        elif 0 <= allele_index < len(all_alleles):
            rendered.append(str(all_alleles[allele_index]))
        else:
            rendered.append(".")

    phased = bool(sample.phased)
    delimiter = "|" if phased else "/"
    return delimiter.join(rendered)


def normalise_info_value(value: Any) -> Any:
    if isinstance(value, tuple):
        return [normalise_info_value(item) for item in value]
    return value


def batch_variants(
    records: list[dict[str, Any]], batch_size: int = DEFAULT_BATCH_SIZE
) -> list[list[dict[str, Any]]]:
    return [records[i : i + batch_size] for i in range(0, len(records), batch_size)]


def build_vep_variant_string(record: dict[str, Any]) -> str:
    chrom = str(record["chrom"]).removeprefix("chr")
    pos = int(record["pos"])
    return f"{chrom} {pos} {pos} {record['ref']}/{record['alt']} +"


def consequence_rank(term: str) -> int:
    return SEVERITY_INDEX.get(term, len(SEVERITY_INDEX) + 100)


def choose_best_annotation(entry: dict[str, Any]) -> dict[str, Any]:
    candidates: list[dict[str, Any]] = []

    for feature_type, consequences in (
        ("transcript", entry.get("transcript_consequences", [])),
        ("regulatory", entry.get("regulatory_feature_consequences", [])),
        ("motif", entry.get("motif_feature_consequences", [])),
        ("intergenic", entry.get("intergenic_consequences", [])),
    ):
        for consequence in consequences:
            terms = consequence.get("consequence_terms", []) or [
                entry.get("most_severe_consequence", "")
            ]
            best_term = min(terms, key=consequence_rank)
            candidates.append(
                {
                    "feature_type": feature_type,
                    "feature_id": consequence.get("transcript_id")
                    or consequence.get("regulatory_feature_id")
                    or consequence.get("motif_feature_id")
                    or "",
                    "gene": consequence.get("gene_symbol")
                    or consequence.get("gene_id")
                    or "",
                    "gene_id": consequence.get("gene_id", ""),
                    "consequence": best_term,
                    "impact": consequence.get("impact", entry.get("impact", "")),
                    "hgvsc": consequence.get("hgvsc", ""),
                    "hgvsp": consequence.get("hgvsp", ""),
                    "canonical": int(bool(consequence.get("canonical"))),
                    "mane": int(
                        bool(
                            consequence.get("mane_select")
                            or consequence.get("mane_plus_clinical")
                        )
                    ),
                }
            )

    if candidates:
        candidates.sort(
            key=lambda item: (
                consequence_rank(item["consequence"]),
                -item["canonical"],
                -item["mane"],
                item["gene"] == "",
            )
        )
        return candidates[0]

    return {
        "feature_type": "",
        "feature_id": "",
        "gene": "",
        "gene_id": "",
        "consequence": entry.get("most_severe_consequence", ""),
        "impact": entry.get("impact", ""),
        "hgvsc": "",
        "hgvsp": "",
        "canonical": 0,
        "mane": 0,
    }


def collect_consequences(entry: dict[str, Any]) -> list[str]:
    consequences: list[str] = []
    for key in (
        "transcript_consequences",
        "regulatory_feature_consequences",
        "motif_feature_consequences",
        "intergenic_consequences",
    ):
        for consequence in entry.get(key, []):
            for term in ensure_list(consequence.get("consequence_terms")):
                if term and term not in consequences:
                    consequences.append(str(term))
    if not consequences and entry.get("most_severe_consequence"):
        consequences.append(str(entry["most_severe_consequence"]))
    return consequences


def extract_clinvar_fields(entry: dict[str, Any]) -> dict[str, str]:
    significances: list[str] = []
    accessions: list[str] = []

    for colocated in entry.get("colocated_variants", []):
        for value in ensure_list(colocated.get("clin_sig")):
            normalised = str(value).strip().replace(" ", "_")
            if normalised:
                significances.append(normalised)
        for key in ("clinvar_accession", "clinvar_id", "id"):
            raw = colocated.get(key)
            if raw is None:
                continue
            if key == "id" and not str(raw).startswith("rs"):
                continue
            accessions.append(str(raw))

    return {
        "clinvar_significance": join_unique(significances),
        "clinvar_accessions": join_unique(accessions),
    }


def extract_frequency_fields(entry: dict[str, Any], alt: str) -> dict[str, Any]:
    alt_key = alt.upper()
    gnomad_values: list[float] = []
    all_values: list[float] = []
    population_freqs: dict[str, float] = {}

    for colocated in entry.get("colocated_variants", []):
        freqs = colocated.get("frequencies") or {}
        if alt_key in freqs and isinstance(freqs[alt_key], dict):
            for label, value in freqs[alt_key].items():
                numeric = to_float(value)
                if numeric is None:
                    continue
                all_values.append(numeric)
                population_freqs[label] = numeric
                if "gnomad" in label.lower():
                    gnomad_values.append(numeric)

        for label, value in colocated.items():
            if not isinstance(value, (int, float, str)):
                continue
            numeric = to_float(value)
            if numeric is None:
                continue
            label_lower = label.lower()
            if "af" not in label_lower:
                continue
            all_values.append(numeric)
            population_freqs[label] = numeric
            if "gnomad" in label_lower:
                gnomad_values.append(numeric)

    highest_population = ""
    min_af = None
    spread = None
    if population_freqs:
        highest_population = max(population_freqs.items(), key=lambda item: item[1])[0]
        min_af = min(population_freqs.values())
        spread = max(population_freqs.values()) - min_af

    return {
        "gnomad_af": min(gnomad_values) if gnomad_values else None,
        "global_af": population_freqs.get("gnomad")
        or population_freqs.get("af")
        or population_freqs.get("gnomad_af"),
        "max_af": max(all_values) if all_values else None,
        "min_af": min_af,
        "highest_frequency_population": highest_population,
        "population_frequency_spread": spread,
        "population_freq_summary": summarize_population_frequencies(population_freqs),
        "population_outlier_flag": spread is not None and spread >= 0.05,
    }


def ensure_list(value: Any) -> list[Any]:
    if value is None:
        return []
    if isinstance(value, list):
        return value
    if isinstance(value, tuple):
        return list(value)
    return [value]


def join_unique(values: list[str]) -> str:
    unique = []
    for value in values:
        if value and value not in unique:
            unique.append(value)
    return ";".join(unique)


def to_float(value: Any) -> float | None:
    try:
        if value in (None, "", "."):
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def has_pathogenic_clinvar(significance: str) -> bool:
    if not significance:
        return False
    parts = {
        part.strip().lower()
        for part in significance.replace(";", "/").split("/")
        if part.strip()
    }
    return bool(parts & CLINVAR_PATHOGENIC_TERMS)


def classify_clinvar_bucket(significance: str) -> str:
    if not significance:
        return "other_or_uncertain"

    lowered = significance.lower()
    if "conflict" in lowered:
        return "conflicting"
    if has_pathogenic_clinvar(significance):
        return "pathogenic_or_likely_pathogenic"
    if "benign" in lowered:
        return "benign_or_likely_benign"
    if "risk_factor" in lowered or "risk factor" in lowered:
        return "risk_factor"
    if "drug_response" in lowered or "drug response" in lowered:
        return "drug_response"
    return "other_or_uncertain"


def impact_tier(impact: str) -> str:
    impact = (impact or "").upper()
    return impact if impact in IMPACT_SCORES else "UNKNOWN"


def is_rare_variant(gnomad_af: float | None) -> bool:
    return gnomad_af is not None and gnomad_af < 0.001


def summarize_population_frequencies(population_freqs: dict[str, float]) -> str:
    if not population_freqs:
        return ""
    ordered = sorted(population_freqs.items(), key=lambda item: (-item[1], item[0]))
    return ";".join(f"{label}={value:.6g}" for label, value in ordered)


def review_reasons_for(annotation: dict[str, Any]) -> list[str]:
    reasons: list[str] = []
    if annotation.get("clinically_relevant"):
        reasons.append("rare_pathogenic_clinvar")
    if annotation.get("impact_tier") == "HIGH":
        reasons.append("high_impact")
    if annotation.get("impact_tier") == "MODERATE":
        reasons.append("moderate_impact")
    if annotation.get("population_outlier_flag"):
        reasons.append("population_frequency_outlier")
    bucket = annotation.get("clinvar_bucket")
    if bucket in {"risk_factor", "drug_response", "conflicting"}:
        reasons.append(bucket)
    return reasons


def compute_priority(annotation: dict[str, Any]) -> tuple[int, str, list[str]]:
    score = 0
    score += IMPACT_SCORES.get(annotation.get("impact_tier", ""), 0) * 10
    score += CLINVAR_BUCKET_SCORES.get(annotation.get("clinvar_bucket", ""), 0) * 15
    if annotation.get("is_rare"):
        score += 20
    if annotation.get("population_outlier_flag"):
        score += 10
    severity_rank = annotation.get("consequence_severity_rank")
    if isinstance(severity_rank, int):
        score += max(0, 40 - min(severity_rank, 40))

    reasons = review_reasons_for(annotation)

    if annotation.get("clinically_relevant"):
        bucket = "Tier 1"
    elif annotation.get("impact_tier") in {"HIGH", "MODERATE"} and (
        annotation.get("is_rare")
        or annotation.get("clinvar_bucket") in {"conflicting", "other_or_uncertain"}
    ):
        bucket = "Tier 2"
    elif annotation.get("clinvar_bucket") in {
        "risk_factor",
        "drug_response",
    } or annotation.get("impact_tier") in {"LOW", "MODERATE"}:
        bucket = "Tier 3"
    else:
        bucket = "Tier 4"

    return score, bucket, reasons


def is_clinically_relevant(annotation: dict[str, Any]) -> bool:
    gnomad_af = annotation.get("gnomad_af")
    if gnomad_af is None:
        return False
    return gnomad_af < 0.001 and has_pathogenic_clinvar(
        annotation.get("clinvar_significance", "")
    )


def normalize_annotation(
    record: dict[str, Any], entry: dict[str, Any]
) -> dict[str, Any]:
    best = choose_best_annotation(entry)
    clinvar = extract_clinvar_fields(entry)
    freqs = extract_frequency_fields(entry, record["alt"])
    all_consequences = collect_consequences(entry)
    clinvar_bucket = classify_clinvar_bucket(clinvar["clinvar_significance"])
    picked_impact_tier = impact_tier(best["impact"])
    picked_severity_rank = consequence_rank(best["consequence"])
    existing_variation = join_unique(
        [
            str(v)
            for v in ensure_list(entry.get("id")) + ensure_list(entry.get("input"))
            if v
        ]
    )

    annotation = {
        "input_id": record["input_id"],
        "chrom": record["chrom"],
        "pos": record["pos"],
        "ref": record["ref"],
        "alt": record["alt"],
        "genotype": record.get("genotype", "."),
        "gene": best["gene"],
        "gene_id": best["gene_id"],
        "feature_type": best["feature_type"],
        "feature_id": best["feature_id"],
        "consequence": best["consequence"],
        "all_consequences": join_unique(all_consequences),
        "impact": best["impact"],
        "impact_tier": picked_impact_tier,
        "consequence_severity_rank": picked_severity_rank,
        "hgvsc": best["hgvsc"],
        "hgvsp": best["hgvsp"],
        "clinvar_significance": clinvar["clinvar_significance"],
        "clinvar_bucket": clinvar_bucket,
        "clinvar_accessions": clinvar["clinvar_accessions"],
        "is_rare": False,
        "gnomad_af": freqs["gnomad_af"],
        "global_af": freqs["global_af"],
        "max_af": freqs["max_af"],
        "min_af": freqs["min_af"],
        "highest_frequency_population": freqs["highest_frequency_population"],
        "population_frequency_spread": freqs["population_frequency_spread"],
        "population_freq_summary": freqs["population_freq_summary"],
        "population_outlier_flag": freqs["population_outlier_flag"],
        "variant_class": entry.get("variant_class", ""),
        "existing_variation": existing_variation,
        "priority_score": 0,
        "priority_bucket": "Tier 4",
        "review_reasons": "",
        "clinically_relevant": False,
        "vep_status": "ok",
        "warning": join_unique([str(w) for w in ensure_list(entry.get("warnings"))]),
    }
    annotation["is_rare"] = is_rare_variant(annotation["gnomad_af"])
    annotation["clinically_relevant"] = is_clinically_relevant(annotation)
    priority_score, priority_bucket, reasons = compute_priority(annotation)
    annotation["priority_score"] = priority_score
    annotation["priority_bucket"] = priority_bucket
    annotation["review_reasons"] = join_unique(reasons)
    return annotation


def annotate_variants(
    records: list[dict[str, Any]],
    client: VEPClient,
    assembly: str = DEFAULT_ASSEMBLY,
    batch_size: int = DEFAULT_BATCH_SIZE,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    annotations: list[dict[str, Any]] = []
    failures: list[dict[str, Any]] = []
    batch_count = 0
    failed_batches = 0

    for batch in batch_variants(records, batch_size=batch_size):
        batch_count += 1
        variant_strings = [build_vep_variant_string(record) for record in batch]
        try:
            response = client.annotate_batch(variant_strings, assembly=assembly)
        except Exception as exc:
            failed_batches += 1
            message = str(exc)
            for record in batch:
                failure = {
                    "input_id": record["input_id"],
                    "chrom": record["chrom"],
                    "pos": record["pos"],
                    "ref": record["ref"],
                    "alt": record["alt"],
                    "error": message,
                }
                failures.append(failure)
                annotations.append(placeholder_annotation(record, "error", message))
            continue

        for index, record in enumerate(batch):
            if index >= len(response):
                message = "VEP response did not include an annotation for this variant"
                failures.append(
                    {
                        "input_id": record["input_id"],
                        "chrom": record["chrom"],
                        "pos": record["pos"],
                        "ref": record["ref"],
                        "alt": record["alt"],
                        "error": message,
                    }
                )
                annotations.append(placeholder_annotation(record, "missing", message))
                continue
            annotations.append(normalize_annotation(record, response[index]))

    metadata = {
        "assembly": assembly,
        "batch_size": batch_size,
        "batches_sent": batch_count,
        "failed_batches": failed_batches,
    }
    return annotations, failures, metadata


def placeholder_annotation(
    record: dict[str, Any], status: str, warning: str
) -> dict[str, Any]:
    return {
        "input_id": record["input_id"],
        "chrom": record["chrom"],
        "pos": record["pos"],
        "ref": record["ref"],
        "alt": record["alt"],
        "genotype": record.get("genotype", "."),
        "gene": "",
        "gene_id": "",
        "feature_type": "",
        "feature_id": "",
        "consequence": "",
        "all_consequences": "",
        "impact": "",
        "impact_tier": "UNKNOWN",
        "consequence_severity_rank": "",
        "hgvsc": "",
        "hgvsp": "",
        "clinvar_significance": "",
        "clinvar_bucket": "other_or_uncertain",
        "clinvar_accessions": "",
        "is_rare": False,
        "gnomad_af": "",
        "global_af": "",
        "max_af": "",
        "min_af": "",
        "highest_frequency_population": "",
        "population_frequency_spread": "",
        "population_freq_summary": "",
        "population_outlier_flag": False,
        "variant_class": "",
        "existing_variation": "",
        "priority_score": 0,
        "priority_bucket": "Tier 4",
        "review_reasons": "",
        "clinically_relevant": False,
        "vep_status": status,
        "warning": warning,
    }


def write_tsv(output_dir: Path, annotations: list[dict[str, Any]]) -> Path:
    tables_dir = output_dir / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)
    output_path = tables_dir / "annotated_variants.tsv"

    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=ANNOTATED_COLUMNS, delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        for row in annotations:
            writer.writerow(
                {key: render_scalar(row.get(key)) for key in ANNOTATED_COLUMNS}
            )

    return output_path


def render_scalar(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    return str(value)


def summarize_annotations(annotations: list[dict[str, Any]]) -> dict[str, Any]:
    clinically_relevant = [row for row in annotations if row.get("clinically_relevant")]
    annotated_ok = [row for row in annotations if row.get("vep_status") == "ok"]
    annotated_missing = [
        row for row in annotations if row.get("vep_status") == "missing"
    ]
    annotated_error = [row for row in annotations if row.get("vep_status") == "error"]

    def safe_priority_score(row: dict[str, Any]) -> float:
        score = to_float(row.get("priority_score"))
        return score if score is not None else 0.0

    top_variants = sorted(
        annotations,
        key=lambda row: (
            -safe_priority_score(row),
            consequence_rank(str(row.get("consequence", ""))),
            str(row.get("input_id", "")),
        ),
    )
    pathogenic = [
        row
        for row in annotations
        if row.get("clinvar_bucket") == "pathogenic_or_likely_pathogenic"
    ]
    genes = Counter(row["gene"] for row in annotations if row.get("gene"))
    consequences = Counter(
        row["consequence"] for row in annotations if row.get("consequence")
    )
    impact_counts = Counter(
        row["impact_tier"] for row in annotations if row.get("impact_tier")
    )
    clinvar_counts = Counter(
        row["clinvar_bucket"] for row in annotations if row.get("clinvar_bucket")
    )
    tier_counts = Counter(
        row["priority_bucket"] for row in annotations if row.get("priority_bucket")
    )
    population_outliers = [
        row for row in annotations if row.get("population_outlier_flag")
    ]
    notable_non_flagged = [
        row
        for row in top_variants
        if not row.get("clinically_relevant")
        and row.get("vep_status") == "ok"
        and row.get("priority_bucket") in {"Tier 1", "Tier 2", "Tier 3"}
    ]
    warnings = [row for row in annotations if row.get("warning")]

    return {
        "clinically_relevant": clinically_relevant,
        "top_variants": top_variants[:10],
        "pathogenic": pathogenic,
        "genes": genes,
        "consequences": consequences,
        "impact_counts": impact_counts,
        "clinvar_counts": clinvar_counts,
        "tier_counts": tier_counts,
        "population_outliers": population_outliers,
        "annotated_ok": annotated_ok,
        "annotated_missing": annotated_missing,
        "annotated_error": annotated_error,
        "notable_non_flagged": notable_non_flagged[:5],
        "warnings": warnings,
    }


def format_af(value: Any) -> str:
    numeric = to_float(value)
    if numeric is None:
        return "NA"
    return f"{numeric:.6g}"


def format_reason_label(reason: str) -> str:
    mapping = {
        "rare_pathogenic_clinvar": "rare pathogenic ClinVar hit",
        "high_impact": "high predicted impact",
        "moderate_impact": "moderate predicted impact",
        "population_frequency_outlier": "population frequency outlier",
        "risk_factor": "ClinVar risk factor",
        "drug_response": "ClinVar drug response",
        "conflicting": "conflicting ClinVar evidence",
    }
    return mapping.get(reason, reason.replace("_", " "))


def format_review_reasons(value: Any) -> str:
    reasons = [
        format_reason_label(reason) for reason in str(value or "").split(";") if reason
    ]
    return ", ".join(reasons) if reasons else "No specific prioritization note"


def summarize_population_context(row: dict[str, Any]) -> str:
    highest_population = row.get("highest_frequency_population") or "NA"
    highest_af = format_af(row.get("max_af"))
    min_af = format_af(row.get("min_af"))
    spread = format_af(row.get("population_frequency_spread"))
    return (
        f"highest {highest_population}={highest_af}; min AF {min_af}; spread {spread}"
    )


def explain_non_flagged_reason(row: dict[str, Any]) -> str:
    if row.get("clinically_relevant"):
        return "Flagged as clinically relevant"
    if row.get("vep_status") != "ok":
        return f"Annotation status: {row.get('vep_status') or 'unknown'}"
    if not has_pathogenic_clinvar(str(row.get("clinvar_significance", ""))):
        return "Not flagged because ClinVar does not include a pathogenic label"
    gnomad_af = to_float(row.get("gnomad_af"))
    if gnomad_af is None:
        return "Not flagged because gnomAD AF was unavailable for the clinical relevance rule"
    if gnomad_af >= 0.001:
        return f"Not flagged because gnomAD AF {gnomad_af:.6g} is above the 0.001 rarity cutoff"
    return "Not flagged because the clinical relevance rule was not fully met"


def plot_summary_figure(summary: dict[str, Any], output_dir: Path) -> Path | None:
    """Generate a four-panel summary figure when matplotlib is available."""
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return None

    figures_dir = output_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    figure_path = figures_dir / "variant_annotation_summary.png"

    tier_counts = summary["tier_counts"]
    impact_counts = summary["impact_counts"]
    clinvar_counts = summary["clinvar_counts"]
    top_variants = summary["top_variants"][:5]
    clinically_relevant = summary["clinically_relevant"]
    annotated_ok = summary["annotated_ok"]

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    fig.patch.set_facecolor("white")

    tier_order = ["Tier 1", "Tier 2", "Tier 3", "Tier 4"]
    tier_values = [tier_counts.get(label, 0) for label in tier_order]
    tier_colors = ["#b91c1c", "#ea580c", "#2563eb", "#94a3b8"]
    axes[0, 0].bar(tier_order, tier_values, color=tier_colors, edgecolor="white")
    axes[0, 0].set_title("Priority Buckets")
    axes[0, 0].set_ylabel("Variant count")
    axes[0, 0].grid(axis="y", alpha=0.25)

    impact_order = ["HIGH", "MODERATE", "LOW", "MODIFIER", "UNKNOWN"]
    impact_values = [impact_counts.get(label, 0) for label in impact_order]
    impact_colors = ["#991b1b", "#f97316", "#facc15", "#60a5fa", "#cbd5e1"]
    axes[0, 1].bar(impact_order, impact_values, color=impact_colors, edgecolor="white")
    axes[0, 1].set_title("Impact Tiers")
    axes[0, 1].set_ylabel("Variant count")
    axes[0, 1].tick_params(axis="x", rotation=20)
    axes[0, 1].grid(axis="y", alpha=0.25)

    clinvar_order = [
        "pathogenic_or_likely_pathogenic",
        "conflicting",
        "risk_factor",
        "drug_response",
        "other_or_uncertain",
        "benign_or_likely_benign",
    ]
    clinvar_values = [clinvar_counts.get(label, 0) for label in clinvar_order]
    clinvar_labels = [
        "Path/LPath",
        "Conflicting",
        "Risk",
        "Drug",
        "Uncertain",
        "Benign",
    ]
    clinvar_colors = ["#b91c1c", "#7c3aed", "#c2410c", "#0f766e", "#64748b", "#2563eb"]
    axes[1, 0].barh(
        clinvar_labels, clinvar_values, color=clinvar_colors, edgecolor="white"
    )
    axes[1, 0].set_title("ClinVar Summary")
    axes[1, 0].set_xlabel("Variant count")
    axes[1, 0].invert_yaxis()
    axes[1, 0].grid(axis="x", alpha=0.25)

    top_scores = [to_float(row.get("priority_score")) or 0.0 for row in top_variants]
    top_labels = [
        f"{row['input_id']} | {row.get('gene') or 'NA'}" for row in top_variants
    ]
    bucket_colors = {
        "Tier 1": "#b91c1c",
        "Tier 2": "#ea580c",
        "Tier 3": "#2563eb",
        "Tier 4": "#94a3b8",
    }
    top_colors = [
        bucket_colors.get(row.get("priority_bucket"), "#94a3b8") for row in top_variants
    ]
    axes[1, 1].barh(top_labels, top_scores, color=top_colors, edgecolor="white")
    axes[1, 1].set_title("Top Prioritized Variants")
    axes[1, 1].set_xlabel("Priority score")
    axes[1, 1].invert_yaxis()
    axes[1, 1].grid(axis="x", alpha=0.25)
    for index, row in enumerate(top_variants):
        score = top_scores[index]
        axes[1, 1].text(
            score + max(1.0, score * 0.01),
            index,
            row.get("priority_bucket") or "NA",
            va="center",
            fontsize=8,
        )

    total_variants = (
        len(summary["annotated_ok"])
        + len(summary["annotated_missing"])
        + len(summary["annotated_error"])
    )
    fig.suptitle("Variant Annotation Summary", fontsize=16, fontweight="bold")
    fig.text(
        0.99,
        0.98,
        "\n".join(
            [
                f"Parsed: {total_variants}",
                f"Annotated: {len(annotated_ok)}",
                f"Clinically relevant: {len(clinically_relevant)}",
                f"Top score: {int(max(top_scores)) if top_scores else 0}",
            ]
        ),
        ha="right",
        va="top",
        fontsize=9,
        bbox={
            "boxstyle": "round,pad=0.4",
            "facecolor": "#f8fafc",
            "edgecolor": "#cbd5e1",
        },
    )
    fig.tight_layout(rect=(0, 0, 0.98, 0.95))
    fig.savefig(figure_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return figure_path


def generate_markdown_report(
    input_path: Path,
    output_dir: Path,
    annotations: list[dict[str, Any]],
    failures: list[dict[str, Any]],
    metadata: dict[str, Any],
    tsv_path: Path,
    figure_path: Path | None = None,
) -> str:
    summary = summarize_annotations(annotations)
    clinically_relevant = summary["clinically_relevant"]
    pathogenic = summary["pathogenic"]
    genes = summary["genes"]
    consequences = summary["consequences"]
    impact_counts = summary["impact_counts"]
    clinvar_counts = summary["clinvar_counts"]
    tier_counts = summary["tier_counts"]
    population_outliers = summary["population_outliers"]
    top_variants = summary["top_variants"]
    annotated_ok = summary["annotated_ok"]
    annotated_missing = summary["annotated_missing"]
    annotated_error = summary["annotated_error"]
    notable_non_flagged = summary["notable_non_flagged"]
    warnings = summary["warnings"]

    header = generate_report_header(
        title="Variant Annotation Report",
        skill_name=SKILL_NAME,
        input_files=[input_path],
        extra_metadata={
            "Assembly": metadata["assembly"],
            "Batches": str(metadata["batches_sent"]),
            "Output": str(output_dir),
        },
    )

    lines = [
        "## Overview",
        "",
        f"- Variants parsed: {len(annotations)}",
        f"- Variants annotated successfully: {len(annotated_ok)}",
        f"- Clinically relevant findings: {len(clinically_relevant)}",
        f"- Pathogenic / likely pathogenic ClinVar calls: {len(pathogenic)}",
        f"- Population-frequency outliers: {len(population_outliers)}",
        f"- Failed variant batches: {metadata['failed_batches']}",
        f"- Annotated TSV: `{tsv_path.relative_to(output_dir)}`",
        "",
        "## Executive Summary",
        "",
    ]

    if clinically_relevant:
        lines.append(
            f"- {len(clinically_relevant)} variants met the strict clinical relevance rule: a pathogenic / likely pathogenic ClinVar label plus gnomAD AF below 0.001."
        )
    else:
        lines.append(
            "- No variants met the strict clinical relevance rule of pathogenic ClinVar evidence plus gnomAD AF below 0.001."
        )

    if pathogenic:
        lines.append(
            f"- {len(pathogenic)} variants carried pathogenic or likely pathogenic ClinVar labels, but many remained outside the clinically relevant set because they were too common or lacked gnomAD frequency support."
        )
    else:
        lines.append(
            "- No variants carried pathogenic or likely pathogenic ClinVar labels."
        )

    if top_variants:
        top = top_variants[0]
        lines.append(
            f"- Top-ranked variant: `{top['input_id']}` in `{top.get('gene') or 'NA'}` ({top.get('priority_bucket') or 'NA'}, score {top.get('priority_score') or 0}) driven by {format_review_reasons(top.get('review_reasons'))}."
        )

    if population_outliers:
        lines.append(
            f"- {len(population_outliers)} variants showed notable frequency spread across populations; review these carefully before treating rarity as globally representative."
        )

    if annotated_missing or annotated_error or warnings:
        lines.append(
            f"- Data quality note: {len(annotated_missing)} variants were missing VEP output, {len(annotated_error)} had request-level errors, and {len(warnings)} carried VEP warnings."
        )
    else:
        lines.append(
            "- Data quality note: all parsed variants returned without request-level errors or VEP warnings."
        )

    lines.extend(
        [
            "",
            "## Data Quality / Coverage",
            "",
            f"- VEP status `ok`: {len(annotated_ok)}",
            f"- VEP status `missing`: {len(annotated_missing)}",
            f"- VEP status `error`: {len(annotated_error)}",
            f"- Variants with warning metadata: {len(warnings)}",
            "",
            "## Top Prioritized Variants",
            "",
        ]
    )

    if top_variants:
        lines.extend(
            [
                "| Variant | Gene | Genotype | Change | Priority | Rationale |",
                "|---|---|---|---|---|---|",
            ]
        )
        for row in top_variants[:5]:
            change = (
                row.get("hgvsp") or row.get("hgvsc") or row.get("consequence") or "NA"
            )
            lines.append(
                f"| {row['input_id']} | {row.get('gene') or 'NA'} | {row.get('genotype') or 'NA'} | {change} | "
                f"{row.get('priority_bucket') or 'NA'} ({row.get('priority_score') or 0}) | {format_review_reasons(row.get('review_reasons'))} |"
            )
    else:
        lines.append("No prioritized variants were available.")

    lines.extend(
        [
            "",
            "## Clinically Relevant Findings",
            "",
        ]
    )

    if clinically_relevant:
        lines.extend(
            [
                "| Variant | Gene | Genotype | Change | ClinVar | Evidence |",
                "|---|---|---|---|---|---|",
            ]
        )
        for row in clinically_relevant:
            variant_label = f"{row['input_id']} ({row['chrom']}:{row['pos']} {row['ref']}>{row['alt']})"
            change = (
                row.get("hgvsp") or row.get("hgvsc") or row.get("consequence") or "NA"
            )
            lines.append(
                f"| {variant_label} | {row.get('gene') or 'NA'} | {row.get('genotype') or 'NA'} | {change} | "
                f"{row.get('clinvar_significance') or 'NA'} | gnomAD AF {format_af(row.get('gnomad_af'))}; {format_review_reasons(row.get('review_reasons'))} |"
            )
    else:
        lines.append(
            "No variants met the rare (`gnomAD AF < 0.001`) and ClinVar pathogenicity filter."
        )

    lines.extend(
        [
            "",
            "## Notable Non-Flagged Variants",
            "",
        ]
    )
    if notable_non_flagged:
        lines.extend(
            [
                "| Variant | Gene | Priority | Why it was not flagged |",
                "|---|---|---|---|",
            ]
        )
        for row in notable_non_flagged:
            lines.append(
                f"| {row['input_id']} | {row.get('gene') or 'NA'} | {row.get('priority_bucket') or 'NA'} ({row.get('priority_score') or 0}) | {explain_non_flagged_reason(row)} |"
            )
    else:
        lines.append(
            "No additional high-priority non-flagged variants stood out in this run."
        )

    lines.extend(
        [
            "",
            "## Annotation Summary",
            "",
            "### Priority Buckets",
            "",
        ]
    )
    if tier_counts:
        for bucket, count in sorted(tier_counts.items()):
            lines.append(f"- {bucket}: {count}")
    else:
        lines.append("- No priority buckets were assigned.")

    lines.extend(
        [
            "",
            "### Impact Tiers",
            "",
        ]
    )
    if impact_counts:
        for tier in ("HIGH", "MODERATE", "LOW", "MODIFIER", "UNKNOWN"):
            if tier in impact_counts:
                lines.append(f"- {tier}: {impact_counts[tier]}")
    else:
        lines.append("- No impact tiers were returned.")

    lines.extend(
        [
            "",
            "### ClinVar Summary",
            "",
        ]
    )
    if clinvar_counts:
        for bucket, count in clinvar_counts.most_common():
            lines.append(f"- {bucket}: {count}")
    else:
        lines.append("- No ClinVar categories were returned.")

    lines.extend(
        [
            "",
            "### Most Frequent Consequences",
            "",
        ]
    )
    if consequences:
        for consequence, count in consequences.most_common(10):
            lines.append(f"- {consequence}: {count}")
    else:
        lines.append("- No consequence annotations were returned.")

    lines.extend(
        [
            "",
            "### Genes Hit",
            "",
        ]
    )
    if genes:
        for gene, count in genes.most_common(10):
            lines.append(f"- {gene}: {count}")
    else:
        lines.append("- No gene-level annotations were returned.")

    lines.extend(
        [
            "",
            "### Population Frequency Context",
            "",
        ]
    )
    if population_outliers:
        for row in population_outliers[:5]:
            lines.append(f"- {row['input_id']}: {summarize_population_context(row)}")
    else:
        lines.append("- No strong population-frequency outliers were detected.")

    lines.extend(
        [
            "",
            "## Interpretation Caveats",
            "",
            "- `Clinically relevant` is a narrow reporting rule in this skill, not a clinical diagnosis or ACMG-style classification.",
            "- ClinVar strings may contain mixed or conflicting labels; review the underlying accessions before acting on a summary bucket.",
            "- Population rarity depends on dataset coverage and ancestry context, so frequency outliers should be interpreted with care.",
            "- Variants with missing gnomAD values are not promoted to the clinically relevant set by the current rule, even if ClinVar evidence is strong.",
        ]
    )

    if figure_path is not None:
        lines.extend(
            [
                "",
                "## Summary Figure",
                "",
                f"- Summary figure: `figures/{figure_path.name}`",
                "- The figure is supplementary; the text tables above remain the primary machine-readable summary.",
            ]
        )

    lines.extend(
        [
            "",
            "## Methods",
            "",
            "- Parsed the VCF with `pysam.VariantFile`.",
            f"- Submitted variants to Ensembl VEP in batches of {metadata['batch_size']} using `{VEP_ENDPOINT}`.",
            "- Prioritized the most severe consequence per variant across transcript and other feature annotations.",
            "- Normalized ClinVar labels into buckets and extracted population-aware frequency context where available.",
            "- Assigned a numeric priority score plus a human-readable tier to rank variants for reporting.",
            "- Marked variants as clinically relevant when ClinVar included a pathogenic label and gnomAD AF was below 0.001.",
            "",
            "## Notes",
            "",
        ]
    )
    if failures:
        note_bits: list[str] = []
        if annotated_missing:
            note_bits.append(f"{len(annotated_missing)} missing VEP annotations")
        if annotated_error:
            note_bits.append(f"{len(annotated_error)} request-level errors")
        failure_note = (
            ", ".join(note_bits) if note_bits else f"{len(failures)} failed annotations"
        )
        lines.append(f"- Observed {failure_note}; see `result.json` for details.")
    else:
        lines.append("- All variant batches completed without request-level errors.")
    if warnings:
        lines.append(
            f"- {len(warnings)} variants carried VEP warning metadata; review `result.json` if any prioritized calls seem incomplete."
        )

    return header + "\n".join(lines) + generate_report_footer()


def parse_input(input_path: Path) -> list[dict[str, Any]]:
    """Read and validate the input VCF."""
    return parse_vcf(input_path)


def analyse(
    records: list[dict[str, Any]],
    args: argparse.Namespace,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    """Run the core variant annotation workflow."""
    client = VEPClient(cache_dir=args.cache_dir, use_cache=not args.no_cache)
    return annotate_variants(
        records,
        client=client,
        assembly=args.assembly,
        batch_size=args.batch_size,
    )


def generate_report(
    input_path: Path,
    output_dir: Path,
    annotations: list[dict[str, Any]],
    failures: list[dict[str, Any]],
    metadata: dict[str, Any],
    tsv_path: Path,
    figure_path: Path | None = None,
) -> Path:
    """Write the markdown report and return its path."""
    report_text = generate_markdown_report(
        input_path=input_path,
        output_dir=output_dir,
        annotations=annotations,
        failures=failures,
        metadata=metadata,
        tsv_path=tsv_path,
        figure_path=figure_path,
    )
    report_path = output_dir / "report.md"
    report_path.write_text(report_text)
    return report_path


def write_reproducibility_bundle(
    output_dir: Path,
    args: argparse.Namespace,
    input_path: Path,
) -> Path:
    repro_dir = output_dir / "reproducibility"
    repro_dir.mkdir(parents=True, exist_ok=True)

    parts = ["python", Path(__file__).name]
    if args.demo:
        parts.append("--demo")
    else:
        parts.extend(["--input", str(input_path)])
    parts.extend(["--output", str(output_dir)])
    parts.extend(["--assembly", args.assembly, "--batch-size", str(args.batch_size)])
    if args.no_cache:
        parts.append("--no-cache")
    elif args.cache_dir:
        parts.extend(["--cache-dir", str(args.cache_dir)])

    commands_path = repro_dir / "commands.sh"
    command_text = (
        "#!/usr/bin/env bash\n"
        f"# Generated: {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}\n"
        f"{' '.join(shlex.quote(part) for part in parts)}\n"
    )
    commands_path.write_text(command_text)

    root_commands_path = output_dir / "commands.sh"
    root_commands_path.write_text(command_text)
    return commands_path


def run_pipeline(args: argparse.Namespace) -> dict[str, Any]:
    input_path = DEFAULT_DEMO_INPUT if args.demo else Path(args.input)
    if not input_path.exists():
        raise FileNotFoundError(f"Input VCF not found: {input_path}")

    output_dir = Path(args.output)
    if output_dir.exists() and any(output_dir.iterdir()):
        print(
            f"Warning: output directory already exists and contains files: {output_dir}"
        )
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Parsing VCF with pysam: {input_path}")
    records = parse_input(input_path)
    print(f"  Parsed {len(records)} variant records")

    print("Annotating variants via Ensembl VEP...")
    annotations, failures, metadata = analyse(records, args)

    print("Writing annotated TSV...")
    tsv_path = write_tsv(output_dir, annotations)

    print("Writing summary figure...")
    summary_data = summarize_annotations(annotations)
    figure_path = plot_summary_figure(summary_data, output_dir)

    print("Writing markdown report...")
    report_path = generate_report(
        input_path=input_path,
        output_dir=output_dir,
        annotations=annotations,
        failures=failures,
        metadata=metadata,
        tsv_path=tsv_path,
        figure_path=figure_path,
    )

    print("Writing reproducibility bundle...")
    commands_path = write_reproducibility_bundle(output_dir, args, input_path)

    clinically_relevant = summary_data["clinically_relevant"]
    pathogenic = summary_data["pathogenic"]
    genes = sorted({row["gene"] for row in annotations if row.get("gene")})
    top_variants = summary_data["top_variants"]
    impact_counts = summary_data["impact_counts"]
    tier_counts = summary_data["tier_counts"]
    clinvar_counts = summary_data["clinvar_counts"]
    input_checksum = sha256_file(input_path)

    print("Writing result.json...")
    write_result_json(
        output_dir=output_dir,
        skill=SKILL_NAME,
        version=SKILL_VERSION,
        input_checksum=input_checksum,
        summary={
            "variants_in_vcf": len(records),
            "variants_annotated": len(annotations),
            "batches_sent": metadata["batches_sent"],
            "failed_batches": metadata["failed_batches"],
            "clinically_relevant_findings": len(clinically_relevant),
            "tier_1_variants": tier_counts.get("Tier 1", 0),
            "high_impact_variants": impact_counts.get("HIGH", 0),
            "pathogenic_or_likely_pathogenic": len(pathogenic),
            "genes_hit": len(genes),
            "assembly": metadata["assembly"],
        },
        data={
            "annotations": annotations,
            "top_variants": top_variants,
            "flagged_findings": clinically_relevant,
            "summary_tables": {
                "priority_buckets": dict(tier_counts),
                "impact_tiers": dict(impact_counts),
                "clinvar_buckets": dict(clinvar_counts),
            },
            "failed_variants": failures,
            "api_metadata": metadata,
            "files": {
                "report": str(report_path.name),
                "annotated_tsv": str(tsv_path.relative_to(output_dir)),
                "commands": str(commands_path.relative_to(output_dir)),
                "root_commands": "commands.sh",
                "summary_figure": str(figure_path.relative_to(output_dir))
                if figure_path
                else "",
            },
            "disclaimer": DISCLAIMER,
        },
    )

    return {
        "input_path": str(input_path),
        "output_dir": str(output_dir),
        "report_path": str(report_path),
        "annotated_tsv": str(tsv_path),
        "variants": len(annotations),
        "clinically_relevant": len(clinically_relevant),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Annotate VCF variants with Ensembl VEP and summarize clinically relevant findings."
    )
    parser.add_argument("--input", type=Path, help="Input VCF path")
    parser.add_argument(
        "--output",
        "-o",
        type=Path,
        default=Path("variant_annotation_report"),
        help="Output directory",
    )
    parser.add_argument(
        "--demo", action="store_true", help="Run the bundled synthetic VCF demo"
    )
    parser.add_argument(
        "--assembly",
        default=DEFAULT_ASSEMBLY,
        help="Genome assembly to use (default: GRCh38)",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=DEFAULT_BATCH_SIZE,
        help="Variants per VEP POST request",
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=DEFAULT_CACHE_DIR,
        help="Local cache directory for VEP responses",
    )
    parser.add_argument(
        "--no-cache", action="store_true", help="Disable local VEP response caching"
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if not args.demo and not args.input:
        print("ERROR: Provide --input <variants.vcf> or --demo", file=sys.stderr)
        sys.exit(1)
    if args.batch_size <= 0:
        print("ERROR: --batch-size must be greater than 0", file=sys.stderr)
        sys.exit(1)

    try:
        result = run_pipeline(args)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)

    print()
    print(f"Report written to {result['report_path']}")
    print(f"Annotated TSV written to {result['annotated_tsv']}")
    print(f"Clinically relevant findings: {result['clinically_relevant']}")


if __name__ == "__main__":
    main()
