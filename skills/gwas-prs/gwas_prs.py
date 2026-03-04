#!/usr/bin/env python3
"""
GWAS Polygenic Risk Score Calculator — PGS Catalog integration.

Calculates Polygenic Risk Scores (PRS) from DTC genetic data (23andMe/AncestryDNA)
using the PGS Catalog REST API and pre-curated scoring files.

Usage:
    python gwas_prs.py --input <23andme_file> --trait "type 2 diabetes" --output <dir>
    python gwas_prs.py --input <23andme_file> --pgs-id PGS000013 --output <dir>
    python gwas_prs.py --demo --output /tmp/prs_demo
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import requests

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

PGS_API_BASE = "https://www.pgscatalog.org/rest"
RATE_LIMIT_INTERVAL = 0.55  # seconds between requests (stay under 2 req/sec)
CACHE_TTL = 86400  # 24 hours
USER_AGENT = "ClawBio-GWAS-PRS/0.1.0"

DISCLAIMER = (
    "ClawBio is a research and educational tool. It is not a medical device "
    "and does not provide clinical diagnoses. Consult a healthcare "
    "professional before making any medical decisions."
)

# Risk category thresholds (percentile-based)
RISK_CATEGORIES = [
    (25, "Low"),
    (75, "Average"),
    (95, "Elevated"),
    (100, "High"),
]

SKILL_DIR = Path(__file__).resolve().parent
DATA_DIR = SKILL_DIR / "data"

# ---------------------------------------------------------------------------
# Curated scores — offline demo scoring with reference distributions
# ---------------------------------------------------------------------------

CURATED_SCORES: dict[str, dict] = {
    "PGS000013": {
        "name": "PGS000013",
        "trait": "Type 2 diabetes",
        "variants_count": 8,
        "publication": "Vassy et al. (2014) Ann Intern Med",
        "reference_distribution": {
            "mean": 1.12,
            "sd": 0.30,
            "population": "EUR",
        },
    },
    "PGS000011": {
        "name": "PGS000011",
        "trait": "Atrial fibrillation",
        "variants_count": 12,
        "publication": "Tada et al. (2014) Circ Cardiovasc Genet",
        "reference_distribution": {
            "mean": 0.65,
            "sd": 0.23,
            "population": "EUR",
        },
    },
    "PGS000004": {
        "name": "PGS000004",
        "trait": "Coronary artery disease",
        "variants_count": 46,
        "publication": "Abraham et al. (2016) Eur Heart J",
        "reference_distribution": {
            "mean": 2.84,
            "sd": 0.28,
            "population": "EUR",
        },
    },
    "PGS000001": {
        "name": "PGS000001",
        "trait": "Breast cancer",
        "variants_count": 77,
        "publication": "Mavaddat et al. (2015) J Natl Cancer Inst",
        "reference_distribution": {
            "mean": 4.23,
            "sd": 0.54,
            "population": "EUR",
        },
    },
    "PGS000057": {
        "name": "PGS000057",
        "trait": "Prostate cancer",
        "variants_count": 147,
        "publication": "Schumacher et al. (2018) Nat Genet",
        "reference_distribution": {
            "mean": 7.11,
            "sd": 0.56,
            "population": "EUR",
        },
    },
    "PGS000039": {
        "name": "PGS000039",
        "trait": "BMI",
        "variants_count": 97,
        "publication": "Locke et al. (2015) Nature",
        "reference_distribution": {
            "mean": 2.89,
            "sd": 0.25,
            "population": "EUR",
        },
    },
}


# ---------------------------------------------------------------------------
# PGS Catalog API Client
# ---------------------------------------------------------------------------


class PGSCatalogClient:
    """Rate-limited, caching HTTP client for the PGS Catalog REST API."""

    def __init__(self, cache_dir: Path, use_cache: bool = True):
        self.cache_dir = cache_dir
        self.use_cache = use_cache
        self._last_request_time = 0.0
        self.session = requests.Session()
        self.session.headers.update({
            "Accept": "application/json",
            "User-Agent": USER_AGENT,
        })
        if use_cache:
            cache_dir.mkdir(parents=True, exist_ok=True)

    # --- Rate limiting ---

    def _throttle(self):
        """Enforce minimum interval between API requests."""
        elapsed = time.time() - self._last_request_time
        if elapsed < RATE_LIMIT_INTERVAL:
            time.sleep(RATE_LIMIT_INTERVAL - elapsed)
        self._last_request_time = time.time()

    # --- Caching ---

    def _cache_key(self, endpoint: str, params: dict) -> str:
        """Generate a short SHA256 cache key from endpoint + params."""
        raw = f"{endpoint}|{json.dumps(params, sort_keys=True)}"
        return hashlib.sha256(raw.encode()).hexdigest()[:16]

    def _get_cached(self, key: str) -> Optional[dict]:
        """Return cached response if present and not expired."""
        path = self.cache_dir / f"{key}.json"
        if not path.exists():
            return None
        try:
            data = json.loads(path.read_text())
            if time.time() - data.get("_cached_at", 0) < CACHE_TTL:
                return data.get("response")
        except (json.JSONDecodeError, KeyError):
            pass
        return None

    def _set_cached(self, key: str, response_data):
        """Write response data to cache with timestamp."""
        path = self.cache_dir / f"{key}.json"
        path.write_text(json.dumps({
            "_cached_at": time.time(),
            "response": response_data,
        }, indent=2))

    # --- Core request ---

    def _request(self, endpoint: str, params: dict | None = None) -> dict:
        """Core GET request with caching, rate limiting, and 429 retry."""
        params = params or {}
        cache_key = self._cache_key(endpoint, params)

        if self.use_cache:
            cached = self._get_cached(cache_key)
            if cached is not None:
                return cached

        self._throttle()
        url = f"{PGS_API_BASE}/{endpoint}"
        resp = self.session.get(url, params=params, timeout=30)

        if resp.status_code == 429:
            retry_after = float(resp.headers.get("Retry-After", "5"))
            print(f"  Rate limited (429), waiting {retry_after}s...")
            time.sleep(retry_after)
            resp = self.session.get(url, params=params, timeout=30)

        resp.raise_for_status()
        data = resp.json()

        if self.use_cache:
            self._set_cached(cache_key, data)

        return data

    # --- Public API methods ---

    def search_traits(self, term: str) -> list[dict]:
        """Search PGS Catalog traits by keyword.

        Returns list of trait dicts with 'id' (EFO ID) and 'label'.
        """
        data = self._request("trait/search", {"term": term})
        results = data.get("results", [])
        traits = []
        for item in results:
            traits.append({
                "id": item.get("id", ""),
                "label": item.get("label", ""),
                "description": item.get("description", ""),
                "associated_pgs_ids": item.get("associated_pgs_ids", []),
            })
        return traits

    def get_score_metadata(self, pgs_id: str) -> dict:
        """Get metadata for a specific PGS score.

        Returns dict with score metadata including trait, variant count, etc.
        """
        data = self._request(f"score/{pgs_id}")
        return data

    def search_scores_by_trait(self, trait_id: str) -> list[dict]:
        """Search scores associated with an EFO trait ID.

        Handles API pagination (the API returns 'next' URL).
        Returns list of score metadata dicts.
        """
        all_scores = []
        endpoint = "score/search"
        params = {"trait_id": trait_id}

        while True:
            data = self._request(endpoint, params)
            results = data.get("results", [])
            all_scores.extend(results)

            next_url = data.get("next")
            if not next_url:
                break

            # Parse the next URL to extract endpoint and params
            # The 'next' field is a full URL, so we need to request it directly
            self._throttle()
            resp = self.session.get(next_url, timeout=30)
            if resp.status_code == 429:
                retry_after = float(resp.headers.get("Retry-After", "5"))
                time.sleep(retry_after)
                resp = self.session.get(next_url, timeout=30)
            resp.raise_for_status()
            data = resp.json()
            results = data.get("results", [])
            all_scores.extend(results)
            next_url = data.get("next")
            if not next_url:
                break
            # Continue loop with the new next URL — re-encode as params
            # to keep the caching/throttle flow, we break and use direct URL
            # Actually, we need to keep paginating with direct URLs
            while next_url:
                self._throttle()
                resp = self.session.get(next_url, timeout=30)
                if resp.status_code == 429:
                    retry_after = float(resp.headers.get("Retry-After", "5"))
                    time.sleep(retry_after)
                    resp = self.session.get(next_url, timeout=30)
                resp.raise_for_status()
                data = resp.json()
                all_scores.extend(data.get("results", []))
                next_url = data.get("next")
            break

        return all_scores

    def download_scoring_file(
        self, pgs_id: str, build: str = "GRCh37"
    ) -> Path:
        """Download the harmonized scoring file from EBI FTP.

        URL pattern:
            https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/{pgs_id}/
            ScoringFiles/Harmonized/{pgs_id}_hmPOS_{build}.txt.gz

        The gzipped file is cached locally. Returns the local file path.
        """
        filename = f"{pgs_id}_hmPOS_{build}.txt.gz"
        local_path = self.cache_dir / filename

        if self.use_cache and local_path.exists():
            # Check file age for cache TTL
            age = time.time() - local_path.stat().st_mtime
            if age < CACHE_TTL:
                return local_path

        url = (
            f"https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/"
            f"{pgs_id}/ScoringFiles/Harmonized/{filename}"
        )
        self._throttle()
        print(f"  Downloading {filename}...")
        resp = self.session.get(url, timeout=120, stream=True)

        if resp.status_code == 429:
            retry_after = float(resp.headers.get("Retry-After", "5"))
            print(f"  Rate limited (429), waiting {retry_after}s...")
            time.sleep(retry_after)
            resp = self.session.get(url, timeout=120, stream=True)

        resp.raise_for_status()

        self.cache_dir.mkdir(parents=True, exist_ok=True)
        with open(local_path, "wb") as f:
            for chunk in resp.iter_content(chunk_size=65536):
                f.write(chunk)

        return local_path


# ---------------------------------------------------------------------------
# Genotype parser — 23andMe / AncestryDNA / generic
# ---------------------------------------------------------------------------


def detect_format(lines: list[str]) -> str:
    """Detect genotype file format from the first 20 lines.

    Returns: '23andme', 'ancestrydna', 'generic', or 'unknown'.
    """
    for line in lines[:20]:
        if line.startswith("# rsid"):
            return "23andme"
        if "RSID" in line and "CHROMOSOME" in line:
            return "ancestrydna"
        if "rsid" in line and "chromosome" in line:
            return "generic"
    # Fallback: inspect data lines
    for line in lines[:20]:
        if line.startswith("#") or line.strip() == "":
            continue
        parts = line.split("\t")
        if len(parts) >= 4 and parts[0].startswith("rs"):
            return "23andme"
        parts = line.split(",")
        if len(parts) >= 4 and parts[0].startswith("rs"):
            return "ancestrydna"
    return "unknown"


def parse_genotype_file(path: str | Path) -> tuple[str, int, dict[str, str]]:
    """Parse a DTC genotype file, returning ALL rsIDs.

    Supports .gz files, 23andMe (tab-separated) and AncestryDNA
    (comma-separated) formats.

    Returns:
        (format_name, total_snps, {rsid: genotype_string})
    """
    path = Path(path)
    if str(path).endswith(".gz"):
        with gzip.open(path, "rt", errors="replace") as fh:
            content = fh.read()
    else:
        content = path.read_text(errors="replace")

    lines = content.split("\n")
    fmt = detect_format(lines)

    snps: dict[str, str] = {}
    for line in lines:
        if line.startswith("#") or line.strip() == "":
            continue
        # Skip header lines
        if "rsid" in line.lower() and "chromosome" in line.lower():
            continue

        parts = line.split("\t") if "\t" in line else line.split(",")
        if len(parts) >= 4:
            rsid = parts[0].strip()
            if not rsid.startswith("rs"):
                continue
            # AncestryDNA has allele1, allele2 in columns 3 and 4
            if len(parts) == 5:
                genotype = parts[3].strip() + parts[4].strip()
            else:
                genotype = parts[3].strip()
            if genotype and genotype not in ("--", "00", ""):
                snps[rsid] = genotype.upper()

    return fmt, len(snps), snps


# ---------------------------------------------------------------------------
# Scoring file parser — PGS Catalog harmonized TSV
# ---------------------------------------------------------------------------


def parse_scoring_file(filepath: str | Path) -> list[dict]:
    """Parse a PGS Catalog harmonized scoring file (TSV, may be gzipped).

    The file has comment lines starting with '#' (including metadata like
    #pgs_id=...), then a header line, then data rows.

    Key columns: hm_rsID (or rsID), effect_allele, effect_weight,
    optionally other_allele, allelefrequency_effect.

    Returns:
        List of dicts with keys: rsid, effect_allele, effect_weight,
        other_allele (optional), allele_freq (optional).
    """
    filepath = Path(filepath)

    if str(filepath).endswith(".gz"):
        with gzip.open(filepath, "rt", errors="replace") as fh:
            raw_lines = fh.readlines()
    else:
        raw_lines = filepath.read_text(errors="replace").splitlines(keepends=True)

    # Separate comment lines from data
    header_line = None
    data_lines = []
    for line in raw_lines:
        stripped = line.strip()
        if stripped.startswith("#"):
            continue
        if header_line is None:
            header_line = stripped
            continue
        if stripped:
            data_lines.append(stripped)

    if header_line is None:
        return []

    # Parse header — find relevant column indices
    headers = header_line.split("\t")
    col_map: dict[str, int] = {}
    for i, h in enumerate(headers):
        col_map[h.strip().lower()] = i

    # Identify rsID column (prefer hm_rsID over rsID)
    rsid_col = col_map.get("hm_rsid", col_map.get("rsid"))
    effect_allele_col = col_map.get("effect_allele")
    weight_col = col_map.get("effect_weight")
    other_allele_col = col_map.get("other_allele", col_map.get("hm_inferOtherAllele"))
    freq_col = col_map.get("allelefrequency_effect")

    if rsid_col is None or effect_allele_col is None or weight_col is None:
        print(f"  WARNING: Scoring file missing required columns. "
              f"Found: {list(col_map.keys())}")
        return []

    variants = []
    for line in data_lines:
        cols = line.split("\t")
        if len(cols) <= max(rsid_col, effect_allele_col, weight_col):
            continue

        rsid = cols[rsid_col].strip()
        # Skip rows without a valid rsID
        if not rsid or not rsid.startswith("rs"):
            continue

        effect_allele = cols[effect_allele_col].strip().upper()
        try:
            weight = float(cols[weight_col].strip())
        except (ValueError, IndexError):
            continue

        variant: dict = {
            "rsid": rsid,
            "effect_allele": effect_allele,
            "effect_weight": weight,
        }

        # Optional: other allele
        if other_allele_col is not None and other_allele_col < len(cols):
            oa = cols[other_allele_col].strip().upper()
            if oa:
                variant["other_allele"] = oa

        # Optional: allele frequency
        if freq_col is not None and freq_col < len(cols):
            try:
                freq = float(cols[freq_col].strip())
                variant["allele_freq"] = freq
            except (ValueError, IndexError):
                pass

        variants.append(variant)

    return variants


# ---------------------------------------------------------------------------
# PRS Calculator
# ---------------------------------------------------------------------------


def compute_allele_dosage(genotype: str, effect_allele: str) -> int:
    """Compute allele dosage: count of effect allele copies in genotype.

    Args:
        genotype: Genotype string, e.g. 'AG', 'AA', 'A' (hemizygous).
        effect_allele: Single-letter effect allele, e.g. 'A'.

    Returns:
        Dosage: 0, 1, or 2 (0 or 1 for hemizygous single-letter genotypes).
    """
    effect_allele = effect_allele.upper()
    genotype = genotype.upper()

    # Handle hemizygous (single allele, e.g. X chromosome in males)
    if len(genotype) == 1:
        return 1 if genotype == effect_allele else 0

    # Standard diploid genotype (two characters)
    count = 0
    for allele in genotype:
        if allele == effect_allele:
            count += 1
    return count


def calculate_prs(
    genotype_dict: dict[str, str],
    scoring_variants: list[dict],
) -> dict:
    """Calculate Polygenic Risk Score = sum(weight * dosage).

    Args:
        genotype_dict: {rsid: genotype_string} from parsed genotype file.
        scoring_variants: List of variant dicts from parse_scoring_file().

    Returns:
        Dict with: raw_score, variants_used, variants_total,
        variants_missing, overlap_fraction, per_variant details.
    """
    raw_score = 0.0
    variants_used = 0
    variants_missing = 0
    per_variant: list[dict] = []

    for sv in scoring_variants:
        rsid = sv["rsid"]
        effect_allele = sv["effect_allele"]
        weight = sv["effect_weight"]

        if rsid in genotype_dict:
            genotype = genotype_dict[rsid]
            dosage = compute_allele_dosage(genotype, effect_allele)
            contribution = weight * dosage
            raw_score += contribution
            variants_used += 1
            per_variant.append({
                "rsid": rsid,
                "effect_allele": effect_allele,
                "genotype": genotype,
                "dosage": dosage,
                "weight": weight,
                "contribution": contribution,
                "status": "scored",
            })
        else:
            variants_missing += 1
            per_variant.append({
                "rsid": rsid,
                "effect_allele": effect_allele,
                "genotype": None,
                "dosage": None,
                "weight": weight,
                "contribution": 0.0,
                "status": "missing",
            })

    variants_total = len(scoring_variants)
    overlap = variants_used / variants_total if variants_total > 0 else 0.0

    return {
        "raw_score": raw_score,
        "variants_used": variants_used,
        "variants_total": variants_total,
        "variants_missing": variants_missing,
        "overlap_fraction": overlap,
        "per_variant": per_variant,
    }


# ---------------------------------------------------------------------------
# Percentile Estimator
# ---------------------------------------------------------------------------


def _percentile_from_z(z: float) -> float:
    """Compute the CDF of the standard normal distribution at z.

    Uses math.erf for an exact closed-form calculation.

    Returns:
        Percentile as a value between 0 and 100.
    """
    return 50.0 * (1.0 + math.erf(z / math.sqrt(2.0)))


def _assign_risk_category(percentile: float) -> str:
    """Map a percentile (0-100) to a risk category label."""
    for threshold, label in RISK_CATEGORIES:
        if percentile <= threshold:
            return label
    return "High"


def estimate_percentile(
    raw_score: float,
    pgs_id: str,
    scoring_variants: list[dict],
) -> dict:
    """Estimate percentile rank using a tiered approach.

    Tier 1: Use curated reference distribution (mean/SD) if available.
    Tier 2: Estimate from allele frequencies in the scoring file.
    Tier 3: Return None (unavailable).

    Returns:
        Dict with: percentile (float or None), method (str),
        risk_category (str or None), z_score (float or None),
        reference_population (str or None).
    """
    result: dict = {
        "percentile": None,
        "method": "unavailable",
        "risk_category": None,
        "z_score": None,
        "reference_population": None,
    }

    # Tier 1: Curated reference distribution
    curated = CURATED_SCORES.get(pgs_id)
    if curated and "reference_distribution" in curated:
        ref = curated["reference_distribution"]
        mean = ref["mean"]
        sd = ref["sd"]
        if sd > 0:
            z = (raw_score - mean) / sd
            pct = _percentile_from_z(z)
            result["percentile"] = round(pct, 1)
            result["method"] = "curated_reference"
            result["risk_category"] = _assign_risk_category(pct)
            result["z_score"] = round(z, 3)
            result["reference_population"] = ref.get("population", "EUR")
            return result

    # Tier 2: Estimate from allele frequencies in scoring file
    freqs_available = [
        sv for sv in scoring_variants
        if "allele_freq" in sv and 0 < sv["allele_freq"] < 1
    ]
    if len(freqs_available) >= 3:
        # Expected score under HWE: E[PRS] = sum(2 * freq * weight)
        # Variance: Var[PRS] = sum(2 * freq * (1-freq) * weight^2)
        expected = 0.0
        variance = 0.0
        for sv in freqs_available:
            freq = sv["allele_freq"]
            w = sv["effect_weight"]
            expected += 2.0 * freq * w
            variance += 2.0 * freq * (1.0 - freq) * (w ** 2)

        if variance > 0:
            sd = math.sqrt(variance)
            z = (raw_score - expected) / sd
            pct = _percentile_from_z(z)
            result["percentile"] = round(pct, 1)
            result["method"] = "allele_frequency_estimate"
            result["risk_category"] = _assign_risk_category(pct)
            result["z_score"] = round(z, 3)
            result["reference_population"] = "estimated (allele freq)"
            return result

    # Tier 3: Cannot estimate
    result["method"] = "unavailable"
    return result


# ---------------------------------------------------------------------------
# Report Generator
# ---------------------------------------------------------------------------


def generate_report(
    results: list[dict],
    input_info: dict,
    args: argparse.Namespace,
) -> str:
    """Generate a Markdown report from PRS results.

    Args:
        results: List of per-score result dicts, each containing:
            pgs_id, trait, prs (from calculate_prs), percentile_info
            (from estimate_percentile), metadata, scoring_variants.
        input_info: Dict with format, total_snps, filepath.
        args: CLI args namespace.

    Returns:
        Markdown string.
    """
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    lines = [
        "# GWAS Polygenic Risk Score Report",
        "",
        f"**Date**: {now}",
        f"**Tool**: ClawBio GWAS-PRS v0.1.0",
        f"**Input file**: {input_info.get('filepath', 'N/A')}",
        f"**Format detected**: {input_info.get('format', 'N/A')}",
        f"**Total SNPs in file**: {input_info.get('total_snps', 'N/A'):,}",
        f"**Genome build**: {args.build}",
        "",
        "---",
        "",
    ]

    # ----- Summary table -----
    lines.append("## Summary")
    lines.append("")
    lines.append(
        "| PGS ID | Trait | Raw Score | Percentile | Risk Category | Overlap |"
    )
    lines.append(
        "|--------|-------|-----------|------------|---------------|---------|"
    )

    for r in results:
        pgs_id = r["pgs_id"]
        trait = r.get("trait", "Unknown")
        prs = r["prs"]
        pct_info = r["percentile_info"]
        raw = f"{prs['raw_score']:.4f}"
        pct = (
            f"{pct_info['percentile']:.1f}%"
            if pct_info["percentile"] is not None
            else "N/A"
        )
        risk = pct_info.get("risk_category") or "N/A"
        overlap = f"{prs['overlap_fraction'] * 100:.1f}%"
        lines.append(
            f"| {pgs_id} | {trait} | {raw} | {pct} | {risk} | {overlap} |"
        )

    lines.append("")

    # ----- Overlap warnings -----
    low_overlap = [r for r in results if r["prs"]["overlap_fraction"] < 0.70]
    if low_overlap:
        lines.append("### Overlap Warnings")
        lines.append("")
        lines.append(
            "> **Warning**: The following scores have < 70% SNP overlap with "
            "your genotype file. Results should be interpreted with caution."
        )
        lines.append("")
        for r in low_overlap:
            prs = r["prs"]
            lines.append(
                f"- **{r['pgs_id']}** ({r.get('trait', '')}): "
                f"{prs['variants_used']}/{prs['variants_total']} variants "
                f"({prs['overlap_fraction'] * 100:.1f}% overlap)"
            )
        lines.append("")

    # ----- Per-score detail sections -----
    lines.append("---")
    lines.append("")
    lines.append("## Score Details")
    lines.append("")

    for r in results:
        pgs_id = r["pgs_id"]
        trait = r.get("trait", "Unknown")
        prs = r["prs"]
        pct_info = r["percentile_info"]
        metadata = r.get("metadata", {})

        lines.append(f"### {pgs_id} — {trait}")
        lines.append("")

        # Metadata table
        lines.append("| Property | Value |")
        lines.append("|----------|-------|")
        lines.append(f"| **PGS ID** | {pgs_id} |")
        lines.append(f"| **Trait** | {trait} |")
        if metadata.get("publication"):
            lines.append(f"| **Publication** | {metadata['publication']} |")
        lines.append(
            f"| **Variants in score** | {prs['variants_total']} |"
        )
        lines.append(
            f"| **Variants matched** | {prs['variants_used']} |"
        )
        lines.append(
            f"| **Variants missing** | {prs['variants_missing']} |"
        )
        lines.append(
            f"| **Overlap** | {prs['overlap_fraction'] * 100:.1f}% |"
        )
        lines.append(f"| **Raw PRS** | {prs['raw_score']:.6f} |")
        if pct_info["percentile"] is not None:
            lines.append(
                f"| **Percentile** | {pct_info['percentile']:.1f}% |"
            )
            lines.append(
                f"| **Risk category** | {pct_info['risk_category']} |"
            )
            lines.append(
                f"| **Z-score** | {pct_info['z_score']:.3f} |"
            )
            lines.append(
                f"| **Estimation method** | {pct_info['method']} |"
            )
            if pct_info.get("reference_population"):
                lines.append(
                    f"| **Reference population** | "
                    f"{pct_info['reference_population']} |"
                )
        else:
            lines.append("| **Percentile** | N/A (no reference available) |")
        lines.append("")

        # Variant breakdown (show scored variants, limit to top contributors)
        scored = [
            v for v in prs["per_variant"] if v["status"] == "scored"
        ]
        if scored:
            # Sort by absolute contribution (descending)
            scored_sorted = sorted(
                scored, key=lambda v: abs(v["contribution"]), reverse=True
            )
            show_count = min(len(scored_sorted), 20)
            lines.append(
                f"**Top {show_count} contributing variants** "
                f"(of {len(scored_sorted)} scored):"
            )
            lines.append("")
            lines.append(
                "| rsID | Genotype | Effect Allele | Dosage | Weight | "
                "Contribution |"
            )
            lines.append(
                "|------|----------|---------------|--------|--------|"
                "--------------|"
            )
            for v in scored_sorted[:show_count]:
                lines.append(
                    f"| {v['rsid']} | {v['genotype']} | "
                    f"{v['effect_allele']} | {v['dosage']} | "
                    f"{v['weight']:.4f} | {v['contribution']:.4f} |"
                )
            if len(scored_sorted) > show_count:
                lines.append(
                    f"| ... | | | | | "
                    f"({len(scored_sorted) - show_count} more variants) |"
                )
            lines.append("")
        else:
            lines.append("*No variants from this score were found in your "
                         "genotype file.*")
            lines.append("")

    # ----- Methods -----
    lines.append("---")
    lines.append("")
    lines.append("## Methods")
    lines.append("")
    lines.append(
        "Polygenic Risk Scores (PRS) were calculated as the weighted sum of "
        "effect allele dosages across all available variants:"
    )
    lines.append("")
    lines.append("```")
    lines.append("PRS = sum(effect_weight_i * dosage_i)")
    lines.append("```")
    lines.append("")
    lines.append(
        "Where `dosage_i` is the number of copies (0, 1, or 2) of the effect "
        "allele at variant *i*, and `effect_weight_i` is the log odds ratio "
        "(or beta coefficient) from the original GWAS."
    )
    lines.append("")
    lines.append(
        "Percentile estimates use one of two methods: (1) curated EUR "
        "reference distributions (mean and SD from published cohorts), or "
        "(2) population-level allele frequency estimates from the PGS Catalog "
        "scoring file to compute expected score and variance under "
        "Hardy-Weinberg equilibrium."
    )
    lines.append("")
    lines.append("**Scoring files**: PGS Catalog "
                 "(https://www.pgscatalog.org/)")
    lines.append("")
    lines.append("**Genome build**: " + args.build)
    lines.append("")

    # ----- Limitations -----
    lines.append("## Limitations")
    lines.append("")
    lines.append(
        "- PRS calculated from DTC genotype data typically covers only a "
        "subset of the variants in the full score, reducing predictive power."
    )
    lines.append(
        "- Percentile estimates assume a European (EUR) reference population "
        "unless otherwise noted. PRS accuracy and transferability vary "
        "significantly across ancestries."
    )
    lines.append(
        "- DTC genotype arrays may not include all variants in a PGS, "
        "especially rare variants or indels."
    )
    lines.append(
        "- Risk categories are simplified and do not account for family "
        "history, lifestyle, or environmental factors."
    )
    lines.append("")

    # ----- Disclaimer -----
    lines.append("## Disclaimer")
    lines.append("")
    lines.append(f"*{DISCLAIMER}*")
    lines.append("")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description="GWAS Polygenic Risk Score Calculator (PGS Catalog)"
    )
    parser.add_argument(
        "--input",
        help="Path to 23andMe/AncestryDNA genotype file",
    )
    parser.add_argument(
        "--trait",
        help="Search PGS Catalog by trait name (e.g., 'type 2 diabetes')",
    )
    parser.add_argument(
        "--pgs-id",
        help="Specific PGS Catalog score ID (e.g., PGS000013)",
    )
    parser.add_argument(
        "--output", "-o",
        help="Output directory for report",
    )
    parser.add_argument(
        "--demo",
        action="store_true",
        help="Run demo with synthetic patient + curated scores",
    )
    parser.add_argument(
        "--min-overlap",
        type=float,
        default=0.5,
        help="Minimum SNP overlap fraction (default: 0.5)",
    )
    parser.add_argument(
        "--max-variants",
        type=int,
        default=50000,
        help="Skip scores with more variants (default: 50000)",
    )
    parser.add_argument(
        "--build",
        default="GRCh37",
        choices=["GRCh37", "GRCh38"],
        help="Genome build (default: GRCh37)",
    )
    parser.add_argument(
        "--no-cache",
        action="store_true",
        help="Skip cache",
    )
    parser.add_argument(
        "--cache-dir",
        default=str(Path.home() / ".clawbio" / "pgs_cache"),
        help="Cache directory (default: ~/.clawbio/pgs_cache/)",
    )

    args = parser.parse_args()

    # Validate: need at least one of --demo, --trait, or --pgs-id
    if not any([args.demo, args.trait, args.pgs_id]):
        parser.print_help()
        print(
            "\nError: provide at least one of --demo, --trait, or --pgs-id"
        )
        sys.exit(1)

    # Validate: --trait and --pgs-id require --input (unless --demo)
    if not args.demo and not args.input:
        parser.print_help()
        print("\nError: --input is required when not using --demo")
        sys.exit(1)

    cache_dir = Path(args.cache_dir)
    use_cache = not args.no_cache
    client = PGSCatalogClient(cache_dir=cache_dir, use_cache=use_cache)

    # -----------------------------------------------------------
    # Step 1: Determine which scoring files to use
    # -----------------------------------------------------------
    scoring_files: list[dict] = []
    # Each entry: {pgs_id, trait, filepath, metadata}

    if args.demo:
        # Demo mode: use built-in demo patient + curated scores
        demo_patient = SKILL_DIR / "demo_patient_prs.txt"
        if not demo_patient.exists():
            print(f"Error: demo patient file not found at {demo_patient}")
            sys.exit(1)
        args.input = str(demo_patient)

        print("GWAS-PRS: running demo with curated scores (no API calls)")
        print(f"  Demo patient: {demo_patient}")
        print()

        for pgs_id, meta in CURATED_SCORES.items():
            scoring_path = DATA_DIR / f"{pgs_id}_hmPOS_{args.build}.txt"
            # Also check for gzipped version
            scoring_path_gz = DATA_DIR / f"{pgs_id}_hmPOS_{args.build}.txt.gz"
            if scoring_path.exists():
                fpath = scoring_path
            elif scoring_path_gz.exists():
                fpath = scoring_path_gz
            else:
                print(f"  WARNING: Scoring file not found for {pgs_id}, "
                      f"skipping")
                continue
            scoring_files.append({
                "pgs_id": pgs_id,
                "trait": meta["trait"],
                "filepath": fpath,
                "metadata": {
                    "publication": meta.get("publication", ""),
                    "variants_count": meta.get("variants_count", 0),
                },
            })
            print(f"  Loaded curated score: {pgs_id} ({meta['trait']})")

        print()

    elif args.pgs_id:
        # Specific PGS ID mode
        pgs_id = args.pgs_id.strip().upper()
        if not pgs_id.startswith("PGS"):
            pgs_id = "PGS" + pgs_id.lstrip("0")

        print(f"GWAS-PRS: fetching score {pgs_id}")
        print()

        # Check if we have it pre-downloaded in data/
        local_path = DATA_DIR / f"{pgs_id}_hmPOS_{args.build}.txt"
        local_path_gz = DATA_DIR / f"{pgs_id}_hmPOS_{args.build}.txt.gz"

        if local_path.exists():
            filepath = local_path
            print(f"  Using pre-downloaded scoring file: {filepath}")
        elif local_path_gz.exists():
            filepath = local_path_gz
            print(f"  Using pre-downloaded scoring file: {filepath}")
        else:
            # Fetch metadata from API
            try:
                print(f"  Fetching metadata for {pgs_id}...")
                meta = client.get_score_metadata(pgs_id)
            except requests.HTTPError as e:
                print(f"Error: could not fetch {pgs_id}: {e}")
                sys.exit(1)

            trait_name = "Unknown"
            trait_efo = meta.get("trait_efo", [])
            if trait_efo:
                trait_name = trait_efo[0].get("label", "Unknown")
            elif meta.get("trait_reported"):
                trait_name = ", ".join(meta["trait_reported"])

            variant_count = meta.get("variants_number", 0)
            print(f"  Trait: {trait_name}")
            print(f"  Variants: {variant_count}")

            if variant_count > args.max_variants:
                print(
                    f"  Skipping: {variant_count} variants exceeds "
                    f"--max-variants {args.max_variants}"
                )
                sys.exit(0)

            # Download scoring file
            try:
                filepath = client.download_scoring_file(
                    pgs_id, build=args.build
                )
            except requests.HTTPError as e:
                print(f"Error downloading scoring file for {pgs_id}: {e}")
                sys.exit(1)

        # Build metadata
        trait = "Unknown"
        publication = ""
        if pgs_id in CURATED_SCORES:
            trait = CURATED_SCORES[pgs_id]["trait"]
            publication = CURATED_SCORES[pgs_id].get("publication", "")
        else:
            # Try API metadata if we fetched it
            try:
                meta_data = client.get_score_metadata(pgs_id)
                trait_efo = meta_data.get("trait_efo", [])
                if trait_efo:
                    trait = trait_efo[0].get("label", "Unknown")
                elif meta_data.get("trait_reported"):
                    trait = ", ".join(meta_data["trait_reported"])
                pub = meta_data.get("publication", {})
                if pub:
                    publication = (
                        f"{pub.get('firstauthor', '')} "
                        f"({pub.get('date_publication', '')}) "
                        f"{pub.get('journal', '')}"
                    )
            except Exception:
                pass

        scoring_files.append({
            "pgs_id": pgs_id,
            "trait": trait,
            "filepath": filepath,
            "metadata": {
                "publication": publication,
            },
        })
        print()

    elif args.trait:
        # Trait search mode
        print(f"GWAS-PRS: searching for trait '{args.trait}'")
        print()

        try:
            traits = client.search_traits(args.trait)
        except requests.HTTPError as e:
            print(f"Error searching traits: {e}")
            sys.exit(1)

        if not traits:
            print(f"  No traits found matching '{args.trait}'")
            sys.exit(0)

        # Show matching traits
        print(f"  Found {len(traits)} matching trait(s):")
        for t in traits[:5]:
            pgs_count = len(t.get("associated_pgs_ids", []))
            print(f"    - {t['id']}: {t['label']} ({pgs_count} scores)")
        if len(traits) > 5:
            print(f"    ... and {len(traits) - 5} more")
        print()

        # Use the first (most relevant) trait
        best_trait = traits[0]
        trait_id = best_trait["id"]
        trait_label = best_trait["label"]
        print(f"  Using best match: {trait_id} ({trait_label})")
        print()

        # Search scores for this trait
        try:
            scores = client.search_scores_by_trait(trait_id)
        except requests.HTTPError as e:
            print(f"Error searching scores for {trait_id}: {e}")
            sys.exit(1)

        if not scores:
            print(f"  No PGS scores found for trait {trait_id}")
            sys.exit(0)

        print(f"  Found {len(scores)} score(s) for {trait_label}:")
        for s in scores[:10]:
            sid = s.get("id", "?")
            vcount = s.get("variants_number", "?")
            sname = s.get("name", "")
            print(f"    - {sid}: {sname} ({vcount} variants)")
        if len(scores) > 10:
            print(f"    ... and {len(scores) - 10} more")
        print()

        # Filter by max-variants and download
        for s in scores:
            sid = s.get("id", "")
            vcount = s.get("variants_number", 0)
            if vcount > args.max_variants:
                print(
                    f"  Skipping {sid}: {vcount} variants exceeds "
                    f"--max-variants {args.max_variants}"
                )
                continue

            # Check for pre-downloaded file
            local_path = DATA_DIR / f"{sid}_hmPOS_{args.build}.txt"
            local_path_gz = DATA_DIR / f"{sid}_hmPOS_{args.build}.txt.gz"

            if local_path.exists():
                filepath = local_path
            elif local_path_gz.exists():
                filepath = local_path_gz
            else:
                try:
                    filepath = client.download_scoring_file(
                        sid, build=args.build
                    )
                except requests.HTTPError as e:
                    print(f"  WARNING: could not download {sid}: {e}")
                    continue

            pub = s.get("publication", {})
            publication = ""
            if pub:
                publication = (
                    f"{pub.get('firstauthor', '')} "
                    f"({pub.get('date_publication', '')}) "
                    f"{pub.get('journal', '')}"
                )

            scoring_files.append({
                "pgs_id": sid,
                "trait": trait_label,
                "filepath": filepath,
                "metadata": {
                    "publication": publication,
                    "variants_count": vcount,
                },
            })

        if not scoring_files:
            print("  No scores passed filtering. Try increasing --max-variants.")
            sys.exit(0)

        print(f"\n  Will score {len(scoring_files)} PGS score(s)")
        print()

    # -----------------------------------------------------------
    # Step 2: Parse genotype file
    # -----------------------------------------------------------
    input_path = Path(args.input)
    if not input_path.exists():
        print(f"Error: input file not found: {input_path}")
        sys.exit(1)

    print(f"Parsing genotype file: {input_path.name}")
    fmt, total_snps, genotype_dict = parse_genotype_file(input_path)
    print(f"  Format: {fmt}")
    print(f"  Total SNPs parsed: {total_snps:,}")
    print()

    input_info = {
        "filepath": str(input_path),
        "format": fmt,
        "total_snps": total_snps,
    }

    # -----------------------------------------------------------
    # Step 3: Calculate PRS for each scoring file
    # -----------------------------------------------------------
    all_results: list[dict] = []

    for sf in scoring_files:
        pgs_id = sf["pgs_id"]
        trait = sf["trait"]
        filepath = sf["filepath"]
        metadata = sf["metadata"]

        print(f"Scoring {pgs_id} ({trait})...")

        # Parse scoring file
        scoring_variants = parse_scoring_file(filepath)
        if not scoring_variants:
            print(f"  WARNING: No valid variants in scoring file, skipping")
            continue

        print(f"  Variants in score: {len(scoring_variants)}")

        # Calculate PRS
        prs = calculate_prs(genotype_dict, scoring_variants)
        print(
            f"  Variants matched: {prs['variants_used']}/{prs['variants_total']} "
            f"({prs['overlap_fraction'] * 100:.1f}%)"
        )

        # Check minimum overlap
        if prs["overlap_fraction"] < args.min_overlap:
            print(
                f"  Skipping: overlap {prs['overlap_fraction'] * 100:.1f}% "
                f"below --min-overlap {args.min_overlap * 100:.0f}%"
            )
            continue

        print(f"  Raw PRS: {prs['raw_score']:.6f}")

        # Estimate percentile
        pct_info = estimate_percentile(
            prs["raw_score"], pgs_id, scoring_variants
        )
        if pct_info["percentile"] is not None:
            print(
                f"  Percentile: {pct_info['percentile']:.1f}% "
                f"({pct_info['risk_category']}) "
                f"[{pct_info['method']}]"
            )
        else:
            print("  Percentile: N/A (no reference distribution available)")

        all_results.append({
            "pgs_id": pgs_id,
            "trait": trait,
            "prs": prs,
            "percentile_info": pct_info,
            "metadata": metadata,
            "scoring_variants": scoring_variants,
        })
        print()

    if not all_results:
        print("No scores produced results. Check input file and overlap.")
        sys.exit(0)

    # -----------------------------------------------------------
    # Step 4: Generate report
    # -----------------------------------------------------------
    # Strip per_variant from results before report to keep report manageable
    # (the generate_report function accesses prs["per_variant"] directly)

    report = generate_report(all_results, input_info, args)

    if args.output:
        output_dir = Path(args.output)
        output_dir.mkdir(parents=True, exist_ok=True)

        # Write markdown report
        report_path = output_dir / "prs_report.md"
        report_path.write_text(report)
        print(f"Report written to {report_path}")

        # Write JSON results (without per_variant details for compactness)
        json_results = []
        for r in all_results:
            jr = {
                "pgs_id": r["pgs_id"],
                "trait": r["trait"],
                "raw_score": r["prs"]["raw_score"],
                "variants_used": r["prs"]["variants_used"],
                "variants_total": r["prs"]["variants_total"],
                "overlap_fraction": r["prs"]["overlap_fraction"],
                "percentile": r["percentile_info"]["percentile"],
                "risk_category": r["percentile_info"]["risk_category"],
                "z_score": r["percentile_info"]["z_score"],
                "method": r["percentile_info"]["method"],
                "reference_population": r["percentile_info"].get(
                    "reference_population"
                ),
            }
            json_results.append(jr)

        json_path = output_dir / "prs_results.json"
        json_path.write_text(json.dumps(json_results, indent=2))
        print(f"JSON results written to {json_path}")

        # Write per-variant CSV for detailed analysis
        csv_lines = [
            "pgs_id,rsid,effect_allele,genotype,dosage,weight,contribution,status"
        ]
        for r in all_results:
            for v in r["prs"]["per_variant"]:
                gt = v["genotype"] if v["genotype"] else ""
                dosage = str(v["dosage"]) if v["dosage"] is not None else ""
                csv_lines.append(
                    f"{r['pgs_id']},{v['rsid']},{v['effect_allele']},"
                    f"{gt},{dosage},{v['weight']:.6f},"
                    f"{v['contribution']:.6f},{v['status']}"
                )
        csv_path = output_dir / "prs_variants.csv"
        csv_path.write_text("\n".join(csv_lines) + "\n")
        print(f"Variant details written to {csv_path}")

        print(f"\nFull output in {output_dir}/")
    else:
        # Print report to stdout
        print()
        print(report)

    # Print final summary
    print()
    print("=" * 60)
    print("PRS SUMMARY")
    print("=" * 60)
    for r in all_results:
        pct = r["percentile_info"]
        pct_str = (
            f"{pct['percentile']:.1f}% ({pct['risk_category']})"
            if pct["percentile"] is not None
            else "N/A"
        )
        print(
            f"  {r['pgs_id']:12s} {r['trait']:<30s} "
            f"Score: {r['prs']['raw_score']:.4f}  "
            f"Percentile: {pct_str}"
        )
    print()
    print(DISCLAIMER)


if __name__ == "__main__":
    main()
