#!/usr/bin/env python3
"""
GWAS Polygenic Risk Score Calculator — PGS Catalog integration.

Calculates Polygenic Risk Scores (PRS) from DTC genetic data (23andMe/AncestryDNA)
using the PGS Catalog REST API and pre-curated scoring files.

Usage:
    python gwas_prs.py --input <23andme_file> --trait "type 2 diabetes" --output <dir>
    python gwas_prs.py --input <23andme_file> --pgs-id PGS000031 --output <dir>
    python gwas_prs.py --input <23andme_file> --panel-id CLAWBIO-T2D-8 --output <dir>
    python gwas_prs.py --demo --output /tmp/prs_demo
"""

from __future__ import annotations

import argparse
import os
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
# Shared library imports
# ---------------------------------------------------------------------------

_SKILL_DIR = Path(__file__).resolve().parent
_PROJECT_ROOT = _SKILL_DIR.parent.parent
for _import_path in (_PROJECT_ROOT, _SKILL_DIR):
    if str(_import_path) not in sys.path:
        sys.path.insert(0, str(_import_path))

from prs_evidence import assess as assess_evidence  # noqa: E402
from repro_bundle import create_reproducibility_bundle  # noqa: E402

from clawbio.common.checksums import sha256_hex  # noqa: E402
from clawbio.common.parsers import genotypes_to_simple, parse_genetic_file  # noqa: E402
from clawbio.common.report import DISCLAIMER as _SHARED_DISCLAIMER  # noqa: E402
from clawbio.common.report import write_result_json  # noqa: E402

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

PGS_API_BASE = "https://www.pgscatalog.org/rest"
RATE_LIMIT_INTERVAL = 0.55  # seconds between requests (stay under 2 req/sec)
CACHE_TTL = 86400  # 24 hours
USER_AGENT = "ClawBio-GWAS-PRS/0.3.0"

DISCLAIMER = _SHARED_DISCLAIMER

# Risk category thresholds (percentile-based)
RISK_CATEGORIES = [
    (25, "Low"),
    (75, "Average"),
    (95, "Elevated"),
    (100, "High"),
]

SKILL_DIR = _SKILL_DIR
DATA_DIR = SKILL_DIR / "data"

# ---------------------------------------------------------------------------
# Curated demo panels: offline scoring with reference distributions
# ---------------------------------------------------------------------------
#
# These six entries are ClawBio-curated illustrative panels of well-established
# trait-associated loci. They are NOT PGS Catalog scoring files. The
# ``CLAWBIO-*`` id is the dictionary key and canonical identity. Historical PGS
# accessions remain metadata so old output can be recognised without presenting
# the panel as that Catalog score (issue #356).
#
# `pgs_catalog_id` is per panel, not a blanket None. For five of the six the
# accession belongs to an entirely DIFFERENT published score, so it is None.
# PGS000001 is the exception: that accession really is Mavaddat 2015 PRS77_BC,
# the same trait and the same variant count as this panel, so recording None
# there would discard a true mapping. The bundled file is still not the
# published score; it shares 60 of 77 rsIDs and 31 of those 60 weights differ
# by more than 0.02, so treat it as a lossy derivative, not a copy.
#
# `publication` is the LOCUS reference, not the source of the weights. The
# weights here are approximate and are not the published betas: Vassy 2014 is a
# 62-locus score against this 8-locus panel, and Abraham 2016 is 49,310 SNPs
# against 46. clawbio-bench's own ground truth attributes the T2D betas to
# Morris 2012 / Mahajan 2018. Do not let a reader infer the weights came from
# the cited paper.
#
# Every `pmid` below was resolved against PubMed on 2026-08-21. Three earlier
# citations named papers that had nothing to do with the trait: a
# ventilator-pneumonia nursing paper, the ExAC paper, and a body-fat
# distribution meta-analysis. Do not edit a citation without re-resolving it;
# `tests/test_score_provenance.py` pins them.

CURATED_SCORES: dict[str, dict] = {
    "CLAWBIO-T2D-8": {
        "name": "Curated demo panel: type 2 diabetes (8 loci)",
        "curated_panel_id": "CLAWBIO-T2D-8",
        "legacy_pgs_id": "PGS000013",
        "pgs_catalog_id": None,
        "trait_id": "EFO_0001360",
        "trait": "Type 2 diabetes",
        "variants_count": 8,
        "publication": "Vassy JL et al. (2014) Diabetes",
        "pmid": "24520119",
        "reference_distribution": {
            "mean": 1.12,
            "sd": 0.30,
            "population": "EUR",
        },
    },
    "CLAWBIO-AF-12": {
        "name": "Curated demo panel: atrial fibrillation (12 loci)",
        "curated_panel_id": "CLAWBIO-AF-12",
        "legacy_pgs_id": "PGS000011",
        "pgs_catalog_id": None,
        "trait_id": "EFO_0000275",
        "trait": "Atrial fibrillation",
        "variants_count": 12,
        "publication": "Tada H et al. (2014) Stroke",
        "pmid": "25123217",
        "reference_distribution": {
            "mean": 0.65,
            "sd": 0.23,
            "population": "EUR",
        },
    },
    "CLAWBIO-CAD-46": {
        "name": "Curated demo panel: coronary artery disease (46 loci)",
        "curated_panel_id": "CLAWBIO-CAD-46",
        "legacy_pgs_id": "PGS000004",
        "pgs_catalog_id": None,
        "trait_id": "EFO_0001645",
        "trait": "Coronary artery disease",
        "variants_count": 46,
        "publication": "Abraham G et al. (2016) Eur Heart J",
        "pmid": "27655226",
        "reference_distribution": {
            "mean": 2.84,
            "sd": 0.28,
            "population": "EUR",
        },
    },
    "CLAWBIO-BC-77": {
        "name": "Curated demo panel: breast cancer (77 loci)",
        "curated_panel_id": "CLAWBIO-BC-77",
        "legacy_pgs_id": "PGS000001",
        "pgs_catalog_id": "PGS000001",
        "trait_id": "EFO_0000305",
        "trait": "Breast cancer",
        "variants_count": 77,
        "publication": "Mavaddat N et al. (2015) J Natl Cancer Inst",
        "pmid": "25855707",
        "reference_distribution": {
            "mean": 4.23,
            "sd": 0.54,
            "population": "EUR",
        },
    },
    "CLAWBIO-PC-147": {
        "name": "Curated demo panel: prostate cancer (147 loci)",
        "curated_panel_id": "CLAWBIO-PC-147",
        "legacy_pgs_id": "PGS000057",
        "pgs_catalog_id": None,
        "trait_id": "EFO_0001663",
        "trait": "Prostate cancer",
        "variants_count": 147,
        "publication": "Schumacher FR et al. (2018) Nat Genet",
        "pmid": "29892016",
        "reference_distribution": {
            "mean": 7.11,
            "sd": 0.56,
            "population": "EUR",
        },
    },
    "CLAWBIO-BMI-97": {
        "name": "Curated demo panel: body mass index (97 loci)",
        "curated_panel_id": "CLAWBIO-BMI-97",
        "legacy_pgs_id": "PGS000039",
        "pgs_catalog_id": None,
        "trait_id": "EFO_0004340",
        "trait": "BMI",
        "variants_count": 97,
        "publication": "Locke AE et al. (2015) Nature",
        "pmid": "25673413",
        "reference_distribution": {
            "mean": 2.89,
            "sd": 0.25,
            "population": "EUR",
        },
    },
}

# The pinned clawbio_bench revision in CI still invokes the T2D panel as
# ``--pgs-id PGS000013`` and searches that field in ``prs_results.json``.
# Keep exactly this one explicit compatibility alias until the benchmark can
# request ``CLAWBIO-T2D-8``. The other five historical accessions follow the
# real PGS Catalog path.
LEGACY_PGS_PANEL_COMPAT = {
    "PGS000013": "CLAWBIO-T2D-8",
}

# What the aliased accession really is, so a refusal can say so.
LEGACY_PGS_ACCESSION_FACTS = {
    "PGS000013": "Khera AV et al. 2018, coronary artery disease, 6,630,150 variants",
}

# Issue #356, remaining half. The alias above fires only when the caller opts
# in by setting this variable to a true value. The benchmark workflow sets it
# (.github/workflows/bench-leaderboard.yml); nothing else should. Without it,
# ``--pgs-id PGS000013`` is treated as what it is, a PGS Catalog accession,
# and fetched from the Catalog like any other. Scoring the wrong panel in
# answer to a specific request is never the default.
LEGACY_ALIAS_ENV = "CLAWBIO_ALLOW_LEGACY_PGS_ALIAS"
_TRUE_VALUES = {"1", "true", "yes", "on"}


def legacy_alias_enabled(environ=None) -> bool:
    """True only when ``CLAWBIO_ALLOW_LEGACY_PGS_ALIAS`` holds a true value."""
    env = os.environ if environ is None else environ
    return str(env.get(LEGACY_ALIAS_ENV, "")).strip().lower() in _TRUE_VALUES


def explain_legacy_accession(pgs_id: str, panel_id: str) -> None:
    """Tell the caller what the accession is and how to get the panel instead."""
    facts = LEGACY_PGS_ACCESSION_FACTS.get(pgs_id, "a published PGS Catalog score")
    print(
        f"GWAS-PRS: {pgs_id} is a PGS Catalog accession ({facts}), "
        f"not the curated panel {panel_id}."
    )
    print("  ClawBio does not substitute the panel for the accession.")
    print(f"  For the curated panel use --panel-id {panel_id}.")
    print("  See ClawBio issue #356.")
    print()


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
# Genotype parser — delegates to clawbio.common.parsers
# ---------------------------------------------------------------------------


def detect_format(lines: list[str]) -> str:
    """Detect file format from a list of header lines.

    Thin wrapper around clawbio.common.parsers.detect_format for backward
    compatibility with tests that pass lines rather than a file path.
    """
    from clawbio.common.parsers import detect_format as _detect_fmt
    import tempfile, os
    # Write lines to a temp file so the shared detector can inspect them
    fd, tmp = tempfile.mkstemp(suffix=".txt")
    try:
        with os.fdopen(fd, "w") as f:
            f.write("\n".join(lines))
        return _detect_fmt(Path(tmp))
    finally:
        os.unlink(tmp)


def load_genotypes(path: str | Path) -> tuple[str, int, dict[str, str]]:
    """Parse a DTC genotype file via the shared parser.

    Supports .gz files, 23andMe (tab-separated) and AncestryDNA
    (comma-separated) formats.

    Returns:
        (format_name, total_snps, {rsid: genotype_string})
    """
    from clawbio.common.parsers import detect_format as _detect_fmt

    path = Path(path)
    fmt = _detect_fmt(path)
    records = parse_genetic_file(str(path))
    genotypes = genotypes_to_simple(records)
    return fmt, len(genotypes), genotypes


# Backward-compatible alias
parse_genotype_file = load_genotypes


# ---------------------------------------------------------------------------
# Scoring file parser — PGS Catalog harmonized TSV
# ---------------------------------------------------------------------------

CURATED_PANEL_MARKER = "#clawbio_panel=curated_demo"


def curated_panel_path(panel_id: str, build: str = "GRCh37") -> Path:
    """Return the bundled file for a canonical ``CLAWBIO-*`` panel id."""
    normalized = panel_id.strip().upper()
    if normalized not in CURATED_SCORES:
        raise KeyError(normalized)
    return DATA_DIR / f"{normalized}_{build}.txt"


def is_curated_demo_panel(filepath: str | Path) -> bool:
    """True if `filepath` is a bundled ClawBio demo panel, not a catalog score.

    The demo panels live at the same path the download cache uses, so without
    this check a request for a real accession is answered from the local
    lookalike and never reaches the API. That is issue #356: `--pgs-id
    PGS000013` returned an 8-variant type 2 diabetes panel in place of the
    real PGS000013 (Khera 2018, coronary artery disease, 6,630,150 variants).

    Reads only the leading comment block, so it is cheap on large real files
    and safe on a missing or unreadable path.
    """
    path = Path(filepath)
    try:
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "rt") as handle:  # type: ignore[operator]
            for line in handle:
                if not line.startswith("#"):
                    return False
                if line.strip() == CURATED_PANEL_MARKER:
                    return True
    except (OSError, EOFError, UnicodeDecodeError):
        return False
    return False


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
    score_id: str,
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
    curated = CURATED_SCORES.get(score_id)
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
            score_id, pgs_id, trait, prs (from calculate_prs), percentile_info
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
        f"**Tool**: ClawBio GWAS-PRS v0.2.0",
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
        "| Score ID | Trait | Raw Score | Percentile | Risk Category | Overlap |"
    )
    lines.append(
        "|--------|-------|-----------|------------|---------------|---------|"
    )

    for r in results:
        score_id = r["score_id"]
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
            f"| {score_id} | {trait} | {raw} | {pct} | {risk} | {overlap} |"
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
                f"- **{r['score_id']}** ({r.get('trait', '')}): "
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
        score_id = r["score_id"]
        pgs_id = r.get("pgs_id")
        trait = r.get("trait", "Unknown")
        prs = r["prs"]
        pct_info = r["percentile_info"]
        metadata = r.get("metadata", {})

        lines.append(f"### {score_id} — {trait}")
        lines.append("")

        # Metadata table
        lines.append("| Property | Value |")
        lines.append("|----------|-------|")
        lines.append(f"| **Score ID** | {score_id} |")
        if pgs_id and not r.get("curated_demo_panel"):
            lines.append(f"| **PGS Catalog ID** | {pgs_id} |")
        lines.append(f"| **Trait** | {trait} |")
        if r.get("curated_demo_panel"):
            # Issue #356. Stated here as well as in the provenance footer,
            # because a reader who skims takes the table and stops.
            lines.append(
                f"| **Source** | ClawBio curated demo panel "
                f"`{r.get('curated_panel_id') or score_id}`, "
                "NOT a PGS Catalog score |"
            )
            if r.get("legacy_pgs_id"):
                lines.append(
                    f"| **Historical alias** | `{r['legacy_pgs_id']}` "
                    "(compatibility only; not Catalog identity) |"
                )
        if metadata.get("publication"):
            label = ("Loci reference" if r.get("curated_demo_panel")
                     else "Publication")
            lines.append(f"| **{label}** | {metadata['publication']} |")
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
        assessment = r.get("evidence_assessment")
        if assessment:
            lines.append(f"| **Evidence assessment** | {assessment['status']} |")
            for reason in assessment["reasons"]:
                lines.append(f"| {reason['code']} | {reason['action']} |")
            lines.append("| **Interpretation** | Research percentile only; no individual disease-risk or causal claim |")

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
            lines.append("| **Percentile** | Withheld (see evidence assessment or reference availability) |")
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
    # Issue #356: do not assert PGS Catalog provenance for a curated panel.
    # The warning used to reach stdout only, which is the one channel no
    # downstream consumer, and no reader of this file, ever sees.
    if any(r.get("curated_demo_panel") for r in results):
        lines.append(
            "**Scoring files**: ClawBio curated demo panel, NOT a PGS Catalog "
            "score. The loci are established trait associations but the weights "
            "are approximate and are not the published betas. Curated panels "
            "use `CLAWBIO-*` identifiers and filenames; historical PGS accessions "
            "are provenance metadata only. Do not cite these results as a "
            "published PRS. See ClawBio issue #356."
        )
    else:
        lines.append("**Scoring files**: PGS Catalog "
                     "(https://www.pgscatalog.org/)")
    lines.append("")
    lines.append("**Genome build**: " + args.build)
    lines.append("")

    if any(r.get("evidence_assessment") for r in results):
        lines.append(
            "**Evidence policy 1.0.0**: Research percentiles require caller-supplied, "
            "hash-bound evidence for the exact score and input, complete compatible "
            "variants, independent validation and a normal reference distribution "
            "for the target population/context. Earlier percentile methods above "
            "apply only to the synthetic demo. Evidence is not independently "
            "verified; a supported percentile does not establish clinical utility "
            "or causality. Unsupported raw sums are retained for audit only."
        )
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
        help="Specific PGS Catalog score ID (e.g., PGS000031)",
    )
    parser.add_argument(
        "--panel-id",
        help="Bundled ClawBio curated panel ID (e.g., CLAWBIO-T2D-8)",
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

    parser.add_argument(
        "--evidence-json",
        help="Local hash-bound score/input evidence manifest for research percentiles",
    )
    args = parser.parse_args()
    evidence = None
    if args.evidence_json:
        try:
            evidence = json.loads(Path(args.evidence_json).read_text())
            if (not isinstance(evidence, dict) or evidence.get("schema_version") != 1
                    or not isinstance(evidence.get("input"), dict)
                    or not isinstance(evidence.get("scores"), dict)):
                raise ValueError("expected schema_version 1 with input and scores objects")
        except (OSError, ValueError) as exc:
            parser.error(f"Invalid --evidence-json: {exc}")
    evidence_gate = not args.demo or bool(args.evidence_json)
    if args.output and any((Path(args.output) / name).exists() for name in (
            "prs_report.md", "prs_results.json", "prs_variants.csv", "result.json", "reproducibility")):
        parser.error("Output artifacts already exist; choose a new directory.")


    selectors = [args.demo, args.trait, args.pgs_id, args.panel_id]
    if sum(bool(value) for value in selectors) > 1:
        print(
            "Error: choose exactly one of --demo, --trait, --pgs-id, "
            "or --panel-id"
        )
        sys.exit(1)

    if not any(selectors):
        parser.print_help()
        print(
            "\nError: provide at least one of --demo, --trait, --pgs-id, "
            "or --panel-id"
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
    # Each entry: {score_id, pgs_id, trait, filepath, metadata}

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

        for panel_id, meta in CURATED_SCORES.items():
            fpath = curated_panel_path(panel_id, args.build)
            if not fpath.exists():
                print(f"  WARNING: Scoring file not found for {panel_id}, skipping")
                continue
            scoring_files.append({
                "score_id": panel_id,
                "pgs_id": None,
                "legacy_pgs_id": meta["legacy_pgs_id"],
                "legacy_pgs_compatibility": False,
                "trait": meta["trait"],
                "filepath": fpath,
                "metadata": {
                    "publication": meta.get("publication", ""),
                    "variants_count": meta.get("variants_count", 0),
                },
            })
            print(f"  Loaded curated panel: {panel_id} ({meta['trait']})")

        print()

    elif args.panel_id:
        panel_id = args.panel_id.strip().upper()
        meta = CURATED_SCORES.get(panel_id)
        if meta is None:
            available = ", ".join(sorted(CURATED_SCORES))
            print(
                f"Error: unknown curated panel ID {panel_id!r}. "
                f"Available panels: {available}"
            )
            sys.exit(1)

        filepath = curated_panel_path(panel_id, args.build)
        if not filepath.exists():
            print(f"Error: curated panel file not found: {filepath}")
            sys.exit(1)

        print(f"GWAS-PRS: using curated panel {panel_id}")
        print("  This is an illustrative ClawBio panel, not a PGS Catalog score.")
        print()
        scoring_files.append({
            "score_id": panel_id,
            "pgs_id": None,
            "legacy_pgs_id": meta["legacy_pgs_id"],
            "legacy_pgs_compatibility": False,
            "trait": meta["trait"],
            "filepath": filepath,
            "metadata": {
                "publication": meta.get("publication", ""),
                "variants_count": meta.get("variants_count", 0),
            },
        })

    elif args.pgs_id:
        # Specific PGS ID mode
        pgs_id = args.pgs_id.strip().upper()
        if not pgs_id.startswith("PGS"):
            pgs_id = "PGS" + pgs_id.lstrip("0")

        compat_panel_id = LEGACY_PGS_PANEL_COMPAT.get(pgs_id)
        if compat_panel_id and legacy_alias_enabled():
            meta = CURATED_SCORES[compat_panel_id]
            filepath = curated_panel_path(compat_panel_id, args.build)
            if not filepath.exists():
                print(f"Error: curated panel file not found: {filepath}")
                sys.exit(1)
            print(
                f"GWAS-PRS: benchmark compatibility maps {pgs_id} to "
                f"curated panel {compat_panel_id}"
            )
            print(
                "  WARNING: this is not the PGS Catalog score. Use "
                f"--panel-id {compat_panel_id} for normal curated-panel runs."
            )
            print("  See ClawBio issue #356.")
            print()
            scoring_files.append({
                "score_id": compat_panel_id,
                "pgs_id": pgs_id,
                "legacy_pgs_id": pgs_id,
                "legacy_pgs_compatibility": True,
                "trait": meta["trait"],
                "filepath": filepath,
                "metadata": {
                    "publication": meta.get("publication", ""),
                    "variants_count": meta.get("variants_count", 0),
                },
            })
        else:
            if compat_panel_id:
                explain_legacy_accession(pgs_id, compat_panel_id)
            print(f"GWAS-PRS: fetching score {pgs_id}")
            print()
            local_path = DATA_DIR / f"{pgs_id}_hmPOS_{args.build}.txt"
            local_path_gz = DATA_DIR / f"{pgs_id}_hmPOS_{args.build}.txt.gz"
            for candidate in (local_path, local_path_gz):
                if candidate.exists() and is_curated_demo_panel(candidate):
                    print(
                        f"Error: refusing curated panel at PGS Catalog cache "
                        f"path {candidate}. Use --panel-id instead."
                    )
                    sys.exit(1)

            meta_data = None
            if local_path.exists():
                filepath = local_path
                print(f"  Using pre-downloaded scoring file: {filepath}")
            elif local_path_gz.exists():
                filepath = local_path_gz
                print(f"  Using pre-downloaded scoring file: {filepath}")
            else:
                try:
                    print(f"  Fetching metadata for {pgs_id}...")
                    meta_data = client.get_score_metadata(pgs_id)
                except requests.RequestException as e:
                    print(f"Error: could not fetch {pgs_id}: {e}")
                    if compat_panel_id:
                        print(f"  For the curated panel use --panel-id {compat_panel_id}.")
                    sys.exit(1)

                trait_name = "Unknown"
                trait_efo = meta_data.get("trait_efo", [])
                if trait_efo:
                    trait_name = trait_efo[0].get("label", "Unknown")
                elif meta_data.get("trait_reported"):
                    trait_name = ", ".join(meta_data["trait_reported"])
                variant_count = meta_data.get("variants_number", 0)
                print(f"  Trait: {trait_name}")
                print(f"  Variants: {variant_count}")
                if variant_count > args.max_variants:
                    print(
                        f"  Skipping: {variant_count} variants exceeds "
                        f"--max-variants {args.max_variants}"
                    )
                    sys.exit(0)
                try:
                    filepath = client.download_scoring_file(
                        pgs_id, build=args.build
                    )
                except requests.HTTPError as e:
                    print(f"Error downloading scoring file for {pgs_id}: {e}")
                    sys.exit(1)

            trait = "Unknown"
            publication = ""
            try:
                meta_data = meta_data or client.get_score_metadata(pgs_id)
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
                "score_id": pgs_id,
                "pgs_id": pgs_id,
                "legacy_pgs_id": None,
                "legacy_pgs_compatibility": False,
                "trait": trait,
                "filepath": filepath,
                "metadata": {"publication": publication},
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

            invalid_cache_entry = False
            for candidate in (local_path, local_path_gz):
                if candidate.exists() and is_curated_demo_panel(candidate):
                    print(
                        f"  WARNING: refusing curated panel at PGS Catalog "
                        f"cache path {candidate}; use --panel-id instead."
                    )
                    invalid_cache_entry = True
            if invalid_cache_entry:
                continue

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
                "score_id": sid,
                "pgs_id": sid,
                "legacy_pgs_id": None,
                "legacy_pgs_compatibility": False,
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
    fmt, total_snps, genotype_dict = load_genotypes(input_path)
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
        score_id = sf["score_id"]
        pgs_id = sf.get("pgs_id")
        trait = sf["trait"]
        filepath = sf["filepath"]
        metadata = sf["metadata"]

        print(f"Scoring {score_id} ({trait})...")

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
        if prs["overlap_fraction"] < args.min_overlap and not evidence_gate:
            print(
                f"  Skipping: overlap {prs['overlap_fraction'] * 100:.1f}% "
                f"below --min-overlap {args.min_overlap * 100:.0f}%"
            )
            continue

        print(f"  Raw PRS: {prs['raw_score']:.6f}")

        # Real inputs always require explicit evidence. Only the built-in
        # synthetic demo retains illustrative legacy percentile estimates.
        evidence_assessment = None
        if evidence_gate:
            evidence_assessment = assess_evidence(
                score_id=score_id, scoring_path=filepath, input_path=input_path,
                build=args.build, variants=scoring_variants,
                genotypes=genotype_dict, raw_score=prs["raw_score"],
                evidence=evidence, curated=is_curated_demo_panel(filepath),
            )
            pct_info = {key: evidence_assessment[key] for key in (
                "percentile", "z_score", "risk_category", "method", "reference_population"
            )}
            print("  Evidence assessment: " + evidence_assessment["status"])
            for reason in evidence_assessment["reasons"]:
                print("    " + reason["code"] + ": " + reason["action"])
        else:
            pct_info = estimate_percentile(
                prs["raw_score"], score_id, scoring_variants
            )
        if pct_info["percentile"] is not None:
            print(
                f"  Percentile: {pct_info['percentile']:.1f}% "
                f"({pct_info['risk_category']}) "
                f"[{pct_info['method']}]"
            )
        else:
            print("  Percentile: N/A (no reference distribution available)")

        # Issue #356: carry the panel marker on the result itself, so every
        # durable artefact can state it. Derived from the file that was
        # actually scored, not from the accession, so it stays true whichever
        # path (--demo, --pgs-id, --trait) selected it.
        curated = is_curated_demo_panel(sf["filepath"])
        curated_meta = CURATED_SCORES.get(score_id, {}) if curated else {}
        all_results.append({
            "score_id": score_id,
            "pgs_id": pgs_id,
            "trait": trait,
            "prs": prs,
            "percentile_info": pct_info,
            "evidence_assessment": evidence_assessment,
            "interpretation_scope": "research_percentile" if evidence_gate else "synthetic_demo_only",
            "metadata": metadata,
            "scoring_variants": scoring_variants,
            "scoring_file_path": Path(filepath),
            "curated_demo_panel": curated,
            "curated_panel_id": score_id if curated else None,
            "legacy_pgs_id": curated_meta.get("legacy_pgs_id"),
            "legacy_pgs_compatibility": bool(
                sf.get("legacy_pgs_compatibility")
            ),
            "pgs_catalog_id": (
                curated_meta.get("pgs_catalog_id") if curated else pgs_id
            ),
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
                "score_id": r["score_id"],
                "pgs_id": r["pgs_id"],
                "curated_demo_panel": bool(r["curated_demo_panel"]),
                "curated_panel_id": r["curated_panel_id"],
                "legacy_pgs_id": r["legacy_pgs_id"],
                "legacy_pgs_compatibility": r["legacy_pgs_compatibility"],
                "pgs_catalog_id": r["pgs_catalog_id"],
                "trait": r["trait"],
                "raw_score": r["prs"]["raw_score"],
                "variants_used": r["prs"]["variants_used"],
                "variants_total": r["prs"]["variants_total"],
                "overlap_fraction": r["prs"]["overlap_fraction"],
                "evidence_assessment": r.get("evidence_assessment"),
                "interpretation_scope": r.get("interpretation_scope"),
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
            "score_id,rsid,effect_allele,genotype,dosage,weight,contribution,status"
        ]
        for r in all_results:
            for v in r["prs"]["per_variant"]:
                gt = v["genotype"] if v["genotype"] else ""
                dosage = str(v["dosage"]) if v["dosage"] is not None else ""
                csv_lines.append(
                    f"{r['score_id']},{v['rsid']},{v['effect_allele']},"
                    f"{gt},{dosage},{v['weight']:.6f},"
                    f"{v['contribution']:.6f},{v['status']}"
                )
        csv_path = output_dir / "prs_variants.csv"
        csv_path.write_text("\n".join(csv_lines) + "\n")
        print(f"Variant details written to {csv_path}")

        # Write standardized result.json envelope
        first = all_results[0] if all_results else {}
        first_pct = first.get("percentile_info", {})
        result_json_path = write_result_json(
            output_dir=output_dir,
            skill="gwas-prs",
            version="0.3.0",
            summary={
                "scores_calculated": len(all_results),
                "trait": first.get("trait", ""),
                "score_id": first.get("score_id"),
                "pgs_id": first.get("pgs_id"),
                # Issue #356: a consumer reading this envelope, for example
                # profile-report, must be able to tell a curated panel from a
                # PGS Catalog score without reading our stdout.
                "curated_demo_panel": bool(first.get("curated_demo_panel")),
                "curated_panel_id": first.get("curated_panel_id"),
                "legacy_pgs_id": first.get("legacy_pgs_id"),
                "legacy_pgs_compatibility": first.get(
                    "legacy_pgs_compatibility", False
                ),
                "pgs_catalog_id": first.get("pgs_catalog_id"),
                "raw_score": first.get("prs", {}).get("raw_score"),
                "evidence_assessment": first.get("evidence_assessment"),
                "interpretation_scope": first.get("interpretation_scope"),
                "percentile": first_pct.get("percentile"),
                "risk_category": first_pct.get("risk_category"),
                "overlap_fraction": first.get("prs", {}).get("overlap_fraction"),
            },
            data={
                "input_info": input_info,
                "results": json_results,
            },
            input_checksum=sha256_hex(str(input_path)) if input_path.exists() else "",
        )
        print(f"Result envelope written to {result_json_path}")

        reproducibility_paths = create_reproducibility_bundle(
            output_dir=output_dir,
            input_path=input_path,
            input_info=input_info,
            scoring_files=[
                {
                    "score_id": result["score_id"],
                    "pgs_id": result["pgs_id"],
                    "trait": result["trait"],
                    "filepath": result["scoring_file_path"],
                    "curated_demo_panel": result["curated_demo_panel"],
                    "curated_panel_id": result["curated_panel_id"],
                    "legacy_pgs_id": result["legacy_pgs_id"],
                    "legacy_pgs_compatibility": result["legacy_pgs_compatibility"],
                    "pgs_catalog_id": result["pgs_catalog_id"],
                }
                for result in all_results
            ],
            args=args,
            output_paths=[report_path, json_path, csv_path, result_json_path],
        )
        print(
            "Reproducibility bundle written to "
            f"{reproducibility_paths['commands'].parent}"
        )

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
            f"  {r['score_id']:16s} {r['trait']:<30s} "
            f"Score: {r['prs']['raw_score']:.4f}  "
            f"Percentile: {pct_str}"
        )
    print()
    print(DISCLAIMER)


if __name__ == "__main__":
    main()
