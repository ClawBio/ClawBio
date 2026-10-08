#!/usr/bin/env python3
"""
mock_api_server.py — Deterministic mock API server for offline CI testing.

Simulates ClinVar, Ensembl REST, GWAS Catalog, and ClinPGx endpoints with
known test data. Skills hit this instead of live APIs during CI/nightly sweep.

Inspired by StrongDM's simulated Slack/Jira/Okta pattern (Willison, 2026).

Usage:
    # Start server (background)
    python tests/benchmark/mock_api_server.py &

    # Run skills against it
    CLAWBIO_API_BASE_ENSEMBL=http://localhost:8089/ensembl \
    CLAWBIO_API_BASE_CLINPGX=http://localhost:8089/clinpgx \
    python scripts/nightly_demo_sweep.py

    # Stop server
    kill %1

    # Or use as context manager in tests
    from tests.benchmark.mock_api_server import MockAPIServer
    with MockAPIServer(port=8089) as server:
        # run tests against server.base_url
        pass
"""

from __future__ import annotations

import json
import threading
from http.server import HTTPServer, BaseHTTPRequestHandler
from pathlib import Path
from typing import Any
from urllib.parse import urlparse, parse_qs

GROUND_TRUTH_PATH = Path(__file__).parent / "ad_ground_truth.json"

# ---------------------------------------------------------------------------
# Deterministic test data
# ---------------------------------------------------------------------------

# Load ground truth for consistent variant data
_ground_truth: dict = {}
if GROUND_TRUTH_PATH.exists():
    _ground_truth = json.loads(GROUND_TRUTH_PATH.read_text())

_LEAD_VARIANTS = {
    v["rsid"]: v
    for v in _ground_truth.get("lead_variants", {}).get("variants", [])
}


def _ensembl_variation(rsid: str) -> dict:
    """Deterministic Ensembl /variation/human/{rsid} response."""
    known = _LEAD_VARIANTS.get(rsid, {})
    pos = known.get("pos_grch38", 100000)
    chrom = known.get("chr", "1")
    return {
        "name": rsid,
        "var_class": "SNP",
        "most_severe_consequence": "intron_variant",
        "minor_allele": "A",
        "MAF": 0.15,
        "mappings": [
            {
                "assembly_name": "GRCh38",
                "seq_region_name": chrom,
                "start": pos,
                "end": pos,
                "allele_string": "G/A",
            },
            {
                "assembly_name": "GRCh37",
                "seq_region_name": chrom,
                "start": pos - 50000,
                "end": pos - 50000,
                "allele_string": "G/A",
            },
        ],
        "populations": [
            {"population": "1000GENOMES:phase_3:EUR", "allele": "A", "frequency": 0.18},
            {"population": "1000GENOMES:phase_3:AFR", "allele": "A", "frequency": 0.09},
            {"population": "1000GENOMES:phase_3:EAS", "allele": "A", "frequency": 0.12},
        ],
    }


def _ensembl_vep(rsid: str) -> list:
    """Deterministic Ensembl /vep/human/id/{rsid} response."""
    known = _LEAD_VARIANTS.get(rsid, {})
    gene = known.get("gene", "UNKNOWN")
    return [
        {
            "id": rsid,
            "most_severe_consequence": "missense_variant",
            "transcript_consequences": [
                {
                    "gene_symbol": gene,
                    "gene_id": f"ENSG00000{hash(gene) % 100000:05d}",
                    "consequence_terms": ["missense_variant"],
                    "impact": "MODERATE",
                    "biotype": "protein_coding",
                    "sift_prediction": "deleterious",
                    "polyphen_prediction": "probably_damaging",
                },
            ],
        }
    ]


def _gwas_catalog_associations(rsid: str) -> dict:
    """Deterministic GWAS Catalog REST API v2 /v2/associations?rs_id={rsid}."""
    known = _LEAD_VARIANTS.get(rsid, {})
    pvalue = known.get("pvalue", 5e-8)
    gene = known.get("gene", "UNKNOWN")

    # Decompose p-value into mantissa and exponent
    import math
    if pvalue > 0:
        exp = int(math.floor(math.log10(pvalue)))
        mantissa = round(pvalue / (10 ** exp), 2)
    else:
        mantissa, exp = 0, 0

    return {
        "_embedded": {
            "associations": [
                {
                    "association_id": 1,
                    "p_value": pvalue,
                    "pvalue_mantissa": mantissa,
                    "pvalue_exponent": exp,
                    "snp_effect_allele": [f"{rsid}-A"],
                    "snp_allele": [{"rs_id": rsid, "effect_allele": "A"}],
                    "risk_frequency": "0.15",
                    "or_per_copy_num": 1.25,
                    "or_value": "1.25",
                    "beta": "-",
                    "range": "[1.15-1.36]",
                    "efo_traits": [{"efo_id": "MONDO_0004975", "efo_trait": "Alzheimer's disease"}],
                    "reported_trait": ["Alzheimer's disease"],
                    "accession_id": "GCST90027158",
                    "mapped_genes": [gene],
                },
            ]
        },
        "page": {"size": 20, "totalElements": 1, "totalPages": 1, "number": 0},
    }


def _clinpgx_gene(gene: str) -> dict:
    """Deterministic ClinPGx /v1/gene/{gene} response."""
    return {
        "gene": gene,
        "guidelines": [
            {
                "source": "CPIC",
                "guideline_id": f"cpic-{gene.lower()}-001",
                "drugs": ["warfarin", "codeine"] if gene == "CYP2D6" else ["clopidogrel"],
                "recommendation": f"Test {gene} genotype before prescribing",
                "level": "Strong",
            }
        ],
        "phenotypes": [
            {"phenotype": "Poor Metabolizer", "activity_score": "0"},
            {"phenotype": "Normal Metabolizer", "activity_score": "2"},
        ],
        "variants": [
            {"rsid": "rs3892097", "effect": "no function", "allele": "*4"},
        ],
    }


def _clinpgx_drug(drug: str) -> dict:
    """Deterministic ClinPGx /v1/drug/{drug} response."""
    return {
        "drug": drug,
        "genes": ["CYP2D6", "CYP2C19"],
        "fda_label": True,
        "cpic_guideline": True,
        "dpwg_guideline": True,
    }


# ---------------------------------------------------------------------------
# Route table
# ---------------------------------------------------------------------------

ROUTES: dict[str, Any] = {}


def _route(path_parts: list[str], query: dict) -> tuple[int, dict | list]:
    """Match URL path to a handler and return (status_code, response_body)."""

    # Ensembl: /ensembl/variation/human/{rsid}
    if len(path_parts) >= 4 and path_parts[1] == "ensembl":
        if path_parts[2] == "variation" and path_parts[3] == "human":
            rsid = path_parts[4] if len(path_parts) > 4 else "rs000000"
            return 200, _ensembl_variation(rsid)
        if path_parts[2] == "vep" and path_parts[3] == "human":
            # /ensembl/vep/human/id/{rsid}
            rsid = path_parts[5] if len(path_parts) > 5 else "rs000000"
            return 200, _ensembl_vep(rsid)

    # GWAS Catalog v2: /gwas/rest/api/v2/associations?rs_id={rsid}
    if len(path_parts) >= 3 and path_parts[1] == "gwas":
        if "v2" in path_parts and path_parts[-1] == "associations":
            rsid = (query.get("rs_id") or ["rs000000"])[0]
            return 200, _gwas_catalog_associations(rsid)
        # The legacy v1 API is retired upstream and answers 410.
        return 410, {"error": "gone", "message": "This legacy GWAS Catalog REST API is no longer available."}

    # ClinPGx: /clinpgx/v1/gene/{gene} or /clinpgx/v1/drug/{drug}
    if len(path_parts) >= 4 and path_parts[1] == "clinpgx":
        if path_parts[3] == "gene" and len(path_parts) > 4:
            return 200, _clinpgx_gene(path_parts[4])
        if path_parts[3] == "drug" and len(path_parts) > 4:
            return 200, _clinpgx_drug(path_parts[4])

    # Health check
    if path_parts[-1] == "health":
        return 200, {"status": "ok", "server": "clawbio-mock-api", "version": "1.0.0"}

    return 404, {"error": "not_found", "path": "/".join(path_parts)}


# ---------------------------------------------------------------------------
# HTTP Handler
# ---------------------------------------------------------------------------


class MockAPIHandler(BaseHTTPRequestHandler):
    """Handle GET/POST requests with deterministic responses."""

    def do_GET(self):
        parsed = urlparse(self.path)
        path_parts = [p for p in parsed.path.split("/") if p]
        path_parts.insert(0, "")  # Keep leading empty for consistent indexing
        query = parse_qs(parsed.query)

        status, body = _route(path_parts, query)

        self.send_response(status)
        self.send_header("Content-Type", "application/json")
        self.send_header("X-Mock-Server", "clawbio-benchmark")
        self.end_headers()
        self.wfile.write(json.dumps(body, indent=2).encode())

    def do_POST(self):
        # Read body for potential GraphQL (Open Targets)
        content_length = int(self.headers.get("Content-Length", 0))
        body = self.rfile.read(content_length) if content_length else b""

        parsed = urlparse(self.path)
        path_parts = [p for p in parsed.path.split("/") if p]
        path_parts.insert(0, "")
        query = parse_qs(parsed.query)

        status, resp = _route(path_parts, query)

        self.send_response(status)
        self.send_header("Content-Type", "application/json")
        self.send_header("X-Mock-Server", "clawbio-benchmark")
        self.end_headers()
        self.wfile.write(json.dumps(resp, indent=2).encode())

    def log_message(self, format, *args):
        """Suppress request logging in CI."""
        pass


# ---------------------------------------------------------------------------
# Server wrapper
# ---------------------------------------------------------------------------


class MockAPIServer:
    """Context manager for starting/stopping the mock API server."""

    def __init__(self, host: str = "127.0.0.1", port: int = 8089):
        self.host = host
        self.port = port
        self.server: HTTPServer | None = None
        self.thread: threading.Thread | None = None

    @property
    def base_url(self) -> str:
        return f"http://{self.host}:{self.port}"

    def start(self):
        self.server = HTTPServer((self.host, self.port), MockAPIHandler)
        self.thread = threading.Thread(target=self.server.serve_forever, daemon=True)
        self.thread.start()

    def stop(self):
        if self.server:
            self.server.shutdown()
            self.server = None

    def __enter__(self):
        self.start()
        return self

    def __exit__(self, *args):
        self.stop()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main():
    import argparse

    parser = argparse.ArgumentParser(description="ClawBio mock API server")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8089)
    args = parser.parse_args()

    print(f"Mock API server starting on http://{args.host}:{args.port}")
    print(f"  Ensembl:      http://{args.host}:{args.port}/ensembl/variation/human/rs6733839")
    print(f"  GWAS Catalog:  http://{args.host}:{args.port}/gwas/rest/api/v2/associations?rs_id=rs6733839")
    print(f"  ClinPGx:       http://{args.host}:{args.port}/clinpgx/v1/gene/CYP2D6")
    print(f"  Health:        http://{args.host}:{args.port}/health")
    print("Press Ctrl+C to stop.")

    server = HTTPServer((args.host, args.port), MockAPIHandler)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nShutting down.")
        server.shutdown()


if __name__ == "__main__":
    main()
