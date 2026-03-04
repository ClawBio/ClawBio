"""
test_html_report.py — Tests for the rich-text HTML report generator.

Run with: pytest skills/pharmgx-reporter/tests/test_html_report.py -v
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from pharmgx_reporter import (
    PGX_SNPS,
    GENE_DEFS,
    parse_file,
    call_diplotype,
    call_phenotype,
    lookup_drugs,
)
from html_report import generate_html_report

DEMO = Path(__file__).parent.parent / "demo_patient.txt"


def _build_report_data():
    """Parse demo data and return all args needed for generate_html_report."""
    fmt, total_snps, pgx_snps = parse_file(str(DEMO))
    profiles = {}
    for gene in GENE_DEFS:
        diplotype = call_diplotype(gene, pgx_snps)
        phenotype = call_phenotype(gene, diplotype)
        profiles[gene] = {"diplotype": diplotype, "phenotype": phenotype}
    drug_results = lookup_drugs(profiles)
    return fmt, total_snps, pgx_snps, profiles, drug_results


# ── Structure ────────────────────────────────────────────────────────────────

def test_html_report_is_valid_html():
    fmt, total_snps, pgx_snps, profiles, drug_results = _build_report_data()
    html = generate_html_report(
        str(DEMO), fmt, total_snps, pgx_snps, profiles,
        drug_results, GENE_DEFS, PGX_SNPS,
    )
    assert html.startswith("<!DOCTYPE html>")
    assert "</html>" in html
    assert "<title>" in html


def test_html_report_has_research_banner():
    fmt, total_snps, pgx_snps, profiles, drug_results = _build_report_data()
    html = generate_html_report(
        str(DEMO), fmt, total_snps, pgx_snps, profiles,
        drug_results, GENE_DEFS, PGX_SNPS,
    )
    assert "DEMO" in html
    assert "FOR RESEARCH USE ONLY" in html


def test_html_report_contains_all_sections():
    fmt, total_snps, pgx_snps, profiles, drug_results = _build_report_data()
    html = generate_html_report(
        str(DEMO), fmt, total_snps, pgx_snps, profiles,
        drug_results, GENE_DEFS, PGX_SNPS,
    )
    for section_id in ["summary", "alerts", "genes", "variants", "drugs",
                       "methods", "disclaimer"]:
        assert f'id="{section_id}"' in html, f"Missing section: {section_id}"


def test_html_report_has_sidebar_nav():
    fmt, total_snps, pgx_snps, profiles, drug_results = _build_report_data()
    html = generate_html_report(
        str(DEMO), fmt, total_snps, pgx_snps, profiles,
        drug_results, GENE_DEFS, PGX_SNPS,
    )
    assert '<nav class="sidebar">' in html
    assert 'href="#summary"' in html
    assert 'href="#alerts"' in html
    assert 'href="#drugs"' in html


# ── Content ──────────────────────────────────────────────────────────────────

def test_html_report_contains_avoid_drugs():
    fmt, total_snps, pgx_snps, profiles, drug_results = _build_report_data()
    html = generate_html_report(
        str(DEMO), fmt, total_snps, pgx_snps, profiles,
        drug_results, GENE_DEFS, PGX_SNPS,
    )
    assert "Codeine" in html
    assert "Tramadol" in html
    assert "badge-avoid" in html


def test_html_report_contains_all_gene_profiles():
    fmt, total_snps, pgx_snps, profiles, drug_results = _build_report_data()
    html = generate_html_report(
        str(DEMO), fmt, total_snps, pgx_snps, profiles,
        drug_results, GENE_DEFS, PGX_SNPS,
    )
    for gene in GENE_DEFS:
        assert gene in html, f"Missing gene: {gene}"


def test_html_report_contains_disclaimer():
    fmt, total_snps, pgx_snps, profiles, drug_results = _build_report_data()
    html = generate_html_report(
        str(DEMO), fmt, total_snps, pgx_snps, profiles,
        drug_results, GENE_DEFS, PGX_SNPS,
    )
    assert "research and educational purposes only" in html
    assert "cpicpgx.org" in html


def test_html_report_has_colour_coded_badges():
    fmt, total_snps, pgx_snps, profiles, drug_results = _build_report_data()
    html = generate_html_report(
        str(DEMO), fmt, total_snps, pgx_snps, profiles,
        drug_results, GENE_DEFS, PGX_SNPS,
    )
    assert "badge-avoid" in html
    assert "badge-caution" in html
    assert "badge-ok" in html


def test_html_report_is_self_contained():
    """No external CSS/JS references — everything is inline."""
    fmt, total_snps, pgx_snps, profiles, drug_results = _build_report_data()
    html = generate_html_report(
        str(DEMO), fmt, total_snps, pgx_snps, profiles,
        drug_results, GENE_DEFS, PGX_SNPS,
    )
    assert '<link rel="stylesheet"' not in html
    assert "<script src=" not in html
    assert "<style>" in html
