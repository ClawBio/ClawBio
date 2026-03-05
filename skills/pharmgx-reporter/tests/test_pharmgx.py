"""
test_pharmgx.py — Automated test suite for PharmGx Reporter

Run with: pytest skills/pharmgx-reporter/tests/test_pharmgx.py -v

Uses the FIXED demo patient (demo_patient.txt) with known genotypes
so that all assertions are deterministic and reproducible.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from pharmgx_reporter import (
    PGX_SNPS,
    GENE_DEFS,
    GUIDELINES,
    _EVIDENCE_BADGE_CLASS,
    detect_format,
    parse_file,
    call_diplotype,
    call_phenotype,
    phenotype_to_key,
    lookup_drugs,
    generate_report,
    generate_html_report,
    enrich_with_clinpgx,
    _evidence_cell_html,
    _evidence_level_html,
)

DEMO = Path(__file__).parent.parent / "demo_patient.txt"


# ── Parsing ────────────────────────────────────────────────────────────────────

def test_detect_format_23andme():
    lines = DEMO.read_text().split("\n")
    assert detect_format(lines) == "23andme"


def test_parse_file_finds_all_pgx_snps():
    fmt, total_snps, pgx_snps = parse_file(str(DEMO))
    assert fmt == "23andme"
    assert total_snps == 21  # 21 PGx SNPs present on the 23andMe v5 chip (Corpasome)
    assert len(pgx_snps) == 21, (
        f"Expected 21 PGx SNPs (Corpasome v5 chip coverage), got {len(pgx_snps)}"
    )


def test_parse_file_genotype_values():
    _, _, pgx = parse_file(str(DEMO))
    # CYP2C19 *2 het
    assert pgx["rs4244285"]["genotype"] == "AG"
    # CYP2D6 *4 ref (Corpasome: no *4 variant)
    assert pgx["rs3892097"]["genotype"] == "CC"
    # VKORC1 hom variant (Corpasome: TT = high warfarin sensitivity)
    assert pgx["rs9923231"]["genotype"] == "TT"


# ── Star Allele Calling ───────────────────────────────────────────────────────

def _profiles():
    """Build profiles from demo patient for reuse across tests."""
    _, _, pgx = parse_file(str(DEMO))
    profiles = {}
    for gene in GENE_DEFS:
        diplotype = call_diplotype(gene, pgx)
        phenotype = call_phenotype(gene, diplotype)
        profiles[gene] = {"diplotype": diplotype, "phenotype": phenotype}
    return profiles


def test_cyp2c19_diplotype():
    """Demo patient: rs4244285 AG (*2 het) + rs12248560 CT (*17 het) → *17/*2."""
    p = _profiles()
    assert p["CYP2C19"]["diplotype"] == "*17/*2"


def test_cyp2d6_diplotype():
    """Demo patient: rs16947 AG (*2 het) + rs28371725 CT (*41 het) → *2/*41."""
    p = _profiles()
    assert p["CYP2D6"]["diplotype"] == "*2/*41"


def test_vkorc1_genotype():
    """Demo patient: rs9923231 TT → TT diplotype (hom variant)."""
    p = _profiles()
    assert p["VKORC1"]["diplotype"] == "TT"


def test_slco1b1_genotype():
    """Demo patient: rs4149056 TT → TT diplotype (ref, normal function)."""
    p = _profiles()
    assert p["SLCO1B1"]["diplotype"] == "TT"


def test_cyp3a5_diplotype():
    """Demo patient: rs776746 GG (*3 hom) → *3/*3."""
    p = _profiles()
    assert p["CYP3A5"]["diplotype"] == "*3/*3"


# ── Phenotype Assignment ──────────────────────────────────────────────────────

def test_cyp2c19_intermediate():
    p = _profiles()
    assert p["CYP2C19"]["phenotype"] == "Intermediate Metabolizer"


def test_cyp2d6_intermediate():
    """CYP2D6 *2/*41: *2 normal-function + *41 decreased-function → Intermediate."""
    p = _profiles()
    assert p["CYP2D6"]["phenotype"] == "Intermediate Metabolizer"


def test_vkorc1_high_sensitivity():
    """VKORC1 TT (hom variant) → High Warfarin Sensitivity per CPIC."""
    p = _profiles()
    assert p["VKORC1"]["phenotype"] == "High Warfarin Sensitivity"


def test_slco1b1_normal():
    """SLCO1B1 TT (ref) → Normal Function per CPIC."""
    p = _profiles()
    assert p["SLCO1B1"]["phenotype"] == "Normal Function"


def test_cyp3a5_nonexpressor():
    p = _profiles()
    assert p["CYP3A5"]["phenotype"] == "CYP3A5 Non-expressor"


def test_dpyd_normal():
    """All DPYD SNPs are ref → Normal Metabolizer."""
    p = _profiles()
    assert p["DPYD"]["phenotype"] == "Normal Metabolizer"


def test_tpmt_normal():
    p = _profiles()
    assert p["TPMT"]["phenotype"] == "Normal Metabolizer"


# ── Drug Recommendations ──────────────────────────────────────────────────────

def test_drug_lookup_returns_all_categories():
    p = _profiles()
    results = lookup_drugs(p)
    assert "standard" in results
    assert "caution" in results
    assert "avoid" in results
    total = sum(len(v) for v in results.values())
    assert total > 0


def test_clopidogrel_caution_for_intermediate():
    """CYP2C19 *1/*2 → Intermediate → Clopidogrel should be caution."""
    p = _profiles()
    results = lookup_drugs(p)
    clop = [d for d in results["caution"] if d["drug"] == "Clopidogrel"]
    assert len(clop) == 1, "Clopidogrel should be in caution list"


def test_codeine_caution_for_intermediate_cyp2d6():
    """CYP2D6 *2/*41 → Intermediate Metabolizer → Codeine should be caution."""
    p = _profiles()
    results = lookup_drugs(p)
    codeine = [d for d in results["caution"] if d["drug"] == "Codeine"]
    assert len(codeine) == 1, "Codeine should be in caution list for CYP2D6 IM"


def test_simvastatin_standard_for_normal_slco1b1():
    """SLCO1B1 TT → Normal Function → Simvastatin should be standard."""
    p = _profiles()
    results = lookup_drugs(p)
    simva = [d for d in results["standard"] if d["drug"] == "Simvastatin"]
    assert len(simva) == 1, "Simvastatin should be in standard list for SLCO1B1 Normal Function"


# ── Phenotype Key Mapping ─────────────────────────────────────────────────────

def test_phenotype_key_mapping():
    assert phenotype_to_key("Normal Metabolizer") == "normal_metabolizer"
    assert phenotype_to_key("Poor Metabolizer") == "poor_metabolizer"
    assert phenotype_to_key("High Warfarin Sensitivity") == "high_warfarin_sensitivity"
    assert phenotype_to_key("CYP3A5 Non-expressor") == "poor_metabolizer"
    assert phenotype_to_key("Normal (inferred)") == "normal_metabolizer"


# ── Report Generation ─────────────────────────────────────────────────────────

def test_report_contains_key_sections():
    _, _, pgx = parse_file(str(DEMO))
    p = _profiles()
    results = lookup_drugs(p)
    report = generate_report(str(DEMO), "23andme", 31, pgx, p, results)
    assert "# ClawBio PharmGx Report" in report
    assert "Drug Response Summary" in report
    assert "Gene Profiles" in report
    assert "Detected Variants" in report
    assert "Disclaimer" in report
    assert "Methods" in report
    assert "Reproducibility" in report


def test_report_contains_disclaimer():
    _, _, pgx = parse_file(str(DEMO))
    p = _profiles()
    results = lookup_drugs(p)
    report = generate_report(str(DEMO), "23andme", 31, pgx, p, results)
    assert "NOT a diagnostic device" in report


# ── Data Integrity ─────────────────────────────────────────────────────────────

def test_all_genes_have_phenotype_mappings():
    """Every gene in GENE_DEFS must have at least one phenotype."""
    for gene, gdef in GENE_DEFS.items():
        assert "phenotypes" in gdef, f"{gene} missing phenotypes"
        assert len(gdef["phenotypes"]) >= 2, f"{gene} has fewer than 2 phenotypes"


def test_all_guideline_drugs_reference_valid_genes():
    """Every drug in GUIDELINES must reference a gene in GENE_DEFS."""
    for drug, info in GUIDELINES.items():
        if info.get("special") == "warfarin":
            continue
        gene = info["gene"]
        assert gene in GENE_DEFS, f"{drug} references unknown gene {gene}"


# ── ClinPGx Evidence Enrichment ──────────────────────────────────────────────

def test_enrich_returns_dict():
    """enrich_with_clinpgx returns a dict even when ClinPGx is unavailable."""
    # Pass empty drug results — should return {} without error
    result = enrich_with_clinpgx({"standard": [], "caution": [], "avoid": [], "indeterminate": []})
    assert isinstance(result, {})  if False else True
    assert isinstance(result, dict)


def test_evidence_cell_html_empty():
    """Empty enrichment entry with classification renders fallback summary."""
    html = _evidence_cell_html({}, classification="caution")
    assert "Dose adjustment" in html
    assert "evidence-rec-text" in html


def test_evidence_cell_html_no_data():
    """No enrichment and no classification renders empty."""
    assert _evidence_cell_html({}) == ""


def test_evidence_cell_html_full():
    """Full enrichment entry renders multi-source recs with source acronyms."""
    entry = {
        "evidence_level": "1A",
        "sources": ["CPIC", "DPWG"],
        "verified": True,
        "guideline_name": "Test Guideline",
        "source_recs": [
            {"source": "CPIC", "rec": "Reduce dose by 50% for poor metabolizers.", "strength": "Strong"},
            {"source": "DPWG", "rec": "Use 75% of standard dose.", "strength": ""},
        ],
    }
    html = _evidence_cell_html(entry)
    assert "CPIC" in html
    assert "Reduce dose" in html
    assert "evidence-recs" in html
    assert "DPWG" in html  # second source
    assert "75% of standard dose" in html  # DPWG rec
    assert "title=" in html  # acronym tooltip


def test_evidence_level_html_verified():
    """Evidence level renders badge + checkmark."""
    entry = {"evidence_level": "1A", "verified": True}
    html = _evidence_level_html(entry)
    assert "1A" in html
    assert "badge-evidence-high" in html
    assert "&#10003;" in html


def test_evidence_level_html_unverified():
    """Unverified entry has no checkmark."""
    entry = {"evidence_level": "3", "verified": False}
    html = _evidence_level_html(entry)
    assert "&#10003;" not in html
    assert "badge-evidence-low" in html


def test_evidence_level_html_empty():
    """No enrichment returns empty string."""
    assert _evidence_level_html({}) == ""


def test_extract_phenotype_rec():
    """extract_phenotype_rec extracts matching recommendation from HTML table."""
    from clawbio.common.rec_shortener import extract_phenotype_rec
    html_table = """
    <table>
    <tr><th>Phenotype</th><th>Recommendation</th><th>Classification</th></tr>
    <tr><td>Normal Metabolizer</td><td>Use standard dose.</td><td>Strong</td></tr>
    <tr><td>Intermediate Metabolizer</td><td>Consider dose reduction.</td><td>Moderate</td></tr>
    <tr><td>Poor Metabolizer</td><td>Use alternative drug.</td><td>Strong</td></tr>
    </table>
    """
    rec, strength = extract_phenotype_rec(html_table, "Intermediate Metabolizer")
    assert rec == "Consider dose reduction."
    assert strength == "Moderate"


def test_extract_phenotype_rec_no_match():
    """Returns empty strings when phenotype not found."""
    from clawbio.common.rec_shortener import extract_phenotype_rec
    html_table = """
    <table>
    <tr><th>Phenotype</th><th>Recommendation</th><th>Classification</th></tr>
    <tr><td>Normal Metabolizer</td><td>Use standard dose.</td><td>Strong</td></tr>
    </table>
    """
    rec, strength = extract_phenotype_rec(html_table, "Poor Metabolizer")
    assert rec == ""
    assert strength == ""


def test_evidence_badge_class_mapping():
    """Badge class mapping covers all expected levels."""
    assert _EVIDENCE_BADGE_CLASS["1A"] == "badge-evidence-high"
    assert _EVIDENCE_BADGE_CLASS["1B"] == "badge-evidence-high"
    assert _EVIDENCE_BADGE_CLASS["2A"] == "badge-evidence-moderate"
    assert _EVIDENCE_BADGE_CLASS["2B"] == "badge-evidence-moderate"
    assert _EVIDENCE_BADGE_CLASS["3"] == "badge-evidence-low"
    assert _EVIDENCE_BADGE_CLASS["4"] == "badge-evidence-minimal"


def test_html_report_with_enrichment():
    """Evidence data renders when enrichment is provided."""
    _, _, pgx = parse_file(str(DEMO))
    p = _profiles()
    results = lookup_drugs(p)
    enrichment = {
        "clopidogrel": {
            "evidence_level": "1A", "sources": ["CPIC"], "verified": True,
            "source_recs": [
                {"source": "CPIC", "rec": "Use alternative antiplatelet therapy.", "strength": "Strong"},
            ],
        },
        "codeine": {
            "evidence_level": "1A", "sources": ["CPIC", "DPWG"], "verified": True,
            "source_recs": [
                {"source": "CPIC", "rec": "Use codeine label recommended dosing.", "strength": "Moderate"},
                {"source": "DPWG", "rec": "Monitor for reduced efficacy.", "strength": ""},
            ],
        },
    }
    html = generate_html_report(str(DEMO), "23andme", 21, pgx, p, results,
                                clinpgx_enrichment=enrichment)
    assert "badge-evidence-high" in html
    assert "&#10003;" in html  # checkmark
    assert "alternative antiplatelet" in html
    assert "evidence-recs" in html


def test_html_report_without_enrichment():
    """No evidence data when enrichment is None — still renders fine."""
    _, _, pgx = parse_file(str(DEMO))
    p = _profiles()
    results = lookup_drugs(p)
    html = generate_html_report(str(DEMO), "23andme", 21, pgx, p, results,
                                clinpgx_enrichment=None)
    # The body content should have no evidence badges (CSS classes exist in stylesheet, that's fine)
    body = html.split("<body>")[1]
    assert "badge-evidence-high" not in body
    assert "evidence-rec-source" not in body
