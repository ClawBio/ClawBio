"""
test_pharmgx.py — Automated test suite for PharmGx Reporter

Run with: pytest skills/pharmgx-reporter/tests/test_pharmgx.py -v

Uses the FIXED demo patient (demo_patient.txt) with known genotypes
so that all assertions are deterministic and reproducible.
"""

import re
import subprocess
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
    get_warfarin_rec,
    lookup_drugs,
    lookup_single_drug,
    format_dosage_card,
    generate_report,
    generate_html_report,
    enrich_with_clinpgx,
    write_commands_sh,
    _evidence_cell_html,
    _evidence_level_html,
)

DEMO = Path(__file__).parent.parent / "demo_patient.txt"


# ── Parsing ────────────────────────────────────────────────────────────────────

def test_detect_format_23andme():
    lines = DEMO.read_text().split("\n")
    assert detect_format(lines) == "23andme"


def test_parse_file_finds_all_pgx_snps():
    fmt, total_snps, pgx_snps, _ref = parse_file(str(DEMO))
    assert fmt == "23andme"
    assert total_snps == 23  # 23 PGx SNPs present on the 23andMe v5 chip (Corpasome, incl. MTHFR)
    assert len(pgx_snps) == 23, (
        f"Expected 23 PGx SNPs (Corpasome v5 chip coverage), got {len(pgx_snps)}"
    )


def test_parse_file_genotype_values():
    _, _, pgx, _ = parse_file(str(DEMO))
    # CYP2C19 *2 het
    assert pgx["rs4244285"]["genotype"] == "AG"
    # CYP2D6 *4 ref (Corpasome: no *4 variant)
    assert pgx["rs3892097"]["genotype"] == "CC"
    # VKORC1 hom variant (Corpasome: TT = high warfarin sensitivity)
    assert pgx["rs9923231"]["genotype"] == "TT"


# ── Star Allele Calling ───────────────────────────────────────────────────────

def _profiles():
    """Build profiles from demo patient for reuse across tests."""
    _, _, pgx, _ = parse_file(str(DEMO))
    profiles = {}
    for gene in GENE_DEFS:
        diplotype = call_diplotype(gene, pgx)
        phenotype = call_phenotype(gene, diplotype)
        profiles[gene] = {"diplotype": diplotype, "phenotype": phenotype}
    return profiles


def test_cyp2c19_diplotype():
    """Demo patient: rs4244285 AG (*2 het) + rs12248560 CT (*17 het) → phase ambiguity."""
    p = _profiles()
    assert "Indeterminate (phase ambiguity" in p["CYP2C19"]["diplotype"]


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
    """CYP2C19 *2 het + *17 het: phase ambiguity flagged as Indeterminate."""
    p = _profiles()
    assert "Indeterminate" in p["CYP2C19"]["phenotype"]


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


def test_mthfr_forward_strand_conversion_hom_1298():
    """677GG/1298GG (raw 23andMe fwd-strand) → CPIC 677CC/1298CC."""
    pgx = {
        "rs1801133": {"genotype": "GG", "gene": "MTHFR", "allele": "677T", "effect": "decreased_function"},
        "rs1801131": {"genotype": "GG", "gene": "MTHFR", "allele": "1298C", "effect": "decreased_function"},
    }
    diplotype = call_diplotype("MTHFR", pgx)
    assert diplotype == "677CC/1298CC", f"Expected 677CC/1298CC, got {diplotype}"


def test_mthfr_forward_strand_conversion_corpas():
    """Manuel Corpas fwd-strand rs1801133=GG, rs1801131=GT → CPIC 677CC/1298AC → Intermediate Activity."""
    pgx = {
        "rs1801133": {"genotype": "GG", "gene": "MTHFR", "allele": "677T", "effect": "decreased_function"},
        "rs1801131": {"genotype": "GT", "gene": "MTHFR", "allele": "1298C", "effect": "decreased_function"},
    }
    diplotype = call_diplotype("MTHFR", pgx)
    assert diplotype == "677CC/1298AC", f"Expected 677CC/1298AC, got {diplotype}"
    assert call_phenotype("MTHFR", diplotype) == "Reduced MTHFR enzyme activity (677CT)"


def test_mthfr_level2_reduced():
    """MTHFR 677TT/1298AA: DPWG guideline: Strongly reduced MTHFR enzyme activity."""
    assert call_phenotype("MTHFR", "677TT/1298AA") == "Strongly reduced MTHFR enzyme activity (677TT)"


def test_mthfr_level2_intermediate():
    """MTHFR 677CT/1298AA: DPWG guideline: Reduced MTHFR enzyme activity."""
    assert call_phenotype("MTHFR", "677CT/1298AA") == "Reduced MTHFR enzyme activity (677CT)"


def test_dpyd_normal_or_indeterminate():
    """DPYD with no risk alleles: Normal if fully covered, Indeterminate if partial."""
    p = _profiles()
    pheno = p["DPYD"]["phenotype"]
    assert "Normal" in pheno or "Indeterminate" in pheno, f"Unexpected DPYD phenotype: {pheno}"


def test_tpmt_normal_or_indeterminate():
    """TPMT with no risk alleles: Normal if fully covered, Indeterminate if partial."""
    p = _profiles()
    pheno = p["TPMT"]["phenotype"]
    assert "Normal" in pheno or "Indeterminate" in pheno, f"Unexpected TPMT phenotype: {pheno}"


# ── Drug Recommendations ──────────────────────────────────────────────────────

def test_drug_lookup_returns_all_categories():
    p = _profiles()
    results = lookup_drugs(p)
    assert "standard" in results
    assert "caution" in results
    assert "avoid" in results
    total = sum(len(v) for v in results.values())
    assert total > 0


def test_clopidogrel_indeterminate_for_phase_ambiguous():
    """CYP2C19 *2+*17 phase ambiguity -> Clopidogrel should be indeterminate."""
    p = _profiles()
    results = lookup_drugs(p)
    clop = [d for d in results["indeterminate"] if d["drug"] == "Clopidogrel"]
    assert len(clop) == 1, "Clopidogrel should be in indeterminate list when CYP2C19 is phase-ambiguous"


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


def test_warfarin_avoid_for_high_vkorc1_sensitivity():
    """Demo patient VKORC1 TT → High Warfarin Sensitivity → warfarin should be avoid."""
    p = _profiles()
    results = lookup_drugs(p)
    warfarin = [d for d in results.get("avoid", []) if d["drug"].lower() == "warfarin"]
    assert len(warfarin) == 1, "Warfarin should be in avoid list for High Warfarin Sensitivity"


# ── get_warfarin_rec ──────────────────────────────────────────────────────────

def test_get_warfarin_rec_returns_tuple():
    """get_warfarin_rec must always return a (classification, note) tuple."""
    result = get_warfarin_rec(_profiles())
    assert isinstance(result, tuple) and len(result) == 2


def test_get_warfarin_rec_standard():
    profiles = {
        "CYP2C9": {"phenotype": "Normal Metabolizer"},
        "VKORC1": {"phenotype": "Normal Sensitivity"},
    }
    cls, note = get_warfarin_rec(profiles)
    assert cls == "standard"
    assert note is None


def test_get_warfarin_rec_avoid_high_vkorc1():
    profiles = {
        "CYP2C9": {"phenotype": "Normal Metabolizer"},
        "VKORC1": {"phenotype": "High Warfarin Sensitivity"},
    }
    cls, note = get_warfarin_rec(profiles)
    assert cls == "avoid"
    assert note is None


def test_get_warfarin_rec_avoid_poor_cyp2c9():
    profiles = {
        "CYP2C9": {"phenotype": "Poor Metabolizer"},
        "VKORC1": {"phenotype": "Normal Sensitivity"},
    }
    cls, note = get_warfarin_rec(profiles)
    assert cls == "avoid"
    assert note is None


def test_get_warfarin_rec_indeterminate_missing_cyp2c9():
    """Missing CYP2C9 → indeterminate with a non-None note."""
    cls, note = get_warfarin_rec({"VKORC1": {"phenotype": "Normal Sensitivity"}})
    assert cls == "indeterminate"
    assert note is not None and "CYP2C9" in note


def test_get_warfarin_rec_indeterminate_missing_vkorc1():
    """Missing VKORC1 → indeterminate with a non-None note."""
    cls, note = get_warfarin_rec({"CYP2C9": {"phenotype": "Normal Metabolizer"}})
    assert cls == "indeterminate"
    assert note is not None and "VKORC1" in note


def test_lookup_drugs_warfarin_note_field_when_indeterminate():
    """When warfarin is indeterminate, the entry must carry a 'note' key."""
    p = _profiles()
    p["CYP2C9"] = {"diplotype": "?", "phenotype": ""}
    p["VKORC1"] = {"diplotype": "?", "phenotype": ""}
    results = lookup_drugs(p)
    warfarin_entries = [
        d for cat in results.values() for d in cat if d["drug"].lower() == "warfarin"
    ]
    assert len(warfarin_entries) == 1
    assert "note" in warfarin_entries[0]
    assert warfarin_entries[0]["note"]


def test_lookup_drugs_warfarin_no_note_when_avoid():
    """When warfarin is avoid (no note), the 'note' key should be absent."""
    profiles = {
        "CYP2C9": {"phenotype": "Normal Metabolizer"},
        "VKORC1": {"phenotype": "High Warfarin Sensitivity"},
    }
    # Merge with full profiles so non-warfarin lookups don't KeyError
    p = _profiles()
    p["CYP2C9"] = profiles["CYP2C9"]
    p["VKORC1"] = profiles["VKORC1"]
    results = lookup_drugs(p)
    warfarin_entries = [d for d in results.get("avoid", []) if d["drug"].lower() == "warfarin"]
    assert len(warfarin_entries) == 1
    assert "note" not in warfarin_entries[0]


# ── format_dosage_card note override ─────────────────────────────────────────

def _minimal_result(classification, note=None):
    r = {
        "drug": "Warfarin", "brand": "Coumadin", "class": "Anticoagulant",
        "gene": "CYP2C9 + VKORC1",
        "diplotype": "CYP2C9 *1/*1 / VKORC1 TT",
        "phenotype": "Normal / High Warfarin Sensitivity",
        "classification": classification,
    }
    if note:
        r["note"] = note
    return r


def test_format_dosage_card_uses_note_over_default_text():
    # Note is word-wrapped in the card; check a substring that fits on one line.
    custom_note = "CYP2C9 not genotyped. Clinical testing recommended."
    card = format_dosage_card(_minimal_result("indeterminate", note=custom_note))
    assert "CYP2C9 not genotyped" in card
    assert "Insufficient data" not in card  # default text should be suppressed


def test_format_dosage_card_falls_back_to_cls_text_without_note():
    card = format_dosage_card(_minimal_result("caution"))
    assert "Dose adjustment" in card


# ── Phenotype Key Mapping ─────────────────────────────────────────────────────

def test_phenotype_key_mapping():
    assert phenotype_to_key("Normal Metabolizer") == "normal_metabolizer"
    assert phenotype_to_key("Poor Metabolizer") == "poor_metabolizer"
    assert phenotype_to_key("High Warfarin Sensitivity") == "high_warfarin_sensitivity"
    assert phenotype_to_key("CYP3A5 Non-expressor") == "poor_metabolizer"
    assert phenotype_to_key("Normal (inferred)") == "normal_metabolizer"


# ── Report Generation ─────────────────────────────────────────────────────────

def test_report_contains_key_sections():
    _, _, pgx, _ = parse_file(str(DEMO))
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
    _, _, pgx, _ = parse_file(str(DEMO))
    p = _profiles()
    results = lookup_drugs(p)
    report = generate_report(str(DEMO), "23andme", 31, pgx, p, results)
    assert "NOT a diagnostic device" in report


def test_report_renders_standard_classification_token():
    """Standard-classified drugs must render the literal STANDARD token, not 'OK'.

    Regression for Tessl eval 019e83b8 scenario 1: the STANDARD bucket scored
    0/8 because standard drugs displayed as 'OK', so no STANDARD token appeared
    anywhere in the report (avoid/caution/indeterminate all rendered correctly).
    """
    _, _, pgx, _ = parse_file(str(DEMO))
    p = _profiles()
    results = lookup_drugs(p)
    assert len(results["standard"]) > 0, "demo patient should have >=1 standard drug"
    report = generate_report(str(DEMO), "23andme", 31, pgx, p, results)
    assert "STANDARD" in report, "report must label standard-dosing drugs as STANDARD"


# ── Reproducibility commands.sh ───────────────────────────────────────────────

def test_write_commands_sh_creates_file(tmp_path):
    """write_commands_sh must create reproducibility/commands.sh with the rerun command.

    Regression for Tessl eval 019e83b8 scenario 6: reproducibility/commands.sh
    scored 0/10 because the script never wrote the file the SKILL.md promises.
    """
    path = write_commands_sh(str(tmp_path), str(DEMO))
    assert path.exists(), "commands.sh should be written to disk"
    assert path.name == "commands.sh"
    assert path.parent.name == "reproducibility"
    content = path.read_text()
    assert "pharmgx_reporter.py" in content
    assert DEMO.name in content, "commands.sh should reference the input filename"


# ── Output contract (SKILL.md must match what the skill actually produces) ─────

def _parse_output_contract(skill_md):
    """Extract files promised in the SKILL.md '## Output Structure' tree.

    Directory lines, and any entry whose inline comment contains 'optional', are
    skipped. Returns [] when there is no parseable section.
    """
    if not skill_md.exists():
        return []
    m = re.search(r"##\s*Output Structure\s*\n+```[^\n]*\n(.*?)\n```",
                  skill_md.read_text(), re.S)
    if not m:
        return []
    files, parents = [], {}
    for raw in m.group(1).splitlines():
        if not raw.strip():
            continue
        parts = re.split(r"\s+#", raw, maxsplit=1)
        entry, comment = parts[0], (parts[1] if len(parts) > 1 else "")
        mm = re.match(r"^([\s│├└─]*)(.*)$", entry)
        prefix, name = mm.group(1), mm.group(2).strip()
        if not name:
            continue
        depth = len(prefix) // 4
        if depth == 0:
            continue
        if name.endswith("/"):
            parents[depth] = name.rstrip("/")
            for d in [k for k in parents if k > depth]:
                del parents[d]
            continue
        if "optional" in comment.lower():
            continue
        rel = "/".join(parents[d] for d in sorted(parents) if d < depth)
        files.append(rel + "/" + name if rel else name)
    return files


def test_documented_outputs_are_produced(tmp_path):
    """Every artifact in the SKILL.md Output Structure must be produced by a demo run.

    Deterministic, local capture of the doc/code drift the Tessl eval surfaced:
    reproducibility/commands.sh was documented but never written. This test would
    have failed before that fix.
    """
    skill_dir = Path(__file__).resolve().parent.parent
    promised = _parse_output_contract(skill_dir / "SKILL.md")
    assert promised, "expected a parseable '## Output Structure' section in SKILL.md"
    result = subprocess.run(
        [sys.executable, str(skill_dir / "pharmgx_reporter.py"),
         "--demo", "--output", str(tmp_path), "--no-enrich"],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, f"demo run failed: {result.stderr}"
    missing = [p for p in promised if not (tmp_path / p).exists()]
    assert not missing, (
        "SKILL.md Output Structure promises artifacts the skill did not produce: "
        f"{missing}"
    )


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


# ── CPIC Compliance (clawbio_bench v0.1.4 findings) ─────────────────────────

def test_cyp2c19_rapid_metabolizer_clopidogrel_standard():
    """CYP2C19 Rapid Metabolizer (*1/*17) should map to standard for Clopidogrel."""
    profiles = _profiles()
    profiles["CYP2C19"] = {"diplotype": "*1/*17", "phenotype": "Rapid Metabolizer"}
    results = lookup_drugs(profiles)
    clop = [d for cat in results.values() for d in cat if d["drug"] == "Clopidogrel"]
    assert clop[0]["classification"] == "standard", (
        "CYP2C19 Rapid Metabolizer should get standard for Clopidogrel per CPIC"
    )


def test_cyp2c19_rapid_metabolizer_voriconazole_caution():
    """CYP2C19 Rapid Metabolizer (*1/*17) should map to caution for Voriconazole."""
    profiles = _profiles()
    profiles["CYP2C19"] = {"diplotype": "*1/*17", "phenotype": "Rapid Metabolizer"}
    results = lookup_drugs(profiles)
    vori = [d for cat in results.values() for d in cat if d["drug"] == "Voriconazole"]
    assert vori[0]["classification"] == "caution", (
        "CYP2C19 Rapid Metabolizer should get caution for Voriconazole per CPIC"
    )


def test_diclofenac_no_avoid_for_cyp2c9_pm():
    """Diclofenac should NOT be 'avoid' for CYP2C9 PM (CPIC Theken 2020: no recommendation)."""
    profiles = _profiles()
    profiles["CYP2C9"] = {"diplotype": "*3/*3", "phenotype": "Poor Metabolizer"}
    results = lookup_drugs(profiles)
    diclo = [d for cat in results.values() for d in cat if d["drug"] == "Diclofenac"]
    assert diclo[0]["classification"] != "avoid", (
        "Diclofenac should not be avoid for CYP2C9 PM per CPIC Table S9 (Theken 2020)"
    )


def test_slco1b1_decreased_function_simvastatin_caution():
    """SLCO1B1 Decreased Function (TC) should map to caution for Simvastatin."""
    profiles = _profiles()
    profiles["SLCO1B1"] = {"diplotype": "TC", "phenotype": "Decreased Function"}
    results = lookup_drugs(profiles)
    simva = [d for cat in results.values() for d in cat if d["drug"] == "Simvastatin"]
    assert simva[0]["classification"] == "caution", (
        "SLCO1B1 Decreased Function should get caution for Simvastatin per CPIC 2022"
    )


def test_phenotype_to_key_decreased_function():
    """Decreased Function (CPIC 2022 SLCO1B1 term) must map correctly."""
    assert phenotype_to_key("Decreased Function") == "decreased_function"


def test_cyp2d6_star1_star10_is_normal_metabolizer():
    """CYP2D6 *1/*10 should be Normal Metabolizer per CPIC 2020 (AS >= 1.25)."""
    pheno = call_phenotype("CYP2D6", "*1/*10")
    assert pheno == "Normal Metabolizer", (
        f"CYP2D6 *1/*10 should be NM per CPIC 2020 (Caudle, PMID 31647186), got {pheno}"
    )


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
    _, _, pgx, _ = parse_file(str(DEMO))
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
    _, _, pgx, _ = parse_file(str(DEMO))
    p = _profiles()
    results = lookup_drugs(p)
    html = generate_html_report(str(DEMO), "23andme", 21, pgx, p, results,
                                clinpgx_enrichment=None)
    # The body content should have no evidence badges (CSS classes exist in stylesheet, that's fine)
    body = html.split("<body>")[1]
    assert "badge-evidence-high" not in body
    assert "evidence-rec-source" not in body


# --------------------------------------------------------------------------
# Reference-build detection
#
# Expected coordinates are asserted against dbSNP/Ensembl (see the citation
# block on _PGX_BUILD_POSITIONS), NOT read back from the module under test.
# Deriving them from the module would be tautological: the test would pass for
# any pair of values, including the GRCh37 contamination it exists to catch.
# --------------------------------------------------------------------------

# rsid -> (chrom, GRCh37 pos, GRCh38 pos). Verified against Ensembl REST
# (rest.ensembl.org for GRCh38, grch37.rest.ensembl.org for GRCh37).
_BUILD_TRUTH = {
    "rs4244285":  ("10", 96541616, 94781859),   # CYP2C19*2
    "rs3892097":  ("22", 42524947, 42128945),   # CYP2D6*4
    "rs1799853":  ("10", 96702047, 94942290),   # CYP2C9*2
    "rs9923231":  ("16", 31107689, 31096368),   # VKORC1 -1639G>A
    "rs1801133":  ("1",  11856378, 11796321),   # MTHFR C677T
}


def _positions(build_index, rsids=None, chrom=True):
    """Build a positions dict at GRCh37 (index 1) or GRCh38 (index 2)."""
    out = {}
    for rsid, row in _BUILD_TRUTH.items():
        if rsids is not None and rsid not in rsids:
            continue
        out[rsid] = {"pos": row[build_index]}
        if chrom:
            out[rsid]["chrom"] = row[0]
    return out


def test_build_table_holds_both_builds_with_correct_coordinates():
    """Each panel SNP carries the true GRCh37 and GRCh38 coordinate."""
    from pharmgx_reporter import _PGX_BUILD_POSITIONS
    for rsid, (chrom, pos37, pos38) in _BUILD_TRUTH.items():
        entry = _PGX_BUILD_POSITIONS[rsid]
        assert entry["chrom"] == chrom, f"{rsid} wrong chromosome"
        assert entry["GRCh37"] == pos37, f"{rsid} GRCh37 coordinate wrong"
        assert entry["GRCh38"] == pos38, f"{rsid} GRCh38 coordinate wrong"
        assert pos37 != pos38, f"{rsid} builds must differ to be discriminative"


def test_grch38_input_detected_as_grch38():
    from pharmgx_reporter import detect_reference_genome
    assert detect_reference_genome(_positions(2)) == "GRCh38"


def test_grch37_input_detected_as_grch37_not_unknown():
    """The assertion that fails against a GRCh37-contaminated table."""
    from pharmgx_reporter import detect_reference_genome
    assert detect_reference_genome(_positions(1)) == "GRCh37"


def test_single_snp_is_enough_to_identify_the_build():
    """A sparse targeted panel must not silently fall through to GRCh38."""
    from pharmgx_reporter import detect_reference_genome
    assert detect_reference_genome(_positions(1, {"rs9923231"})) == "GRCh37"
    assert detect_reference_genome(_positions(2, {"rs9923231"})) == "GRCh38"


def test_right_position_on_wrong_chromosome_is_not_a_match():
    from pharmgx_reporter import detect_reference_genome
    pos = _positions(2)
    for rsid in pos:
        pos[rsid]["chrom"] = "21"       # no panel SNP lives on chr21
    assert detect_reference_genome(pos) == "unknown_build"


def test_chromosome_prefixes_and_missing_chrom_are_tolerated():
    from pharmgx_reporter import detect_reference_genome
    prefixed = _positions(2)
    for rsid in prefixed:
        prefixed[rsid]["chrom"] = "chr" + prefixed[rsid]["chrom"]
    assert detect_reference_genome(prefixed) == "GRCh38"
    assert detect_reference_genome(_positions(2, chrom=False)) == "GRCh38"


def test_coordinates_matching_neither_build_are_unknown():
    from pharmgx_reporter import detect_reference_genome
    junk = {rsid: {"chrom": row[0], "pos": 12345} for rsid, row in _BUILD_TRUTH.items()}
    assert detect_reference_genome(junk) == "unknown_build"


def test_no_usable_positions_returns_none():
    from pharmgx_reporter import detect_reference_genome
    assert detect_reference_genome({}) is None
    assert detect_reference_genome({"rs4244285": {"pos": 0}}) is None
    assert detect_reference_genome({"rs9999999": {"pos": 123}}) is None


def test_demo_file_is_grch37_and_still_yields_real_diplotypes():
    """End-to-end: allele calling is rsID-based, so GRCh37 input is callable.

    The demo is genuine GRCh37 23andMe data. Detecting its build must not
    collapse every gene to Indeterminate.
    """
    from pharmgx_reporter import call_diplotype, call_phenotype
    _, _, pgx, ref_genome = parse_file(str(DEMO))
    assert ref_genome == "GRCh37"

    calls = {gene: call_diplotype(gene, pgx) for gene in GENE_DEFS}
    # No gene may be withheld for a build reason.
    assert not [g for g, d in calls.items() if "build" in d.lower()]
    # CYP2C9 *1/*2 is determinate in this sample and is the call the blanket
    # GRCh37 override used to destroy.
    assert "Indeterminate" not in calls["CYP2C9"]
    assert "Indeterminate" not in call_phenotype("CYP2C9", calls["CYP2C9"])
    # Remaining Indeterminates must be scientific (e.g. CYP2C19 *2/*17 phase
    # ambiguity), never coordinate-related.
    determinate = [g for g, d in calls.items() if "Indeterminate" not in d]
    assert len(determinate) >= 8, f"only {len(determinate)} genes called: {calls}"


def test_grch38_lifted_demo_yields_the_same_calls_as_grch37(tmp_path):
    """Build affects only the reported build, never the calls."""
    from pharmgx_reporter import call_diplotype
    lifted = []
    for line in DEMO.read_text().splitlines():
        parts = line.split("\t")
        if len(parts) == 4 and parts[0] in _BUILD_TRUTH:
            parts[2] = str(_BUILD_TRUTH[parts[0]][2])
            line = "\t".join(parts)
        lifted.append(line)
    target = tmp_path / "demo_grch38.txt"
    target.write_text("\n".join(lifted) + "\n")

    _, _, pgx37, ref37 = parse_file(str(DEMO))
    _, _, pgx38, ref38 = parse_file(str(target))
    assert (ref37, ref38) == ("GRCh37", "GRCh38")
    for gene in GENE_DEFS:
        assert call_diplotype(gene, pgx37) == call_diplotype(gene, pgx38), gene


def test_unknown_build_still_forces_indeterminate_end_to_end(tmp_path):
    """Coordinates matching no known build remain a hard stop."""
    corrupt = []
    for line in DEMO.read_text().splitlines():
        parts = line.split("\t")
        if len(parts) == 4 and parts[0] in _BUILD_TRUTH:
            parts[2] = "999999"
            line = "\t".join(parts)
        corrupt.append(line)
    target = tmp_path / "corrupt.txt"
    target.write_text("\n".join(corrupt) + "\n")

    _, _, _, ref_genome = parse_file(str(target))
    assert ref_genome == "unknown_build"

    result = subprocess.run(
        [sys.executable, str(Path(__file__).parent.parent / "pharmgx_reporter.py"),
         "--input", str(target), "--output", str(tmp_path / "out")],
        capture_output=True, text=True,
    )
    combined = (result.stdout + result.stderr).lower()
    assert "unknown_build" in combined
