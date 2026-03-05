"""Unit tests for clawbio.common.rec_shortener (structured table parser)."""

from clawbio.common.rec_shortener import extract_phenotype_rec, extract_all_recs_from_guidelines


# ── extract_phenotype_rec ────────────────────────────────────────────────────

class TestExtractPhenotypeRec:
    """Test extracting recommendations from HTML guideline tables."""

    SAMPLE_TABLE = """
    <table>
    <tr><th>Phenotype</th><th>Implication</th><th>Therapeutic Recommendation</th><th>Classification of Recommendation</th></tr>
    <tr><td>CYP2D6 Ultrarapid Metabolizer</td><td>Increased formation of morphine</td><td>Avoid codeine. Use alternative analgesic.</td><td>Strong</td></tr>
    <tr><td>CYP2D6 Normal Metabolizer</td><td>Normal morphine formation</td><td>Use codeine label recommended dosing.</td><td>Strong</td></tr>
    <tr><td>CYP2D6 Intermediate Metabolizer</td><td>Reduced morphine formation</td><td>Use label recommended dosing. If no response, consider non-tramadol opioid.</td><td>Moderate</td></tr>
    <tr><td>CYP2D6 Poor Metabolizer</td><td>Greatly reduced morphine formation</td><td>Avoid codeine. Use alternative analgesic.</td><td>Strong</td></tr>
    </table>
    """

    def test_intermediate_metabolizer(self):
        rec, strength = extract_phenotype_rec(self.SAMPLE_TABLE, "Intermediate Metabolizer", "CYP2D6")
        assert "label recommended dosing" in rec
        assert strength == "Moderate"

    def test_poor_metabolizer(self):
        rec, strength = extract_phenotype_rec(self.SAMPLE_TABLE, "Poor Metabolizer", "CYP2D6")
        assert "Avoid codeine" in rec
        assert strength == "Strong"

    def test_normal_metabolizer(self):
        rec, strength = extract_phenotype_rec(self.SAMPLE_TABLE, "Normal Metabolizer", "CYP2D6")
        assert "label recommended dosing" in rec

    def test_no_match(self):
        rec, strength = extract_phenotype_rec(self.SAMPLE_TABLE, "Rapid Metabolizer", "CYP2D6")
        assert rec == ""
        assert strength == ""

    def test_empty_html(self):
        assert extract_phenotype_rec("", "Normal Metabolizer") == ("", "")

    def test_empty_phenotype(self):
        assert extract_phenotype_rec(self.SAMPLE_TABLE, "") == ("", "")

    def test_no_table(self):
        assert extract_phenotype_rec("<p>No tables here</p>", "Normal Metabolizer") == ("", "")


class TestExtractAllRecsFromGuidelines:
    """Test batch extraction from guideline list."""

    def test_cpic_prioritised(self):
        guidelines = [
            {
                "source": "DPWG",
                "name": "codeine guideline",
                "textMarkdown": {"html": """
                    <table>
                    <tr><th>Phenotype</th><th>Recommendation</th></tr>
                    <tr><td>IM</td><td>DPWG: monitor closely.</td></tr>
                    </table>
                """},
            },
            {
                "source": "CPIC",
                "name": "codeine guideline",
                "textMarkdown": {"html": """
                    <table>
                    <tr><th>Phenotype</th><th>Recommendation</th><th>Classification</th></tr>
                    <tr><td>Intermediate Metabolizer</td><td>Use label dosing.</td><td>Moderate</td></tr>
                    </table>
                """},
            },
        ]
        rec, strength, source = extract_all_recs_from_guidelines(
            guidelines, "codeine", "Intermediate Metabolizer", "CYP2D6"
        )
        assert rec == "Use label dosing."
        assert source == "CPIC"

    def test_no_match_returns_empty(self):
        rec, strength, source = extract_all_recs_from_guidelines(
            [], "codeine", "Intermediate Metabolizer"
        )
        assert rec == ""
        assert strength == ""
        assert source == ""

    def test_fallback_to_dpwg(self):
        guidelines = [
            {
                "source": "DPWG",
                "name": "codeine guideline",
                "textMarkdown": {"html": """
                    <table>
                    <tr><th>Phenotype</th><th>Recommendation</th></tr>
                    <tr><td>Intermediate Metabolizer</td><td>Monitor for reduced efficacy.</td></tr>
                    </table>
                """},
            },
        ]
        rec, strength, source = extract_all_recs_from_guidelines(
            guidelines, "codeine", "Intermediate Metabolizer"
        )
        assert rec == "Monitor for reduced efficacy."
        assert source == "DPWG"


# ── Phenotype matching edge cases ────────────────────────────────────────────

class TestPhenotypeMatching:
    """Test various phenotype string formats."""

    def test_vkorc1_sensitivity(self):
        html = """
        <table>
        <tr><th>VKORC1 Phenotype</th><th>Recommendation</th></tr>
        <tr><td>High Warfarin Sensitivity</td><td>Reduce warfarin dose significantly.</td></tr>
        <tr><td>Low Warfarin Sensitivity</td><td>Standard warfarin dosing.</td></tr>
        </table>
        """
        rec, _ = extract_phenotype_rec(html, "High Warfarin Sensitivity", "VKORC1")
        assert "Reduce" in rec

    def test_cyp3a5_expressor(self):
        html = """
        <table>
        <tr><th>Phenotype</th><th>Recommendation</th></tr>
        <tr><td>CYP3A5 Expressor</td><td>Increase starting dose.</td></tr>
        <tr><td>CYP3A5 Non-expressor</td><td>Use standard dose.</td></tr>
        </table>
        """
        rec, _ = extract_phenotype_rec(html, "CYP3A5 Non-expressor", "CYP3A5")
        assert "standard dose" in rec

    def test_slco1b1_function(self):
        html = """
        <table>
        <tr><th>Phenotype</th><th>Recommendation</th></tr>
        <tr><td>Normal Function</td><td>Standard statin dosing.</td></tr>
        <tr><td>Decreased Function</td><td>Lower starting dose.</td></tr>
        </table>
        """
        rec, _ = extract_phenotype_rec(html, "Normal Function", "SLCO1B1")
        assert "Standard" in rec
