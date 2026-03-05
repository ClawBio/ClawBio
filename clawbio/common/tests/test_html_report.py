"""Unit tests for clawbio.common.html_report."""

import tempfile
from pathlib import Path

import pytest

from clawbio.common.html_report import HtmlReportBuilder, write_html_report
from clawbio.common.report import DISCLAIMER


def _minimal_builder():
    return HtmlReportBuilder("Test Report", "test-skill")


class TestHtmlSkeleton:
    def test_doctype_and_closing_tags(self):
        out = _minimal_builder().render()
        assert out.startswith("<!DOCTYPE html>")
        assert "</html>" in out
        assert "<title>Test Report</title>" in out

    def test_skill_in_subtitle(self):
        out = _minimal_builder().render()
        assert "test-skill" in out


class TestMetadata:
    def test_metadata_rendering(self):
        out = _minimal_builder().add_metadata({"File": "demo.txt", "Format": "23andMe"}).render()
        assert "<strong>File:</strong>" in out
        assert "demo.txt" in out
        assert "Format" in out

    def test_html_escape_in_metadata(self):
        out = _minimal_builder().add_metadata({"Key": "<script>alert(1)</script>"}).render()
        assert "<script>" not in out
        assert "&lt;script&gt;" in out


class TestSummaryCards:
    def test_cards_with_correct_classes(self):
        cards = [("Standard", 10, "standard"), ("Caution", 3, "caution"), ("Avoid", 1, "avoid")]
        out = _minimal_builder().add_summary_cards(cards).render()
        assert "summary-card standard" in out
        assert "summary-card caution" in out
        assert "summary-card avoid" in out
        assert ">10<" in out
        assert ">3<" in out
        assert ">1<" in out


class TestTable:
    def test_basic_table(self):
        out = _minimal_builder().add_table(["A", "B"], [["1", "2"], ["3", "4"]]).render()
        assert "<th>A</th>" in out
        assert "<td>1</td>" in out
        assert "<td>4</td>" in out

    def test_badge_column(self):
        out = _minimal_builder().add_table(
            ["Drug", "Status"],
            [["Codeine", "avoid"], ["Aspirin", "standard"]],
            badge_col=1,
        ).render()
        assert "badge badge-avoid" in out
        assert "badge badge-standard" in out
        # Should not have raw text "avoid" as plain td
        assert "<td>avoid</td>" not in out


class TestAlertBox:
    def test_avoid_alert(self):
        out = _minimal_builder().add_alert_box("avoid", "Stop!", "Do not use this.").render()
        assert "alert-box-avoid" in out
        assert "Stop!" in out
        assert "Do not use this." in out

    def test_caution_alert(self):
        out = _minimal_builder().add_alert_box("caution", "Warning", "Be careful.").render()
        assert "alert-box-caution" in out


class TestDisclaimer:
    def test_disclaimer_present(self):
        out = _minimal_builder().add_disclaimer().render()
        assert DISCLAIMER.split(".")[0] in out
        assert "disclaimer" in out


class TestWriteHtmlReport:
    def test_writes_to_disk(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = write_html_report(tmpdir, "test.html", "<html></html>")
            assert path.exists()
            assert path.name == "test.html"
            assert path.read_text() == "<html></html>"

    def test_creates_parent_dirs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            nested = Path(tmpdir) / "a" / "b"
            path = write_html_report(nested, "out.html", "content")
            assert path.exists()


class TestHeaderBlock:
    def test_header_renders_gradient_div(self):
        out = _minimal_builder().add_header_block("My Report").render()
        assert "report-header" in out
        assert "My Report" in out

    def test_header_suppresses_default_h1(self):
        b = _minimal_builder()
        b.add_header_block("Custom Title")
        out = b.render()
        # The default <h1>Test Report</h1> should NOT appear outside the header div
        assert out.count("<h1>") == 1  # only inside the header block

    def test_subtitle_rendered(self):
        out = _minimal_builder().add_header_block("Title", "Subtitle text").render()
        assert "Subtitle text" in out
        assert "subtitle" in out


class TestDonutChart:
    def test_svg_present(self):
        out = _minimal_builder().add_donut_chart([("A", 3, "#ff0000"), ("B", 7, "#00ff00")]).render()
        assert "<svg" in out

    def test_segments_rendered(self):
        out = _minimal_builder().add_donut_chart([("Red", 2, "#ff0000"), ("Blue", 8, "#0000ff")]).render()
        assert "#ff0000" in out
        assert "#0000ff" in out

    def test_legend_items(self):
        out = _minimal_builder().add_donut_chart([("Avoid", 1, "#c62828"), ("Standard", 5, "#2e7d32")]).render()
        assert "Avoid" in out
        assert "Standard" in out
        assert "donut-legend-item" in out

    def test_zero_segment_skipped(self):
        out = _minimal_builder().add_donut_chart([("A", 0, "#ff0000"), ("B", 5, "#00ff00")]).render()
        # Zero-count segment should not produce a <circle> stroke
        assert 'stroke="#ff0000"' not in out
        assert 'stroke="#00ff00"' in out


class TestProgressBar:
    def test_bar_renders(self):
        out = _minimal_builder().add_progress_bar("SNPs found", 15, 30).render()
        assert "progress-bar-container" in out
        assert "fill-green" in out

    def test_percentage_calculation(self):
        out = _minimal_builder().add_progress_bar("Test", 15, 30).render()
        assert "50%" in out
        assert "15/30" in out

    def test_zero_max_no_crash(self):
        out = _minimal_builder().add_progress_bar("Empty", 0, 0).render()
        assert "0%" in out


class TestExecutiveSummary:
    def test_renders_grid(self):
        items = [("\u26a0", "1 drug to avoid", "High-risk interaction.")]
        out = _minimal_builder().add_executive_summary(items).render()
        assert "exec-summary-grid" in out
        assert "exec-summary" in out

    def test_items_appear(self):
        items = [("\u2705", "30 standard drugs", "Normal response.")]
        out = _minimal_builder().add_executive_summary(items).render()
        assert "30 standard drugs" in out
        assert "Normal response." in out


class TestTableWrapped:
    def test_mobile_wrapper(self):
        out = _minimal_builder().add_table_wrapped(["A"], [["1"]]).render()
        assert "table-wrap" in out

    def test_row_classes(self):
        out = _minimal_builder().add_table_wrapped(
            ["Drug", "Status"],
            [["Warfarin", "avoid"], ["Aspirin", "standard"]],
            row_classes=["row-avoid", "row-standard"],
        ).render()
        assert "row-avoid" in out
        assert "row-standard" in out

    def test_badge_still_works(self):
        out = _minimal_builder().add_table_wrapped(
            ["Drug", "Status"],
            [["Codeine", "avoid"]],
            badge_col=1,
        ).render()
        assert "badge badge-avoid" in out


class TestFooterBlock:
    def test_footer_renders(self):
        out = _minimal_builder().add_footer_block("PharmGx Reporter", "0.2.0").render()
        assert "report-footer" in out
        assert "ClawBio" in out
        assert "v0.2.0" in out

    def test_suppresses_default_footer(self):
        out = _minimal_builder().add_footer_block("Skill").render()
        assert "footer-brand" in out
        # Default footer uses &middot; — should not be present
        assert "&middot;" not in out


class TestPrintStylesheet:
    def test_print_css_present(self):
        out = _minimal_builder().render()
        assert "@media print" in out


class TestEvidenceCSS:
    def test_evidence_css_classes_present(self):
        out = _minimal_builder().render()
        assert "badge-evidence-high" in out
        assert "badge-evidence-moderate" in out
        assert "badge-evidence-low" in out
        assert "badge-evidence-minimal" in out
        assert "badge-evidence-na" in out
        assert "evidence-verified" in out
        assert "evidence-unverified" in out
        assert "evidence-source" in out
        assert "evidence-recs" in out
        assert "evidence-rec-source" in out
        assert "evidence-rec-text" in out


class TestBackwardCompatibility:
    def test_old_css_variables_aliased(self):
        out = _minimal_builder().render()
        assert "--clawbio-green" in out


class TestMethodChaining:
    def test_chaining(self):
        out = (
            HtmlReportBuilder("Chain", "test")
            .add_metadata({"A": "1"})
            .add_section("Section")
            .add_paragraph("text")
            .add_summary_cards([("X", 5, "standard")])
            .add_alert_box("info", "Note", "body")
            .add_table(["H"], [["R"]])
            .add_disclaimer()
            .render()
        )
        assert "<!DOCTYPE html>" in out
        assert "Chain" in out

    def test_full_premium_chaining(self):
        out = (
            HtmlReportBuilder("Premium", "test")
            .add_disclaimer()
            .add_header_block("Premium Report", "Subtitle")
            .add_executive_summary([("\u2705", "Good", "All clear.")])
            .add_donut_chart([("A", 3, "#ff0000"), ("B", 7, "#00ff00")])
            .add_summary_cards([("X", 5, "standard")])
            .add_progress_bar("Coverage", 8, 10)
            .add_table_wrapped(["H"], [["R"]], row_classes=["row-standard"])
            .add_details("More info", "<p>Hidden</p>")
            .add_footer_block("test-skill", "1.0")
            .render()
        )
        assert "<!DOCTYPE html>" in out
        assert "report-header" in out
        assert "<svg" in out
        assert "report-footer" in out
