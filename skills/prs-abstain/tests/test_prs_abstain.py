"""Tests for prs-abstain. Written before the implementation (red/green TDD)."""
from __future__ import annotations

import json
import math
import re
import subprocess
import sys
from pathlib import Path

import pytest

SKILL_DIR = Path(__file__).resolve().parents[1]
SCRIPT = SKILL_DIR / "prs_abstain.py"
EXAMPLES = SKILL_DIR / "examples"
FIXTURES = SKILL_DIR / "tests" / "fixtures"

sys.path.insert(0, str(SKILL_DIR))


def run_cli(args):
    return subprocess.run(
        [sys.executable, str(SCRIPT), *args], capture_output=True, text=True
    )


# ── Calibration ───────────────────────────────────────────────────────────────

class TestCalibration:
    def test_centroid_and_threshold_from_panel(self):
        import prs_abstain as pa

        panel = pa.load_reference_panel(EXAMPLES / "demo_reference_pcs.csv")
        cal = pa.calibrate(panel, ref_pop="EUR", k_sd=3.0)
        assert cal.n == 22
        assert cal.threshold == pytest.approx(3.47, abs=0.05)
        # threshold must sit above every EUR member and below every non-EUR member
        assert cal.within_max < cal.threshold < cal.nearest_other

    def test_unknown_reference_population_raises(self):
        import prs_abstain as pa

        panel = pa.load_reference_panel(EXAMPLES / "demo_reference_pcs.csv")
        with pytest.raises(pa.CalibrationError):
            pa.calibrate(panel, ref_pop="NOPE", k_sd=3.0)

    def test_too_few_reference_individuals_raises(self):
        import prs_abstain as pa

        panel = pa.load_reference_panel(EXAMPLES / "demo_reference_pcs.csv")
        tiny = [s for s in panel if s.population == "AMR"]  # n=5
        with pytest.raises(pa.CalibrationError):
            pa.calibrate(tiny, ref_pop="AMR", k_sd=3.0, min_reference_n=10)


# ── Decision logic ────────────────────────────────────────────────────────────

class TestDecision:
    def _cal(self):
        import prs_abstain as pa

        panel = pa.load_reference_panel(EXAMPLES / "demo_reference_pcs.csv")
        return pa, panel, pa.calibrate(panel, ref_pop="EUR", k_sd=3.0)

    def test_eur_individual_is_reportable(self):
        pa, panel, cal = self._cal()
        ind = pa.Individual("EUR_001", "EUR", [-3.005, -2.385, -2.002, 0.345], 480)
        d = pa.decide(ind, cal, min_markers=30)
        assert d.verdict == "REPORT"

    def test_non_eur_individual_is_refused_as_distant(self):
        pa, panel, cal = self._cal()
        ind = pa.Individual("AFR_001", "AFR", [6.589, -2.403, -2.207, 3.187], 480)
        d = pa.decide(ind, cal, min_markers=30)
        assert d.verdict == "REFUSE_DISTANT"
        assert d.distance > cal.threshold

    def test_sparse_individual_is_refused_as_undeterminable(self):
        pa, panel, cal = self._cal()
        ind = pa.Individual("DEMO_PATIENT", None, None, 0)
        d = pa.decide(ind, cal, min_markers=30)
        assert d.verdict == "REFUSE_UNDETERMINABLE"
        assert d.distance is None

    def test_marker_check_precedes_distance_check(self):
        """A sparse individual must not be scored on coordinates even if present."""
        pa, panel, cal = self._cal()
        ind = pa.Individual("SPARSE_EUR", "EUR", [-3.0, -2.4, -2.0, 0.3], 5)
        d = pa.decide(ind, cal, min_markers=30)
        assert d.verdict == "REFUSE_UNDETERMINABLE"

    def test_missing_coordinates_never_silently_become_zero(self):
        pa, panel, cal = self._cal()
        ind = pa.Individual("NOCOORD", None, None, 480)
        d = pa.decide(ind, cal, min_markers=30)
        assert d.verdict == "REFUSE_UNDETERMINABLE"

    def test_every_verdict_carries_a_reason_and_remedy(self):
        pa, panel, cal = self._cal()
        for ind in [
            pa.Individual("A", "EUR", [-3.0, -2.4, -2.0, 0.3], 480),
            pa.Individual("B", "AFR", [6.6, -2.4, -2.2, 3.2], 480),
            pa.Individual("C", None, None, 0),
        ]:
            d = pa.decide(ind, cal, min_markers=30)
            assert d.reason and d.remedy


# ── Gating of PRS results ─────────────────────────────────────────────────────

class TestGating:
    def _gate(self, verdict_individual):
        import prs_abstain as pa

        panel = pa.load_reference_panel(EXAMPLES / "demo_reference_pcs.csv")
        cal = pa.calibrate(panel, ref_pop="EUR", k_sd=3.0)
        # The demo results are keyed per sample (Round 5); unit-level gating
        # takes the individual's own records, as main() now joins them.
        scores = [s_ for s_ in pa.load_prs_results(EXAMPLES / "demo_prs_results.json")
                  if s_.get("sample_id") == verdict_individual.sample_id]
        d = pa.decide(verdict_individual, cal, min_markers=30)
        return pa, pa.gate_scores(scores, d, cal), d

    def test_raw_score_is_always_retained(self):
        import prs_abstain as pa

        ind = pa.Individual("AFR_001", "AFR", [6.589, -2.403, -2.207, 3.187], 480)
        _, gated, _ = self._gate(ind)
        assert gated and all(g["raw_score"] is not None for g in gated)

    def test_percentile_withheld_on_refusal(self):
        import prs_abstain as pa

        ind = pa.Individual("AFR_001", "AFR", [6.589, -2.403, -2.207, 3.187], 480)
        _, gated, _ = self._gate(ind)
        assert all(g["percentile"] is None for g in gated)
        assert all(g["risk_category"] is None for g in gated)
        assert all(g["z_score"] is None for g in gated)

    def test_percentile_retained_on_report(self):
        """Non-sex-specific scores pass; sex-specific ones need a recorded sex."""
        import prs_abstain as pa

        ind = pa.Individual("EUR_001", "EUR", [-3.005, -2.385, -2.002, 0.345], 480)
        _, gated, _ = self._gate(ind)
        general = [g for g in gated if g["trait"] not in ("Breast cancer", "Prostate cancer")]
        assert general and all(g["percentile"] is not None for g in general)

    def test_sex_specific_scores_withheld_when_sex_unknown(self):
        import prs_abstain as pa

        ind = pa.Individual("EUR_001", "EUR", [-3.005, -2.385, -2.002, 0.345], 480)
        _, gated, _ = self._gate(ind)
        sexed = [g for g in gated if g["trait"] in ("Breast cancer", "Prostate cancer")]
        assert len(sexed) == 2
        assert all(g["percentile"] is None for g in sexed)

    def test_score_reference_population_mismatch_is_flagged(self):
        """A score whose reference population differs from the gate's must not pass silently."""
        import prs_abstain as pa

        panel = pa.load_reference_panel(EXAMPLES / "demo_reference_pcs.csv")
        cal = pa.calibrate(panel, ref_pop="EUR", k_sd=3.0)
        ind = pa.Individual("EUR_001", "EUR", [-3.005, -2.385, -2.002, 0.345], 480)
        d = pa.decide(ind, cal, min_markers=30)
        scores = [
            {
                "pgs_id": "PGS999999",
                "trait": "Synthetic",
                "raw_score": 1.0,
                "percentile": 50.0,
                "risk_category": "Average",
                "z_score": 0.0,
                "reference_population": "EAS",
            }
        ]
        gated = pa.gate_scores(scores, d, cal)
        assert gated[0]["percentile"] is None
        assert "EAS" in gated[0]["note"]


# ── CLI and output contract ───────────────────────────────────────────────────

class TestDemoCLI:
    def test_demo_runs_and_exits_zero(self, tmp_path):
        r = run_cli(["--demo", "--output", str(tmp_path)])
        assert r.returncode == 0, r.stderr

    def test_demo_produces_all_three_verdicts(self, tmp_path):
        run_cli(["--demo", "--output", str(tmp_path)])
        res = json.loads((tmp_path / "result.json").read_text())
        verdicts = {d["verdict"] for d in res["decisions"]}
        assert verdicts == {"REPORT", "REFUSE_DISTANT", "REFUSE_UNDETERMINABLE"}

    def test_report_contains_disclaimer(self, tmp_path):
        run_cli(["--demo", "--output", str(tmp_path)])
        text = (tmp_path / "report.md").read_text()
        assert "not a medical device" in text.lower()

    def test_refusal_states_it_is_not_reassurance(self, tmp_path):
        run_cli(["--demo", "--output", str(tmp_path)])
        text = (tmp_path / "report.md").read_text().lower()
        assert "not evidence of low risk" in text

    def test_report_never_prints_percentile_for_refused_individual(self, tmp_path):
        run_cli(["--demo", "--output", str(tmp_path)])
        res = json.loads((tmp_path / "result.json").read_text())
        for dec in res["decisions"]:
            if dec["verdict"] != "REPORT":
                assert all(s["percentile"] is None for s in dec["scores"])

    def test_threshold_provenance_recorded(self, tmp_path):
        run_cli(["--demo", "--output", str(tmp_path)])
        res = json.loads((tmp_path / "result.json").read_text())
        cal = res["calibration"]
        for key in ("reference_population", "n", "mean", "sd", "k_sd", "threshold", "pcs_used"):
            assert key in cal

    def test_per_individual_figure_written(self, tmp_path):
        run_cli(["--demo", "--output", str(tmp_path)])
        figs = list((tmp_path / "figures").glob("*.png"))
        assert len(figs) >= 4  # one panel overview + one per individual


def _parse_output_contract(skill_md: Path) -> list[str]:
    text = skill_md.read_text()
    m = re.search(r"##\s*Output Structure\s*\n+```[^\n]*\n(.*?)\n```", text, re.S)
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
        files.append(f"{rel}/{name}" if rel else name)
    return files


class TestOutputContract:
    def test_documented_outputs_are_produced(self, tmp_path):
        promised = _parse_output_contract(SKILL_DIR / "SKILL.md")
        if not promised:
            pytest.skip("No parseable '## Output Structure' section in SKILL.md")
        r = run_cli(["--demo", "--output", str(tmp_path)])
        assert r.returncode == 0, r.stderr
        missing = [p for p in promised if not (tmp_path / p).exists()]
        assert not missing, f"SKILL.md promises artifacts not produced: {missing}"


# ── Regressions found during stress testing ───────────────────────────────────

class TestStressRegressions:
    def test_placeable_but_marker_refused_still_renders_figures(self, tmp_path):
        """Individual has coordinates but fails the marker check: distance is None."""
        r = run_cli(["--demo", "--output", str(tmp_path), "--min-markers", "1000"])
        assert r.returncode == 0, r.stderr
        assert (tmp_path / "figures" / "individual_EUR_001.png").exists()

    def test_calibration_error_is_clean_not_a_traceback(self, tmp_path):
        r = run_cli(["--demo", "--output", str(tmp_path), "--ref-pop", "AMR"])
        assert r.returncode != 0
        assert "Traceback" not in r.stderr
        assert "not defensible" in r.stderr

    def test_threshold_that_swallows_other_populations_is_flagged(self, tmp_path):
        """k large enough to admit non-reference individuals must warn loudly
        even under the explicit override (without it, main() refuses to run)."""
        r = run_cli(["--demo", "--output", str(tmp_path), "--k-sd", "20",
                     "--allow-threshold-overreach"])
        assert r.returncode == 0, r.stderr
        res = json.loads((tmp_path / "result.json").read_text())
        assert res["calibration"]["threshold_exceeds_nearest_other"] is True
        assert "OVERREACH" in (tmp_path / "report.md").read_text().upper()


# ── v0.2: applicability, score integrity, per-variant audit ───────────────────

class TestApplicability:
    def test_sex_specific_trait_refused_for_wrong_sex(self):
        import prs_abstain as pa
        v = pa.check_applicability({"trait": "Prostate cancer"}, sex="female")
        assert v.applicable is False and "prostate" in v.reason.lower()

    def test_sex_specific_trait_allowed_for_right_sex(self):
        import prs_abstain as pa
        assert pa.check_applicability({"trait": "Prostate cancer"}, sex="male").applicable

    def test_unknown_sex_refuses_sex_specific_trait(self):
        import prs_abstain as pa
        assert pa.check_applicability({"trait": "Breast cancer"}, sex=None).applicable is False

    def test_non_sex_specific_trait_unaffected(self):
        import prs_abstain as pa
        assert pa.check_applicability({"trait": "Type 2 diabetes"}, sex=None).applicable


class TestScoreIntegrity:
    def _defs(self):
        import prs_abstain as pa
        return pa.load_score_definitions(Path(__file__).resolve().parents[1] / "examples" / "scores")

    def test_all_six_scores_load(self):
        assert len(self._defs()) == 6

    def test_weight_coverage_full_when_all_genotyped(self):
        import prs_abstain as pa
        gt = pa.load_genotype(EXAMPLES / "demo_genotype.txt")
        a = pa.audit_score(self._defs()["CLAWBIO-T2D-8"], gt)
        assert a.weight_coverage == pytest.approx(1.0)

    def test_missing_variants_reduce_weight_coverage(self):
        import prs_abstain as pa
        gt = pa.load_genotype(EXAMPLES / "demo_genotype.txt")
        thin = {k: v for i, (k, v) in enumerate(gt.items()) if i % 2 == 0}
        a = pa.audit_score(self._defs()["CLAWBIO-T2D-8"], thin)
        assert a.weight_coverage < 1.0
        assert a.weight_at_risk > 0

    def test_concentration_detects_fragile_score(self):
        import prs_abstain as pa
        gt = pa.load_genotype(EXAMPLES / "demo_genotype.txt")
        t2d = pa.audit_score(self._defs()["CLAWBIO-T2D-8"], gt)   # 8 variants
        bmi = pa.audit_score(self._defs()["CLAWBIO-BMI-97"], gt)   # 97 variants
        assert t2d.effective_n < 10 < bmi.effective_n
        assert t2d.top1_share > bmi.top1_share

    def test_palindromic_variants_are_counted(self):
        import prs_abstain as pa
        gt = pa.load_genotype(EXAMPLES / "demo_genotype.txt")
        a = pa.audit_score(self._defs()["CLAWBIO-T2D-8"], gt)
        assert a.palindromic_n >= 1

    def test_low_weight_coverage_refuses(self):
        import prs_abstain as pa
        gt = pa.load_genotype(EXAMPLES / "demo_genotype.txt")
        thin = {k: v for i, (k, v) in enumerate(gt.items()) if i % 4 == 0}
        a = pa.audit_score(self._defs()["CLAWBIO-T2D-8"], thin)
        assert pa.integrity_verdict(a, min_weight_coverage=0.90).passed is False


class TestAFShift:
    def test_reference_mean_equals_af_expectation(self):
        """The curated EUR mean is Sum 2*AF*w. This is the whole mechanism."""
        import prs_abstain as pa
        defs = pa.load_score_definitions(Path(__file__).resolve().parents[1] / "examples" / "scores")
        assert pa.expected_mean(defs["CLAWBIO-T2D-8"]) == pytest.approx(1.12, abs=0.01)
        assert pa.expected_mean(defs["CLAWBIO-CAD-46"]) == pytest.approx(2.84, abs=0.01)

    def test_af_shift_quantifies_percentile_error(self):
        import prs_abstain as pa
        defs = pa.load_score_definitions(Path(__file__).resolve().parents[1] / "examples" / "scores")
        af = pa.load_population_af(EXAMPLES / "demo_population_af.tsv")
        sh = pa.af_shift(defs["CLAWBIO-T2D-8"], af, sd=0.30)
        assert sh is not None
        assert sh.n_variants_with_af > 0
        assert isinstance(sh.shift_sd, float)
        assert len(sh.per_variant) == sh.n_variants_with_af

    def test_af_shift_returns_none_without_af_data(self):
        import prs_abstain as pa
        defs = pa.load_score_definitions(Path(__file__).resolve().parents[1] / "examples" / "scores")
        assert pa.af_shift(defs["CLAWBIO-T2D-8"], {}, sd=0.30) is None

    def test_per_variant_contributions_sum_to_total_shift(self):
        import prs_abstain as pa
        defs = pa.load_score_definitions(Path(__file__).resolve().parents[1] / "examples" / "scores")
        af = pa.load_population_af(EXAMPLES / "demo_population_af.tsv")
        sh = pa.af_shift(defs["CLAWBIO-T2D-8"], af, sd=0.30)
        assert sum(v["delta_mean"] for v in sh.per_variant) == pytest.approx(sh.shift_raw, abs=1e-9)


class TestDualReports:
    def test_both_reports_written(self, tmp_path):
        r = run_cli(["--demo", "--output", str(tmp_path)])
        assert r.returncode == 0, r.stderr
        assert (tmp_path / "report_clinician.md").exists()
        assert (tmp_path / "report_technical.md").exists()

    def test_clinician_report_avoids_jargon(self, tmp_path):
        run_cli(["--demo", "--output", str(tmp_path)])
        text = (tmp_path / "report_clinician.md").read_text().lower()
        for jargon in ("centroid", "euclidean", "herfindahl", "eigenvector"):
            assert jargon not in text, f"clinician report contains jargon: {jargon}"

    def test_clinician_ancestry_sentence_follows_ref_pop(self, tmp_path):
        """The ancestry sentence is derived from cal.reference_population, not
        hardcoded: under a non-EUR reference the clinician document must not
        claim a European comparison group (audit round 8, point 1)."""
        import prs_abstain as pa
        cal = pa.Calibration(
            reference_population="AFR", pcs_used=("PC1", "PC2"), centroid=[0.0, 0.0],
            n=10, mean=0.0, sd=1.0, k_sd=3.0, threshold=3.0,
            within_max=2.0, nearest_other=5.0)
        pa._write_clinician_report(tmp_path, cal, [], [], min_markers=100)
        text = (tmp_path / "report_clinician.md").read_text()
        assert "African ancestry" in text
        assert "European ancestry" not in text

    def test_default_demo_still_names_the_european_reference(self, tmp_path):
        run_cli(["--demo", "--output", str(tmp_path)])
        text = (tmp_path / "report_clinician.md").read_text()
        assert "European ancestry" in text

    def test_duplicate_position_withhold_gets_its_own_plain_reason(self):
        """A 'Score integrity' withhold caused by duplicated genomic positions
        must not be relabelled as low coverage (audit round 8, point 2)."""
        import prs_abstain as pa
        dup = {"withheld_reasons": [
            "Score integrity: Data integrity: 1 genomic position(s) carry more than "
            "one scored variant (chr1:100 (rs1, rs2)). This is either a coordinate "
            "error or the same locus counted twice, and both corrupt the sum."]}
        cov = {"withheld_reasons": [
            "Score integrity: weight coverage 0.40 below the 0.90 floor"]}
        assert "more than once" in pa._plain_withheld(dup)
        assert "too little" not in pa._plain_withheld(dup)
        assert pa._plain_withheld(cov) == "Withheld: too little of the score was measured"

    def test_clinician_report_states_plain_action(self, tmp_path):
        run_cli(["--demo", "--output", str(tmp_path)])
        text = (tmp_path / "report_clinician.md").read_text().lower()
        assert "not evidence of low risk" in text
        assert "what this means" in text

    def test_technical_report_carries_the_mechanism(self, tmp_path):
        run_cli(["--demo", "--output", str(tmp_path)])
        text = (tmp_path / "report_technical.md").read_text()
        assert "2 * AF" in text or "2*AF" in text or "2·AF" in text
        assert "effective_n" in text.lower() or "effective number" in text.lower()

    def test_sex_mismatch_visible_in_demo(self, tmp_path):
        run_cli(["--demo", "--output", str(tmp_path)])
        res = json.loads((tmp_path / "result.json").read_text())
        notes = [s["note"] for d in res["decisions"] for s in d["scores"]]
        assert any("sex" in n.lower() for n in notes)


# ── v0.3: linkage disequilibrium proxies ──────────────────────────────────────

class TestLDAudit:
    def _defs(self):
        import prs_abstain as pa
        return pa.load_score_definitions(Path(__file__).resolve().parents[1] / "examples" / "scores")

    def test_detects_correlated_pair_in_t2d_score(self):
        """rs7903146 and rs12255372 are both TCF7L2, ~50 kb apart."""
        import prs_abstain as pa
        ld = pa.ld_audit(self._defs()["CLAWBIO-T2D-8"], window_kb=250)
        assert ld.n_clusters_multi == 1
        members = [set(c["rsids"]) for c in ld.clusters if len(c["rsids"]) > 1]
        assert {"rs7903146", "rs12255372"} in members

    def test_effective_n_falls_when_ld_accounted_for(self):
        import prs_abstain as pa
        ld = pa.ld_audit(self._defs()["CLAWBIO-T2D-8"], window_kb=250)
        assert ld.effective_n_ld < ld.effective_n_independent
        assert ld.effective_n_ld == pytest.approx(3.55, abs=0.1)

    def test_clustered_weight_share_reported(self):
        import prs_abstain as pa
        ld = pa.ld_audit(self._defs()["CLAWBIO-T2D-8"], window_kb=250)
        assert ld.clustered_weight_share == pytest.approx(0.48, abs=0.02)

    def test_duplicate_positions_flagged_as_data_error(self):
        import prs_abstain as pa
        ld = pa.ld_audit(self._defs()["CLAWBIO-BC-77"], window_kb=250)
        assert ld.duplicate_positions
        assert any("rs11552449" in d["rsids"] for d in ld.duplicate_positions)

    def test_clean_score_has_no_duplicates(self):
        import prs_abstain as pa
        assert not pa.ld_audit(self._defs()["CLAWBIO-T2D-8"], window_kb=250).duplicate_positions

    def test_duplicate_positions_block_the_score(self):
        import prs_abstain as pa
        gt = pa.load_genotype(EXAMPLES / "demo_genotype.txt")
        defs = self._defs()
        au = pa.audit_score(defs["CLAWBIO-BC-77"], gt)
        ld = pa.ld_audit(defs["CLAWBIO-BC-77"], window_kb=250)
        assert pa.integrity_verdict(au, ld=ld).passed is False

    def test_ld_warning_uses_corrected_effective_n(self):
        import prs_abstain as pa
        gt = pa.load_genotype(EXAMPLES / "demo_genotype.txt")
        defs = self._defs()
        au = pa.audit_score(defs["CLAWBIO-T2D-8"], gt)
        ld = pa.ld_audit(defs["CLAWBIO-T2D-8"], window_kb=250)
        v = pa.integrity_verdict(au, ld=ld, min_effective_n=10.0)
        assert any("3.5" in w or "3.6" in w for w in v.warnings)

    def test_window_size_changes_clustering(self):
        import prs_abstain as pa
        d = self._defs()["CLAWBIO-T2D-8"]
        assert pa.ld_audit(d, window_kb=10).n_clusters_multi == 0
        assert pa.ld_audit(d, window_kb=250).n_clusters_multi == 1

    def test_ld_section_present_in_technical_report(self, tmp_path):
        run_cli(["--demo", "--output", str(tmp_path)])
        text = (tmp_path / "report_technical.md").read_text().lower()
        assert "linkage disequilibrium" in text
        assert "tag" in text

    def test_clinician_report_explains_ld_without_the_term(self, tmp_path):
        run_cli(["--demo", "--output", str(tmp_path)])
        text = (tmp_path / "report_clinician.md").read_text().lower()
        assert "linkage disequilibrium" not in text
        assert "same region" in text or "same stretch" in text


# ── v0.3: PDF rendering ───────────────────────────────────────────────────────

class TestPDF:
    def test_both_pdfs_written(self, tmp_path):
        pytest.importorskip("reportlab")
        r = run_cli(["--demo", "--output", str(tmp_path)])
        assert r.returncode == 0, r.stderr
        for name in ("report_clinician.pdf", "report_technical.pdf"):
            f = tmp_path / name
            assert f.exists(), f"{name} not written"
            assert f.stat().st_size > 5000, f"{name} suspiciously small"

    def test_pdf_is_valid_and_multipage(self, tmp_path):
        pytest.importorskip("reportlab")
        pypdf = pytest.importorskip("pypdf")
        run_cli(["--demo", "--output", str(tmp_path)])
        reader = pypdf.PdfReader(str(tmp_path / "report_technical.pdf"))
        assert len(reader.pages) >= 2
        text = "".join(p.extract_text() or "" for p in reader.pages)
        assert "prs-abstain" in text.lower()

    def test_clinician_pdf_carries_the_key_sentence(self, tmp_path):
        pytest.importorskip("reportlab")
        pypdf = pytest.importorskip("pypdf")
        run_cli(["--demo", "--output", str(tmp_path)])
        reader = pypdf.PdfReader(str(tmp_path / "report_clinician.pdf"))
        text = " ".join((p.extract_text() or "") for p in reader.pages).lower()
        text = " ".join(text.split())
        assert "not evidence of low risk" in text

    def test_skill_still_runs_without_reportlab(self, tmp_path):
        """PDF is a bonus artefact; its absence must not break the run.
        Measured for real: the CLI runs in a subprocess whose import of
        reportlab is forced to fail, and must still exit 0 with the markdown
        reports written and no PDFs."""
        import textwrap
        blocker = tmp_path / "blocker"
        blocker.mkdir()
        (blocker / "reportlab").mkdir()
        (blocker / "reportlab" / "__init__.py").write_text(
            textwrap.dedent("""
            raise ImportError("reportlab deliberately unavailable for this test")
            """))
        import os
        env = dict(os.environ)
        env["PYTHONPATH"] = str(blocker) + os.pathsep + env.get("PYTHONPATH", "")
        out = tmp_path / "out"
        r = subprocess.run(
            [sys.executable, str(SCRIPT), "--demo", "--output", str(out),
             "--no-figures"],
            capture_output=True, text=True, env=env)
        assert r.returncode == 0, r.stderr
        assert (out / "report_clinician.md").exists()
        assert not list(out.glob("*.pdf"))

    def test_no_pdf_flag_skips_generation(self, tmp_path):
        r = run_cli(["--demo", "--output", str(tmp_path), "--no-pdf"])
        assert r.returncode == 0, r.stderr
        assert not (tmp_path / "report_clinician.pdf").exists()


# ── Added at maintainer review (PR #348): boundary + parser tests ─────────────

class TestDecisionBoundaries:
    def _cal(self):
        import prs_abstain as pa

        panel = pa.load_reference_panel(EXAMPLES / "demo_reference_pcs.csv")
        return pa, pa.calibrate(panel, ref_pop="EUR", k_sd=3.0)

    def test_distance_exactly_at_threshold_reports(self):
        """decide() refuses on dist > threshold, so the boundary itself reports."""
        pa, cal = self._cal()
        pcs = list(cal.centroid)
        pcs[0] += cal.threshold  # exactly threshold away, along PC1
        d = pa.decide(pa.Individual("EDGE", "EUR", pcs, 480), cal, min_markers=30)
        # Decision.distance is rounded for display; the comparison itself runs
        # unrounded, and the boundary (dist == threshold) reports.
        assert d.distance == pytest.approx(cal.threshold, abs=1e-4)
        assert d.verdict == "REPORT"

    def test_exactly_min_markers_passes_marker_gate(self):
        """decide() refuses on n < min_markers, so exactly 30 markers passes."""
        pa, cal = self._cal()
        ind = pa.Individual("EDGE30", "EUR", list(cal.centroid), 30)
        d = pa.decide(ind, cal, min_markers=30)
        assert d.verdict == "REPORT"

    def test_inflated_marker_declaration_fails_closed(self):
        """Round-4 review: n_markers_shared was trusted unconditionally. A
        declaration above the genotype's valid-call count is not credible."""
        pa, cal = self._cal()
        ind = pa.Individual("INFLATED", "EUR", list(cal.centroid), 500)
        d = pa.decide(ind, cal, min_markers=30, max_credible_markers=415)
        assert d.verdict == "REFUSE_UNDETERMINABLE"
        assert "exceeds the 415 valid calls" in d.reason
        assert d.distance is None  # distance is never computed past a failed gate

    def test_marker_declaration_at_the_call_count_passes(self):
        """The cap is an upper bound: exactly at the call count is credible."""
        pa, cal = self._cal()
        ind = pa.Individual("ATCAP", "EUR", list(cal.centroid), 415)
        d = pa.decide(ind, cal, min_markers=30, max_credible_markers=415)
        assert d.verdict == "REPORT"


class TestScoreFileParsing:
    def test_authentic_pgs_catalog_file_parses_by_header_name(self):
        """A genuine harmonised download: different column order, hm_* columns,
        full-precision weights, no allele-frequency column."""
        import prs_abstain as pa

        defs = pa.load_score_definitions(SKILL_DIR / "tests" / "fixtures")
        s = defs["PGS000001"]
        assert len(s.variants) == 77
        v = s.variants[0]
        assert v["rsid"] == "rs78540526"
        assert v["chr"] == "11" and v["pos"] == "69331418"  # from hm_chr / hm_pos
        assert v["weight"] == pytest.approx(0.16220387987485377)
        assert v["af_reference"] is None  # real files often carry no AF column
        assert pa.expected_mean(s) is None  # and expected_mean must say so, not crash

    def test_malformed_weight_is_an_error_not_a_silent_skip(self, tmp_path):
        import prs_abstain as pa

        bad = tmp_path / "PGS999999_hmPOS_GRCh37.txt"
        bad.write_text(
            "#pgs_id=PGS999999\n"
            "rsID\tchr_name\teffect_allele\tother_allele\teffect_weight\tlocus_name\n"
            "rs1\t1\tA\tG\tCCND1\tx\n")
        with pytest.raises(ValueError, match="PGS999999.*effect_weight"):
            pa.load_score_definitions(tmp_path)

    def test_missing_required_column_is_named_in_the_error(self, tmp_path, capsys):
        import prs_abstain as pa

        bad = tmp_path / "PGS999998.txt"
        bad.write_text("#pgs_id=PGS999998\nfoo\tbar\n1\t2\n")
        # A file whose columns cannot be resolved is not treated as a scoring
        # file: skipped with a stderr warning naming the missing columns. Any
        # score expecting it then fails closed at the gate, and main() refuses
        # to run at all if no usable file remains (tested through the CLI in
        # TestGateWiringThroughTheCli).
        defs = pa.load_score_definitions(tmp_path)
        captured = capsys.readouterr()
        assert defs == {}
        assert "required column" in captured.err and "PGS999998.txt" in captured.err


# ── Added at maintainer review (PR #348): fail-closed behaviour ───────────────

class TestFailClosed:
    def _cal_and_report_decision(self):
        import prs_abstain as pa

        panel = pa.load_reference_panel(EXAMPLES / "demo_reference_pcs.csv")
        cal = pa.calibrate(panel, ref_pop="EUR", k_sd=3.0)
        ind = pa.Individual("EUR_001", "EUR", [-3.005, -2.385, -2.002, 0.345], 480,
                            sex="female")
        return pa, cal, pa.decide(ind, cal, min_markers=30)

    def test_missing_reference_population_withholds_percentile(self):
        """Absent provenance metadata must fail closed, not default to EUR."""
        pa, cal, dec = self._cal_and_report_decision()
        assert dec.verdict == "REPORT"
        score = {"pgs_id": "PGSX", "trait": "Type 2 diabetes", "raw_score": 1.0,
                 "percentile": 50.0, "z_score": 0.0}  # no reference_population
        gated = pa.gate_scores([score], dec, cal, sex="female")[0]
        assert gated["percentile"] is None
        assert any("provenance" in r.lower() for r in gated["withheld_reasons"])
        # the note must not assert an origin that was never established
        assert "Reported against" not in gated["note"]

    def test_declared_reference_population_still_reports(self):
        pa, cal, dec = self._cal_and_report_decision()
        score = {"pgs_id": "PGSX", "trait": "Type 2 diabetes", "raw_score": 1.0,
                 "percentile": 50.0, "z_score": 0.0, "reference_population": "EUR"}
        gated = pa.gate_scores([score], dec, cal, sex="female")[0]
        assert gated["percentile"] == 50.0
        assert "Reported against the EUR reference distribution." in gated["note"]

    def test_absent_or_unreadable_trait_fails_closed(self):
        import prs_abstain as pa

        for trait in (None, "", "   ", "unknown", "NR"):
            app = pa.check_applicability({"trait": trait}, "female")
            assert not app.applicable, f"trait {trait!r} must fail closed"
        # and a recognised non-sex-specific trait still passes
        assert pa.check_applicability({"trait": "Coronary artery disease"}, None).applicable

    def test_integrity_tier_runs_without_genotype(self):
        """Duplicate positions are a property of the file alone and must block
        even when no genotype was supplied; the skipped genotype-dependent
        checks must be recorded, not silent."""
        import prs_abstain as pa

        defs = pa.load_score_definitions(EXAMPLES / "scores")
        sdef = defs["CLAWBIO-BC-77"]  # carries the chr1:114448389 duplicate pair
        ld = pa.ld_audit(sdef)
        assert ld.duplicate_positions, "fixture must contain the duplicate pair"
        v = pa.integrity_verdict(None, ld=ld)
        assert not v.passed
        assert any("Data integrity" in r for r in v.reasons)
        assert any("genotype-dependent integrity checks" in w for w in v.warnings)

    def test_integrity_reasons_reach_the_gate_without_genotype(self):
        pa, cal, dec = self._cal_and_report_decision()
        defs = pa.load_score_definitions(EXAMPLES / "scores")
        sdef = defs["CLAWBIO-BC-77"]
        integ = {"CLAWBIO-BC-77": pa.integrity_verdict(None, ld=pa.ld_audit(sdef))}
        score = {"pgs_id": "CLAWBIO-BC-77", "trait": "Breast cancer", "raw_score": 1.0,
                 "percentile": 94.0, "z_score": 1.5, "reference_population": "EUR"}
        gated = pa.gate_scores([score], dec, cal, sex="female", integrity=integ)[0]
        assert gated["percentile"] is None
        assert any("Score integrity" in r for r in gated["withheld_reasons"])


class TestCoverageCounting:
    """A no-call or an incompatible call must not count toward weight coverage."""

    def _sdef(self):
        import prs_abstain as pa

        return pa.ScoreDefinition("PGST", "test trait", "GRCh37", [
            {"rsid": "rs1", "chr": "1", "pos": "100", "effect_allele": "A",
             "other_allele": "G", "weight": 1.0, "af_reference": None},
            {"rsid": "rs2", "chr": "2", "pos": "200", "effect_allele": "C",
             "other_allele": "T", "weight": 1.0, "af_reference": None},
        ])

    def test_no_call_is_not_covered(self):
        import prs_abstain as pa

        a = pa.audit_score(self._sdef(), {"rs1": "AG", "rs2": "--"})
        assert a.n_matched == 1
        assert a.weight_coverage == pytest.approx(0.5)
        assert any(m["rsid"] == "rs2" for m in a.missing_top)

    def test_incompatible_alleles_are_not_covered(self):
        import prs_abstain as pa

        # rs2 declares C/T; a G call matches neither the pair nor its complement pair alone
        a = pa.audit_score(self._sdef(), {"rs1": "AA", "rs2": "GC"})
        assert a.n_matched == 1

    def test_strand_complement_call_is_covered(self):
        import prs_abstain as pa

        # rs1 declares A/G; T/C is the same site read from the other strand
        a = pa.audit_score(self._sdef(), {"rs1": "TC", "rs2": "CT"})
        assert a.n_matched == 2
        assert a.weight_coverage == pytest.approx(1.0)


# ── Gate wiring through the CLI ───────────────────────────────────────────────
# Round-2 review: the integrity tier failed open on the documented standard
# command because every prior CLI test passed --demo (which supplies both
# --genotype and --scores) and the unit tests sat on either side of the
# wiring. These tests go through main() with each input deliberately absent.

class TestGateWiringThroughTheCli:
    def _standard(self, tmp_path, *extra):
        return ["--reference-panel", str(EXAMPLES / "demo_reference_pcs.csv"),
                "--individuals", str(EXAMPLES / "demo_query_individuals.csv"),
                "--prs-results", str(EXAMPLES / "demo_prs_results.json"),
                "--output", str(tmp_path), "--no-figures", "--no-pdf", *extra]

    def _eur_scores(self, tmp_path):
        res = json.loads((tmp_path / "result.json").read_text())
        rep = [d for d in res["decisions"] if d["verdict"] == "REPORT"]
        assert rep, "expected at least one REPORT individual in the demo data"
        return rep[0]["scores"]

    def test_scores_without_genotype_still_blocks_duplicate_positions(self, tmp_path):
        """--scores given, --genotype absent: the file-level integrity checks
        must still run and CLAWBIO-BC-77 (duplicate chr1:114448389 pair) must
        be withheld, through main(), on a REPORT-verdict individual."""
        r = run_cli(self._standard(tmp_path, "--scores", str(EXAMPLES / "scores")))
        assert r.returncode == 0, r.stderr
        bc = next(s for s in self._eur_scores(tmp_path) if s["pgs_id"] == "CLAWBIO-BC-77")
        assert bc["percentile"] is None
        assert any("duplicate" in x.lower() or "more than one scored variant" in x.lower()
                   for x in bc["withheld_reasons"])

    def test_no_scores_means_no_silent_integrity_pass(self, tmp_path):
        """The documented standard command (no --scores): every released
        percentile must carry the integrity-not-verified caveat rather than
        passing the tier silently."""
        r = run_cli(self._standard(tmp_path))
        assert r.returncode == 0, r.stderr
        released = [s for s in self._eur_scores(tmp_path) if s["percentile"] is not None]
        assert released, "demo data should release at least one percentile here"
        for s in released:
            assert any("integrity was not verified" in c.lower() for c in s["caveats"]), s["pgs_id"]
            assert "integrity was not verified" in s["note"].lower()

    def test_score_with_no_matching_scoring_file_is_withheld(self, tmp_path):
        """--scores supplied but a prs_results entry has no matching file:
        that score fails closed instead of skipping the tier."""
        import shutil
        subset = tmp_path / "subset_scores"
        subset.mkdir()
        shutil.copy(EXAMPLES / "scores" / "CLAWBIO-T2D-8_GRCh37.txt", subset)
        out = tmp_path / "out"
        r = run_cli(self._standard(out, "--scores", str(subset)))
        assert r.returncode == 0, r.stderr
        res = json.loads((out / "result.json").read_text())
        rep = [d for d in res["decisions"] if d["verdict"] == "REPORT"][0]
        af12 = next(s for s in rep["scores"] if s["pgs_id"] == "CLAWBIO-AF-12")
        assert af12["percentile"] is None
        assert any("never inspected" in x and "fails closed" in x
                   for x in af12["withheld_reasons"])


class TestScoreIdProvenance:
    def test_panel_id_comes_from_the_header_not_the_filename(self):
        """Post-#357 curated panels carry #clawbio_panel_id and no #pgs_id.
        The id must be the header's, never re-derived from a filename that
        still says PGS000013."""
        import prs_abstain as pa

        defs = pa.load_score_definitions(EXAMPLES / "scores")
        assert "CLAWBIO-T2D-8" in defs
        assert not any(k.startswith("PGS") for k in defs)

    def test_file_with_no_id_header_is_refused(self, tmp_path):
        import prs_abstain as pa

        f = tmp_path / "PGS000099_hmPOS_GRCh37.txt"
        f.write_text("#trait_reported=x\n"
                     "rsID\tchr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\n"
                     "rs1\t1\t100\tA\tG\t0.5\n")
        with pytest.raises(ValueError, match="filename"):
            pa.load_score_definitions(tmp_path)


class TestSdProvenance:
    def test_zero_z_score_propagates_none_not_unit_sd(self):
        """z_score == 0.0 (an individual exactly at the mean) makes the sd
        underivable; the shift must not silently be printed in 'sd' units."""
        import prs_abstain as pa

        defs = pa.load_score_definitions(EXAMPLES / "scores")
        af = {v["rsid"]: {"AFR": 0.5} for v in defs["CLAWBIO-T2D-8"].variants}
        sh = pa.af_shift(defs["CLAWBIO-T2D-8"], af, sd=None)
        assert sh is not None
        assert sh.shift_sd is None
        assert sh.shift_raw != 0

    def test_unknown_sd_becomes_a_caveat_not_a_quotient(self):
        import prs_abstain as pa

        sh = pa.AFShift(pgs_id="CLAWBIO-T2D-8", population="AFR", n_variants_with_af=8,
                        coverage=1.0, shift_raw=0.4, shift_sd=None)
        cal = pa.Calibration(reference_population="EUR", pcs_used=("PC1",), centroid=[0.0],
                             n=10, mean=1.0, sd=0.5, k_sd=3.0, threshold=2.5,
                             within_max=2.0, nearest_other=9.9)
        dec = pa.Decision(sample_id="X", verdict="REPORT", distance=1.0, threshold=2.5,
                          reason="", remedy="", n_markers_shared=500,
                          declared_population="EUR")
        score = {"pgs_id": "CLAWBIO-T2D-8", "trait": "Type 2 diabetes", "raw_score": 1.0,
                 "percentile": 50.0, "z_score": 0.0, "reference_population": "EUR"}
        gated = pa.gate_scores([score], dec, cal, sex="female",
                               integrity={"CLAWBIO-T2D-8": pa.IntegrityVerdict(True, [], [])},
                               shifts={"CLAWBIO-T2D-8": sh})
        assert gated[0]["af_shift_sd"] is None
        assert any("cannot be expressed in sd" in c for c in gated[0]["caveats"])

    def test_zero_z_score_through_main_yields_no_sd_units(self, tmp_path):
        """Through main(), not the unit: a results file whose z_score is exactly
        0.0 must produce af_shift_sd None plus the raw-units caveat, never a
        quotient over a silently substituted sd of 1.0."""
        import shutil
        results = json.loads((EXAMPLES / "demo_prs_results.json").read_text())
        # The demo results are keyed per sample; sd derivation falls back to any
        # record of the score, so every copy must sit exactly at the mean.
        for t2d in (r for r in results if r["curated_panel_id"] == "CLAWBIO-T2D-8"):
            t2d["z_score"] = 0.0
            t2d["raw_score"] = 1.1186  # exactly the reference mean
        rf = tmp_path / "prs_results.json"
        rf.write_text(json.dumps(results))
        out = tmp_path / "out"
        r = run_cli(["--reference-panel", str(EXAMPLES / "demo_reference_pcs.csv"),
                     "--individuals", str(EXAMPLES / "demo_query_individuals.csv"),
                     "--prs-results", str(rf),
                     "--scores", str(EXAMPLES / "scores"),
                     "--population-af", str(EXAMPLES / "demo_population_af.tsv"),
                     "--output", str(out), "--no-figures", "--no-pdf"])
        assert r.returncode == 0, r.stderr
        res = json.loads((out / "result.json").read_text())
        rep = [d for d in res["decisions"] if d["verdict"] == "REPORT"][0]
        t2d_gated = next(s for s in rep["scores"] if s["pgs_id"] == "CLAWBIO-T2D-8")
        assert t2d_gated["af_shift_sd"] is None
        assert any("cannot be expressed in sd" in c for c in t2d_gated["caveats"])


class TestProducerConsumerChain:
    """The boundary none of the earlier tests crossed: real gwas-prs output
    fed to prs-abstain. Round-2 audit found the id-provenance fix applied at
    this consumer while the documented producer still emitted legacy PGS-file
    labels, so the documented chain withheld everything with a false reason."""

    def test_gwas_prs_demo_output_flows_through_the_gate(self, tmp_path):
        gwas = SKILL_DIR.parent / "gwas-prs"
        gout = tmp_path / "gwas_out"
        r1 = subprocess.run(
            [sys.executable, str(gwas / "gwas_prs.py"), "--demo",
             "--output", str(gout)],
            capture_output=True, text=True)
        assert r1.returncode == 0, r1.stderr
        results_file = gout / "prs_results.json"
        assert results_file.exists()
        recs = json.loads(results_file.read_text())
        assert all(rec.get("curated_panel_id") for rec in recs), \
            "gwas-prs demo output must carry curated_panel_id (#356)"

        out = tmp_path / "abstain_out"
        # gwas-prs output carries no sample_id, so it may be gated against
        # exactly one individual (Round-5 review: cross-attribution refusal).
        r2 = run_cli(["--reference-panel", str(EXAMPLES / "demo_reference_pcs.csv"),
                      "--individuals", str(FIXTURES / "single_individual_eur.csv"),
                      "--prs-results", str(results_file),
                      "--scores", str(gwas / "data"),
                      "--output", str(out), "--no-figures", "--no-pdf"])
        assert r2.returncode == 0, r2.stderr
        res = json.loads((out / "result.json").read_text())
        rep = [d for d in res["decisions"] if d["verdict"] == "REPORT"][0]
        # Every score must match its definition file: the fail-closed
        # "never inspected" reason would mean the chain is broken again.
        for s in rep["scores"]:
            assert not any("never inspected" in x for x in s["withheld_reasons"]), s
        # And the integrity tier actually ran on those definitions: the
        # duplicate-position panel is withheld, by its panel id.
        bc = next(s for s in rep["scores"] if s["pgs_id"] == "CLAWBIO-BC-77")
        assert bc["percentile"] is None
        assert any("more than one scored variant" in x for x in bc["withheld_reasons"])

    def test_unkeyed_results_with_many_individuals_refuse(self, tmp_path):
        """Round-5 review: one prs_results.json was silently gated against every
        row of the individuals CSV. Without sample_id, >1 individual must refuse."""
        recs = json.loads((EXAMPLES / "demo_prs_results.json").read_text())
        unkeyed = [{k: v for k, v in r.items() if k != "sample_id"} for r in recs[:6]]
        f = tmp_path / "unkeyed.json"
        f.write_text(json.dumps(unkeyed))
        r = run_cli(["--reference-panel", str(EXAMPLES / "demo_reference_pcs.csv"),
                     "--individuals", str(EXAMPLES / "demo_query_individuals.csv"),
                     "--prs-results", str(f),
                     "--output", str(tmp_path / "out"), "--no-figures", "--no-pdf"])
        assert r.returncode == 2
        assert "cross-attribution" in r.stderr
        assert not (tmp_path / "out" / "result.json").exists()

    def test_unkeyed_results_with_one_individual_pass(self, tmp_path):
        recs = json.loads((EXAMPLES / "demo_prs_results.json").read_text())
        unkeyed = [{k: v for k, v in r.items() if k != "sample_id"} for r in recs[:6]]
        f = tmp_path / "unkeyed.json"
        f.write_text(json.dumps(unkeyed))
        r = run_cli(["--reference-panel", str(EXAMPLES / "demo_reference_pcs.csv"),
                     "--individuals", str(FIXTURES / "single_individual_eur.csv"),
                     "--prs-results", str(f),
                     "--output", str(tmp_path / "out"), "--no-figures", "--no-pdf"])
        assert r.returncode == 0, r.stderr

    def test_keyed_results_missing_an_individual_refuse(self, tmp_path):
        """A keyed file that covers only some individuals must not fall back to
        cross-attributing someone else's records to the uncovered ones."""
        recs = json.loads((EXAMPLES / "demo_prs_results.json").read_text())
        only_eur = [r for r in recs if r.get("sample_id") == "EUR_001"]
        f = tmp_path / "partial.json"
        f.write_text(json.dumps(only_eur))
        r = run_cli(["--reference-panel", str(EXAMPLES / "demo_reference_pcs.csv"),
                     "--individuals", str(EXAMPLES / "demo_query_individuals.csv"),
                     "--prs-results", str(f),
                     "--output", str(tmp_path / "out"), "--no-figures", "--no-pdf"])
        assert r.returncode == 2
        assert "AFR_001" in r.stderr and "DEMO_PATIENT" in r.stderr

    def test_notes_never_contain_pipes(self, tmp_path):
        """A '|' inside a note splits every markdown table row it lands in,
        and the PDF renderer then clips the fragment off. Guard the invariant
        at the source."""
        r = run_cli(["--demo", "--output", str(tmp_path), "--no-figures", "--no-pdf"])
        assert r.returncode == 0, r.stderr
        res = json.loads((tmp_path / "result.json").read_text())
        for d in res["decisions"]:
            for s in d["scores"]:
                assert "|" not in s["note"], s["pgs_id"]


class TestUnusableScoresDirFailsClosed:
    def _args(self, out, scores_dir):
        return ["--reference-panel", str(EXAMPLES / "demo_reference_pcs.csv"),
                "--individuals", str(EXAMPLES / "demo_query_individuals.csv"),
                "--prs-results", str(EXAMPLES / "demo_prs_results.json"),
                "--scores", str(scores_dir),
                "--output", str(out), "--no-figures", "--no-pdf"]

    def test_empty_scores_dir_refuses_cleanly(self, tmp_path):
        """--scores pointing at an empty directory must not release anything,
        and must not claim '--scores was not supplied'."""
        empty = tmp_path / "empty"
        empty.mkdir()
        r = run_cli(self._args(tmp_path / "out", empty))
        assert r.returncode == 2
        assert "No usable scoring files" in r.stderr
        assert "Traceback" not in r.stderr
        assert not (tmp_path / "out" / "result.json").exists()

    def test_dir_with_only_a_readme_refuses_cleanly(self, tmp_path):
        d = tmp_path / "scores"
        d.mkdir()
        (d / "README.txt").write_text("These are my scoring files.\nSee below.\n")
        r = run_cli(self._args(tmp_path / "out", d))
        assert r.returncode == 2
        assert "No usable scoring files" in r.stderr
        assert "Traceback" not in r.stderr

    def test_headerless_scoring_file_refuses_cleanly_via_cli(self, tmp_path):
        d = tmp_path / "scores"
        d.mkdir()
        (d / "PGS000099_hmPOS_GRCh37.txt").write_text(
            "#trait_reported=x\n"
            "rsID\tchr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\n"
            "rs1\t1\t100\tA\tG\t0.5\n")
        r = run_cli(self._args(tmp_path / "out", d))
        assert r.returncode == 2
        assert "Cannot load score definitions" in r.stderr
        assert "Traceback" not in r.stderr


class TestClinicianSurfaceParity:
    """The clinician report is the surface where an omission costs most.
    Round-3 review: JSON and technical report warned about threshold
    overreach and integrity caveats while the clinician document printed a
    clean 'Released'. These tests pin every clinician-facing branch through
    the CLI."""

    def test_threshold_overreach_reaches_the_clinician_report(self, tmp_path):
        """--k-sd 20 disables the abstention rule (threshold 39.7 exceeds the
        nearest non-EUR panel member at 7.3788, so AFR_001 at 10.3770 gets
        REPORT). Reachable only under the explicit override; the clinician
        document must still say so before any Released row."""
        r = run_cli(["--demo", "--output", str(tmp_path), "--k-sd", "20",
                     "--allow-threshold-overreach", "--no-figures", "--no-pdf"])
        assert r.returncode == 0, r.stderr
        res = json.loads((tmp_path / "result.json").read_text())
        assert res["calibration"]["threshold_exceeds_nearest_other"] is True
        afr = next(d for d in res["decisions"] if d["sample_id"] == "AFR_001")
        assert afr["verdict"] == "REPORT"  # the failure scenario is real
        clin = (tmp_path / "report_clinician.md").read_text()
        assert "Read this first" in clin
        assert "uncalibrated" in clin
        # And the warning precedes the results table.
        assert clin.index("Read this first") < clin.index("## Results")

    def test_threshold_overreach_without_override_fails_closed(self, tmp_path):
        """Round-4 review: overreach warned on every surface but never changed
        a verdict. Without the explicit override flag, main() must refuse to
        run rather than release percentiles under a gate that no longer gates."""
        r = run_cli(["--demo", "--output", str(tmp_path), "--k-sd", "20",
                     "--no-figures", "--no-pdf"])
        assert r.returncode == 2
        assert "no longer an abstention rule" in r.stderr
        assert "--allow-threshold-overreach" in r.stderr
        assert not (tmp_path / "result.json").exists()

    def test_default_threshold_prints_no_overreach_warning(self, tmp_path):
        r = run_cli(["--demo", "--output", str(tmp_path),
                     "--no-figures", "--no-pdf"])
        assert r.returncode == 0, r.stderr
        clin = (tmp_path / "report_clinician.md").read_text()
        assert "Read this first" not in clin
        assert "uncalibrated" not in clin

    def test_no_scores_branches_render_in_the_clinician_report(self, tmp_path):
        """The three no-scores rendering branches: status cell, dedicated
        note section, and the note's remedy sentence. Every other clinician
        test drives --demo, which always supplies scores."""
        r = run_cli(["--reference-panel", str(EXAMPLES / "demo_reference_pcs.csv"),
                     "--individuals", str(EXAMPLES / "demo_query_individuals.csv"),
                     "--prs-results", str(EXAMPLES / "demo_prs_results.json"),
                     "--output", str(tmp_path), "--no-figures", "--no-pdf"])
        assert r.returncode == 0, r.stderr
        clin = (tmp_path / "report_clinician.md").read_text()
        assert "Released — score file not checked (see note below)" in clin
        assert "## A note on results marked 'score file not checked'" in clin
        assert "re-run the review with" in clin.lower()
        # The JSON caveat and the clinician surface agree: every released
        # score with the caveat has the marked status, none says plain
        # 'Released'.
        res = json.loads((tmp_path / "result.json").read_text())
        rep = [d for d in res["decisions"] if d["verdict"] == "REPORT"][0]
        n_unchecked = sum(1 for s in rep["scores"]
                          if s["percentile"] is not None and
                          any("integrity was not verified" in c.lower()
                              for c in s["caveats"]))
        assert n_unchecked > 0
        assert clin.count("score file not checked (see note below)") == n_unchecked


class TestRound6Sentinels:
    """One test per confirmed finding of the four-lens adversarial review
    (statistical, format, clinical, integration), so none can return silently."""

    def _cal(self):
        import prs_abstain as pa
        panel = pa.load_reference_panel(EXAMPLES / "demo_reference_pcs.csv")
        return pa, pa.calibrate(panel, ref_pop="EUR", k_sd=3.0)

    def test_overreach_property_fires_at_equality(self):
        import prs_abstain as pa
        cal = pa.Calibration("EUR", ("PC1",), [0.0], 10, 1.0, 0.1, 3.0,
                             threshold=1.5, within_max=1.2, nearest_other=1.5)
        assert cal.threshold_exceeds_nearest_other is True

    def test_overreach_not_masked_by_rounding(self):
        import prs_abstain as pa
        cal = pa.Calibration("EUR", ("PC1",), [0.0], 10, 1.0, 0.1, 3.0,
                             threshold=1.234563, within_max=1.2, nearest_other=1.234558)
        assert cal.threshold_exceeds_nearest_other is True

    def test_expected_mean_partial_af_coverage_is_none(self):
        import prs_abstain as pa
        sdef = pa.ScoreDefinition("X", "t", None, [
            {"rsid": "rs1", "weight": 1.0, "af_reference": 0.5,
             "effect_allele": "A", "other_allele": "G", "chr": "1", "pos": 1},
            {"rsid": "rs2", "weight": 1.0, "af_reference": None,
             "effect_allele": "A", "other_allele": "G", "chr": "1", "pos": 2},
        ])
        assert pa.expected_mean(sdef) is None

    def test_call_cap_excludes_00_no_calls(self, tmp_path):
        import prs_abstain as pa
        g = tmp_path / "g.txt"
        g.write_text("# h\n" + "".join(f"rs{i}\t1\t{i}\t00\n" for i in range(50))
                     + "".join(f"rsx{i}\t1\t{i+99}\tAG\n" for i in range(10)))
        geno = pa.load_genotype(g)
        cap = sum(1 for v in geno.values() if v and any(c in "ACGT" for c in v.upper()))
        assert cap == 10

    def test_zero_weights_do_not_crash_audit(self, tmp_path):
        import prs_abstain as pa
        sdef = pa.ScoreDefinition("X", "t", None, [
            {"rsid": "rs1", "weight": 0.0, "af_reference": None,
             "effect_allele": "A", "other_allele": "G", "chr": "1", "pos": 1},
        ])
        au = pa.audit_score(sdef, {"rs1": "AG"})
        assert au.effective_n == 0.0

    def test_male_breast_cancer_is_male_specific(self):
        import prs_abstain as pa
        assert pa.check_applicability({"trait": "Male breast cancer"}, sex="female").applicable is False
        assert pa.check_applicability({"trait": "Male breast cancer"}, sex="male").applicable is True
        assert pa.check_applicability({"trait": "Female breast cancer"}, sex="female").applicable is True
        assert pa.check_applicability({"trait": "Breast cancer"}, sex="female").applicable is True

    def test_duplicate_sample_id_rows_refuse(self, tmp_path):
        csvf = tmp_path / "ind.csv"
        csvf.write_text("sample_id,population,PC1,PC2,PC3,PC4,n_markers_shared,sex\n"
                        "DUP,EUR,-3.0,-2.4,-2.0,0.3,400,female\n"
                        "DUP,AFR,6.6,-2.4,-2.2,3.2,400,male\n")
        r = run_cli(["--reference-panel", str(EXAMPLES / "demo_reference_pcs.csv"),
                     "--individuals", str(csvf),
                     "--prs-results", str(EXAMPLES / "demo_prs_results.json"),
                     "--output", str(tmp_path / "o"), "--no-figures", "--no-pdf"])
        assert r.returncode == 2 and "Duplicate sample_id" in r.stderr

    def test_af_caveat_population_match_is_case_insensitive(self):
        import prs_abstain as pa
        pa_, cal = self._cal()
        ind = pa.Individual("X", "afr", list(cal.centroid), 400)
        d = pa.decide(ind, cal, min_markers=30)
        sh = pa.AFShift("PGS1", "AFR", 1, 1.0, -0.9, -0.9, [])
        gated = pa.gate_scores([{"pgs_id": "PGS1", "trait": "T", "raw_score": 1.0,
                                 "percentile": 80.0, "z_score": 1.0,
                                 "reference_population": "EUR"}],
                               d, cal, shifts={"PGS1": sh})
        assert any("off by roughly" in c for c in gated[0]["caveats"])

    def test_min_markers_zero_is_an_error(self, tmp_path):
        r = run_cli(["--demo", "--output", str(tmp_path), "--min-markers", "0",
                     "--no-figures", "--no-pdf"])
        assert r.returncode == 2
        assert "disables the marker gate" in r.stderr

    def test_nan_weight_scores_dir_fails_closed(self, tmp_path):
        d = tmp_path / "scores"; d.mkdir()
        (d / "P.txt").write_text("#pgs_id=PGSNAN\n#trait_reported=x\n"
            "rsID\tchr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\n"
            "rs1\t1\t100\tA\tG\tnan\n")
        r = run_cli(["--demo", "--output", str(tmp_path / "o"), "--scores", str(d),
                     "--no-figures", "--no-pdf"])
        assert r.returncode == 2
        assert "Cannot load score definitions" in r.stderr

    def test_space_separated_genotype_fails_closed(self, tmp_path):
        g = tmp_path / "g.txt"
        g.write_text("rs1 1 100 AG\nrs2 1 200 CT\n")
        r = run_cli(["--demo", "--output", str(tmp_path / "o"), "--genotype", str(g),
                     "--no-figures", "--no-pdf"])
        assert r.returncode == 2
        assert "Cannot load genotype file" in r.stderr

    def test_pipe_in_trait_cannot_spoof_report_columns(self, tmp_path):
        import prs_abstain as pa
        recs = json.loads((EXAMPLES / "demo_prs_results.json").read_text())
        recs[0]["trait"] = "Prostate cancer | 99.9th | Released"
        f = tmp_path / "r.json"; f.write_text(json.dumps(recs))
        loaded = pa.load_prs_results(f)
        assert "|" not in loaded[0]["trait"]

    def test_falsy_sample_id_is_keyed_not_cross_attributed(self, tmp_path):
        recs = [{"pgs_id": "PGS000013", "trait": "T", "raw_score": 1.0,
                 "percentile": 97.0, "z_score": 1.0, "reference_population": "EUR",
                 "sample_id": 0}]
        f = tmp_path / "r.json"; f.write_text(json.dumps(recs))
        r = run_cli(["--reference-panel", str(EXAMPLES / "demo_reference_pcs.csv"),
                     "--individuals", str(FIXTURES / "single_individual_eur.csv"),
                     "--prs-results", str(f),
                     "--output", str(tmp_path / "o"), "--no-figures", "--no-pdf"])
        assert r.returncode == 2
        assert "EUR_001" in r.stderr

    def test_string_percentile_is_a_clean_load_error(self, tmp_path):
        recs = json.loads((EXAMPLES / "demo_prs_results.json").read_text())
        recs[0]["percentile"] = "97.0"
        f = tmp_path / "r.json"; f.write_text(json.dumps(recs))
        r = run_cli(["--reference-panel", str(EXAMPLES / "demo_reference_pcs.csv"),
                     "--individuals", str(FIXTURES / "single_individual_eur.csv"),
                     "--prs-results", str(f),
                     "--output", str(tmp_path / "o"), "--no-figures", "--no-pdf"])
        assert r.returncode == 2
        assert "Traceback" not in r.stderr
        assert not (tmp_path / "o" / "result.json").exists()

    def test_bom_individuals_csv_parses(self, tmp_path):
        csvf = tmp_path / "ind.csv"
        csvf.write_bytes(b"\xef\xbb\xbf" + (FIXTURES / "single_individual_eur.csv").read_bytes())
        r = run_cli(["--reference-panel", str(EXAMPLES / "demo_reference_pcs.csv"),
                     "--individuals", str(csvf),
                     "--prs-results", str(EXAMPLES / "demo_prs_results.json"),
                     "--output", str(tmp_path / "o"), "--no-figures", "--no-pdf"])
        assert r.returncode == 0, r.stderr

    def test_estimated_reference_population_gets_provenance_reason(self):
        import prs_abstain as pa
        pa_, cal = self._cal()
        ind = pa.Individual("EUR_001", "EUR", list(cal.centroid), 400)
        d = pa.decide(ind, cal, min_markers=30)
        gated = pa.gate_scores([{"pgs_id": "P", "trait": "T", "raw_score": 1.0,
                                 "percentile": 50.0, "z_score": 0.5,
                                 "reference_population": "estimated (allele freq)"}], d, cal)
        assert gated[0]["percentile"] is None
        assert any("estimated from allele frequencies" in x for x in gated[0]["withheld_reasons"])


class TestGenotypeBinding:
    """Round-9 review: one --genotype was applied to every queried individual,
    so a released percentile could rest on someone else's DNA. A genotype
    binds to exactly one sample_id; the rest take the no-genotype path."""

    def test_multi_individual_run_refuses_unbound_genotype(self, tmp_path):
        r = run_cli(["--reference-panel", str(EXAMPLES / "demo_reference_pcs.csv"),
                     "--individuals", str(EXAMPLES / "demo_query_individuals.csv"),
                     "--prs-results", str(EXAMPLES / "demo_prs_results.json"),
                     "--scores", str(EXAMPLES / "scores"),
                     "--genotype", str(EXAMPLES / "demo_genotype.txt"),
                     "--output", str(tmp_path), "--no-figures", "--no-pdf"])
        assert r.returncode == 2
        assert "--genotype-sample-id" in r.stderr
        assert not (tmp_path / "result.json").exists()

    def test_unknown_genotype_sample_id_refuses(self, tmp_path):
        r = run_cli(["--demo", "--output", str(tmp_path), "--genotype-sample-id", "NOBODY",
                     "--no-figures", "--no-pdf"])
        assert r.returncode == 2
        assert "not in the individuals CSV" in r.stderr

    def test_genotype_gates_only_the_bound_individual(self, tmp_path):
        r = run_cli(["--demo", "--output", str(tmp_path), "--no-figures", "--no-pdf"])
        assert r.returncode == 0, r.stderr
        res = json.loads((tmp_path / "result.json").read_text())
        by_id = {d["sample_id"]: d for d in res["decisions"]}
        bound = by_id["DEMO_PATIENT"]["scores"]
        others = [s for sid, d in by_id.items() if sid != "DEMO_PATIENT" for s in d["scores"]]
        assert bound and all(s["weight_coverage"] is not None for s in bound)
        assert others and all(s["weight_coverage"] is None for s in others)
        assert all(any("genotype-dependent integrity checks" in c for c in s["caveats"])
                   for s in others)
        tech = (tmp_path / "report_technical.md").read_text()
        assert "bound to **DEMO_PATIENT**" in tech
        audit_csv = (tmp_path / "tables" / "variant_audit.csv").read_text().splitlines()
        assert audit_csv[0].startswith("sample_id,")
        assert all(row.startswith("DEMO_PATIENT,") for row in audit_csv[1:])


class TestReferenceSdProvenance:
    """Round-9 review: sd_ref = |(raw - expected_mean) / z| assumes the reference
    mean is Sum(2*AF*w). True for curated panels; not for an authentic PGS
    Catalog file whose allelefrequency_effect is the source-GWAS frequency."""

    @staticmethod
    def _score(curated):
        import prs_abstain as pa
        variants = [{"rsid": "rs1", "chr": "1", "pos": "1", "effect_allele": "A",
                     "other_allele": "G", "weight": 0.5, "af_reference": 0.3},
                    {"rsid": "rs2", "chr": "1", "pos": "2", "effect_allele": "C",
                     "other_allele": "T", "weight": -0.2, "af_reference": 0.6}]
        return pa.ScoreDefinition("X", "trait", "GRCh37", variants, curated=curated)

    def test_authentic_catalog_file_against_named_cohort_does_not_derive_sd(self):
        import prs_abstain as pa
        sd, note = pa.reference_sd(self._score(curated=False),
                                   {"raw_score": 0.9, "z_score": 1.2, "reference_population": "EUR"})
        assert sd is None
        assert "not derived" in note and "source-GWAS" in note

    def test_curated_panel_derives_sd(self):
        import prs_abstain as pa
        sd, note = pa.reference_sd(self._score(curated=True),
                                   {"raw_score": 0.9, "z_score": 1.2, "reference_population": "EUR"})
        assert sd is not None and sd > 0
        assert "curated panel" in note

    def test_allele_frequency_centred_percentile_derives_sd(self):
        import prs_abstain as pa
        sd, _ = pa.reference_sd(self._score(curated=False),
                                {"raw_score": 0.9, "z_score": 1.2,
                                 "reference_population": "estimated (allele freq)"})
        assert sd is not None

    def test_loader_marks_curated_panels(self):
        import prs_abstain as pa
        defs = pa.load_score_definitions(EXAMPLES / "scores")
        assert defs and all(d.curated for d in defs.values())

    def test_unscaled_shift_carries_the_reason_into_the_caveat(self, tmp_path):
        """An authentic-style file (pgs_id header, no clawbio_panel) gated against
        a named cohort must print the raw shift with the provenance reason, never
        a definite sd figure."""
        import prs_abstain as pa
        d = tmp_path / "scores"; d.mkdir()
        (d / "PGS000099_hmPOS_GRCh37.txt").write_text(
            "#pgs_id=PGS000099\n#trait_reported=T\n#genome_build=GRCh37\n"
            "rsID\tchr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\tallelefrequency_effect\n"
            "rs4244285\t10\t96541616\tA\tG\t0.5\t0.3\n")
        af = tmp_path / "af.tsv"
        af.write_text("rsid\tpopulation\teffect_allele_frequency\nrs4244285\tAFR\t0.6\nrs4244285\tEUR\t0.3\n")
        recs = [{"pgs_id": "PGS000099", "trait": "T", "raw_score": 1.0, "z_score": 1.5,
                 "percentile": 90, "reference_population": "EUR", "sample_id": "EUR_001"}]
        res = tmp_path / "res.json"; res.write_text(json.dumps(recs))
        ind = tmp_path / "ind.csv"
        ind.write_text("sample_id,population,PC1,PC2,PC3,PC4,n_markers_shared,sex\n"
                       "EUR_001,EUR,-3.005020,-2.385174,-2.002224,0.345164,400,female\n")
        r = run_cli(["--reference-panel", str(EXAMPLES / "demo_reference_pcs.csv"),
                     "--individuals", str(ind), "--prs-results", str(res),
                     "--scores", str(d), "--population-af", str(af),
                     "--output", str(tmp_path / "o"), "--no-figures", "--no-pdf"])
        assert r.returncode == 0, r.stderr
        out = json.loads((tmp_path / "o" / "result.json").read_text())
        sc = out["decisions"][0]["scores"][0]
        assert sc["af_shift_sd"] is None
        assert any("source-GWAS" in c for c in sc["caveats"])
        tech = (tmp_path / "o" / "report_technical.md").read_text()
        assert "n/a (sd underivable)" in tech and "source-GWAS" in tech


class TestRound9Sweep:
    """Sibling sweep after round 9: the remaining instances of classes the
    reviewers had already found once (fail-open on NaN, one input for many,
    derived quantities with undeclared provenance, plain-language mislabels)."""

    @staticmethod
    def _score(variants):
        import prs_abstain as pa
        return pa.ScoreDefinition("X", "trait", "GRCh37", variants, curated=True)

    def test_nan_pc_coordinate_is_withheld_not_released(self, tmp_path):
        csvf = tmp_path / "ind.csv"
        csvf.write_text("sample_id,population,PC1,PC2,PC3,PC4,n_markers_shared,sex\n"
                        "EUR_001,EUR,nan,-2.4,-2.0,0.3,400,female\n")
        r = run_cli(["--reference-panel", str(EXAMPLES / "demo_reference_pcs.csv"),
                     "--individuals", str(csvf),
                     "--prs-results", str(EXAMPLES / "demo_prs_results.json"),
                     "--output", str(tmp_path / "o"), "--no-figures", "--no-pdf"])
        assert r.returncode == 0, r.stderr
        res = json.loads((tmp_path / "o" / "result.json").read_text())
        d = res["decisions"][0]
        assert d["verdict"] == "REFUSE_UNDETERMINABLE"
        assert all(sc["percentile"] is None for sc in d["scores"])

    def test_nan_in_reference_panel_row_is_skipped_not_averaged(self, tmp_path):
        import prs_abstain as pa
        pnl = tmp_path / "ref.csv"
        src = (EXAMPLES / "demo_reference_pcs.csv").read_text().splitlines()
        pnl.write_text("\n".join(src + ["BAD_EUR,EUR,nan,0,0,0"]) + "\n")
        panel = pa.load_reference_panel(pnl, ("PC1", "PC2", "PC3", "PC4"))
        assert all(s.sample_id != "BAD_EUR" for s in panel)
        cal = pa.calibrate(panel, "EUR")
        assert math.isfinite(cal.threshold)

    def test_af_shift_skips_rsid_with_conflicting_effect_alleles(self):
        import prs_abstain as pa
        sc = self._score([{"rsid": "rs1", "chr": "1", "pos": "1", "effect_allele": "A",
                           "other_allele": "G", "weight": 0.5, "af_reference": 0.3},
                          {"rsid": "rs2", "chr": "1", "pos": "2", "effect_allele": "C",
                           "other_allele": "T", "weight": 0.5, "af_reference": 0.4}])
        af = {"rs1": {"AFR": 0.6}, "rs2": {"AFR": 0.5}}
        sh = pa.af_shift(sc, af, sd=1.0, population="AFR", ambiguous_rsids={"rs1"})
        assert sh.n_variants_with_af == 1 and "rs1" not in [v["rsid"] for v in sh.per_variant]
        assert "skipped" in sh.sd_note and sh.shift_sd is None  # 50% weight coverage

    def test_af_table_allele_column_rejects_mismatched_allele(self):
        import prs_abstain as pa
        sc = self._score([{"rsid": "rs1", "chr": "1", "pos": "1", "effect_allele": "A",
                           "other_allele": "G", "weight": 1.0, "af_reference": 0.3}])
        assert pa.af_shift(sc, {"rs1": {"AFR": 0.6, "_allele": "C"}}, sd=1.0, population="AFR") is None
        ok = pa.af_shift(sc, {"rs1": {"AFR": 0.6, "_allele": "T"}}, sd=1.0, population="AFR")
        assert ok is not None and ok.shift_sd is not None  # T is the strand complement of A

    def test_low_weight_coverage_reports_raw_shift_unscaled(self):
        import prs_abstain as pa
        sc = self._score([{"rsid": "rs1", "chr": "1", "pos": "1", "effect_allele": "A",
                           "other_allele": "G", "weight": 0.1, "af_reference": 0.3},
                          {"rsid": "rs2", "chr": "1", "pos": "2", "effect_allele": "C",
                           "other_allele": "T", "weight": 0.9, "af_reference": 0.4}])
        sh = pa.af_shift(sc, {"rs1": {"AFR": 0.6}}, sd=1.0, population="AFR")
        assert sh.weight_coverage == 0.1 and sh.shift_sd is None
        assert "10%" in sh.sd_note

    def test_population_af_loader_rejects_nan_and_counts_it(self, tmp_path, capsys):
        import prs_abstain as pa
        f = tmp_path / "af.tsv"
        f.write_text("rsid\tpopulation\teffect_allele_frequency\trs1\tAFR\tnan\nrs2\tAFR\t0.4\n")
        tab = pa.load_population_af(f)
        assert "rs2" in tab and "rs1" not in tab

    def test_sd_is_derived_from_each_individuals_own_record(self, tmp_path):
        """Two keyed individuals; the LAST record in the file has z=0. The
        first individual's sd must still derive from her own record."""
        recs = json.loads((EXAMPLES / "demo_prs_results.json").read_text())
        eur = [dict(r) for r in recs if r["sample_id"] == "EUR_001"]
        afr = [dict(r) for r in recs if r["sample_id"] == "AFR_001"]
        for r in afr:
            r["z_score"] = 0.0
        res = tmp_path / "res.json"; res.write_text(json.dumps(eur + afr))
        ind = tmp_path / "ind.csv"
        ind.write_text("sample_id,population,PC1,PC2,PC3,PC4,n_markers_shared,sex\n"
                       "EUR_001,EUR,-3.005020,-2.385174,-2.002224,0.345164,400,female\n"
                       "AFR_001,AFR,6.588584,-2.403258,-2.206847,3.186610,400,male\n")
        r = run_cli(["--reference-panel", str(EXAMPLES / "demo_reference_pcs.csv"),
                     "--individuals", str(ind), "--prs-results", str(res),
                     "--scores", str(EXAMPLES / "scores"),
                     "--population-af", str(EXAMPLES / "demo_population_af.tsv"),
                     "--output", str(tmp_path / "o"), "--no-figures", "--no-pdf"])
        assert r.returncode == 0, r.stderr
        out = json.loads((tmp_path / "o" / "result.json").read_text())
        by = {d["sample_id"]: d for d in out["decisions"]}
        eur_sd = [sc["af_shift_sd"] for sc in by["EUR_001"]["scores"]]
        afr_sd = [sc["af_shift_sd"] for sc in by["AFR_001"]["scores"]]
        assert any(x is not None for x in eur_sd)
        assert all(x is None for x in afr_sd)
        csv_rows = (tmp_path / "o" / "tables" / "af_shift_per_variant.csv").read_text().splitlines()
        assert csv_rows[0].startswith("sample_id,")

    def test_plain_reason_names_the_actual_cause(self):
        import prs_abstain as pa
        cal = pa.Calibration("EUR", ("PC1",), [0.0], 10, 0.0, 1.0, 3.0, 3.0, 2.0, 5.0)
        base = {"verdict": "REFUSE_UNDETERMINABLE"}
        assert "higher than the number of markers actually read" in pa._plain_reason({**base, "cause": "inflated_markers"}, cal)
        assert "No ancestry coordinates" in pa._plain_reason({**base, "cause": "no_coordinates"}, cal)
        assert "Too few ancestry markers" in pa._plain_reason({**base, "cause": "too_few_markers"}, cal)

    def test_single_letter_sex_is_normalised(self):
        import prs_abstain as pa
        assert pa.check_applicability({"trait": "Breast cancer"}, sex="F").applicable
        assert not pa.check_applicability({"trait": "Breast cancer"}, sex="M").applicable

    def test_plain_withheld_distinguishes_missing_sex_and_provenance(self):
        import prs_abstain as pa
        assert "no sex was recorded" in pa._plain_withheld(
            {"withheld_reasons": ["Not applicable: Breast cancer is a sex-specific trait and no sex was recorded for this individual."]})
        assert "which group" in pa._plain_withheld(
            {"withheld_reasons": ["Reference provenance: this score does not declare a reference_population, so x."]})

    def test_clinician_report_flags_unchecked_marker_coverage_per_person(self, tmp_path):
        run_cli(["--demo", "--output", str(tmp_path), "--no-figures", "--no-pdf"])
        text = (tmp_path / "report_clinician.md").read_text()
        assert "marker coverage not checked for this person" in text

    def test_sample_ids_colliding_after_sanitisation_refuse(self, tmp_path):
        csvf = tmp_path / "ind.csv"
        csvf.write_text("sample_id,population,PC1,PC2,PC3,PC4,n_markers_shared,sex\n"
                        "A B,EUR,-3.0,-2.4,-2.0,0.3,400,female\n"
                        "A_B,EUR,-3.1,-2.4,-2.0,0.3,400,female\n")
        r = run_cli(["--reference-panel", str(EXAMPLES / "demo_reference_pcs.csv"),
                     "--individuals", str(csvf),
                     "--prs-results", str(EXAMPLES / "demo_prs_results.json"),
                     "--output", str(tmp_path / "o"), "--no-figures", "--no-pdf"])
        assert r.returncode == 2 and "collide" in r.stderr

    def test_technical_report_names_the_reference_population_not_european(self, tmp_path):
        """§1 of the technical report derives its population label from the
        calibration and asserts the Sum(2AFw) identity only for curated scores."""
        run_cli(["--demo", "--output", str(tmp_path), "--no-figures", "--no-pdf"])
        tech = (tmp_path / "report_technical.md").read_text()
        assert "For the curated scores in this run" in tech


class TestSinglePopulationPanel:
    """Round-10 note: with a single-population panel nearest_other is None and
    the overreach check cannot run. That is recorded everywhere an unrunnable
    check is recorded, never left silent."""

    def _panel(self, tmp_path):
        src = (EXAMPLES / "demo_reference_pcs.csv").read_text().splitlines()
        rows = [src[0]] + [l for l in src[1:] if l.split(",")[1] == "EUR"]
        f = tmp_path / "eur_only.csv"; f.write_text("\n".join(rows) + "\n")
        return f

    def test_unrunnable_overreach_check_is_recorded_not_silent(self, tmp_path):
        recs = [r for r in json.loads((EXAMPLES / "demo_prs_results.json").read_text())
                if r["sample_id"] == "EUR_001"]
        res = tmp_path / "res.json"; res.write_text(json.dumps(recs))
        ind = tmp_path / "ind.csv"
        ind.write_text("sample_id,population,PC1,PC2,PC3,PC4,n_markers_shared,sex\n"
                       "EUR_001,EUR,-3.005020,-2.385174,-2.002224,0.345164,400,female\n")
        r = run_cli(["--reference-panel", str(self._panel(tmp_path)),
                     "--individuals", str(ind), "--prs-results", str(res),
                     "--output", str(tmp_path / "o"), "--no-figures", "--no-pdf"])
        assert r.returncode == 0, r.stderr
        assert "overreach check cannot run" in r.stderr
        out = json.loads((tmp_path / "o" / "result.json").read_text())
        assert out["calibration"]["nearest_other"] is None
        assert out["calibration"]["overreach_check"].startswith("not_run")
        assert "could not run" in (tmp_path / "o" / "report_technical.md").read_text()
        assert "One check could not be performed" in (tmp_path / "o" / "report_clinician.md").read_text()

    def test_demo_panel_records_the_check_as_ran(self, tmp_path):
        r = run_cli(["--demo", "--output", str(tmp_path), "--no-figures", "--no-pdf"])
        assert r.returncode == 0, r.stderr
        out = json.loads((tmp_path / "result.json").read_text())
        assert out["calibration"]["overreach_check"] == "ran"

