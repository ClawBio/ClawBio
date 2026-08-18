"""Tests for prs-abstain. Written before the implementation (red/green TDD)."""
from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

import pytest

SKILL_DIR = Path(__file__).resolve().parents[1]
SCRIPT = SKILL_DIR / "prs_abstain.py"
EXAMPLES = SKILL_DIR / "examples"

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
        scores = pa.load_prs_results(EXAMPLES / "demo_prs_results.json")
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
        """k large enough to admit non-reference individuals must warn loudly."""
        r = run_cli(["--demo", "--output", str(tmp_path), "--k-sd", "20"])
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
        a = pa.audit_score(self._defs()["PGS000013"], gt)
        assert a.weight_coverage == pytest.approx(1.0)

    def test_missing_variants_reduce_weight_coverage(self):
        import prs_abstain as pa
        gt = pa.load_genotype(EXAMPLES / "demo_genotype.txt")
        thin = {k: v for i, (k, v) in enumerate(gt.items()) if i % 2 == 0}
        a = pa.audit_score(self._defs()["PGS000013"], thin)
        assert a.weight_coverage < 1.0
        assert a.weight_at_risk > 0

    def test_concentration_detects_fragile_score(self):
        import prs_abstain as pa
        gt = pa.load_genotype(EXAMPLES / "demo_genotype.txt")
        t2d = pa.audit_score(self._defs()["PGS000013"], gt)   # 8 variants
        bmi = pa.audit_score(self._defs()["PGS000039"], gt)   # 97 variants
        assert t2d.effective_n < 10 < bmi.effective_n
        assert t2d.top1_share > bmi.top1_share

    def test_palindromic_variants_are_counted(self):
        import prs_abstain as pa
        gt = pa.load_genotype(EXAMPLES / "demo_genotype.txt")
        a = pa.audit_score(self._defs()["PGS000013"], gt)
        assert a.palindromic_n >= 1

    def test_low_weight_coverage_refuses(self):
        import prs_abstain as pa
        gt = pa.load_genotype(EXAMPLES / "demo_genotype.txt")
        thin = {k: v for i, (k, v) in enumerate(gt.items()) if i % 4 == 0}
        a = pa.audit_score(self._defs()["PGS000013"], thin)
        assert pa.integrity_verdict(a, min_weight_coverage=0.90).passed is False


class TestAFShift:
    def test_reference_mean_equals_af_expectation(self):
        """The curated EUR mean is Sum 2*AF*w. This is the whole mechanism."""
        import prs_abstain as pa
        defs = pa.load_score_definitions(Path(__file__).resolve().parents[1] / "examples" / "scores")
        assert pa.expected_mean(defs["PGS000013"]) == pytest.approx(1.12, abs=0.01)
        assert pa.expected_mean(defs["PGS000004"]) == pytest.approx(2.84, abs=0.01)

    def test_af_shift_quantifies_percentile_error(self):
        import prs_abstain as pa
        defs = pa.load_score_definitions(Path(__file__).resolve().parents[1] / "examples" / "scores")
        af = pa.load_population_af(EXAMPLES / "demo_population_af.tsv")
        sh = pa.af_shift(defs["PGS000013"], af, sd=0.30)
        assert sh is not None
        assert sh.n_variants_with_af > 0
        assert isinstance(sh.shift_sd, float)
        assert len(sh.per_variant) == sh.n_variants_with_af

    def test_af_shift_returns_none_without_af_data(self):
        import prs_abstain as pa
        defs = pa.load_score_definitions(Path(__file__).resolve().parents[1] / "examples" / "scores")
        assert pa.af_shift(defs["PGS000013"], {}, sd=0.30) is None

    def test_per_variant_contributions_sum_to_total_shift(self):
        import prs_abstain as pa
        defs = pa.load_score_definitions(Path(__file__).resolve().parents[1] / "examples" / "scores")
        af = pa.load_population_af(EXAMPLES / "demo_population_af.tsv")
        sh = pa.af_shift(defs["PGS000013"], af, sd=0.30)
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
        ld = pa.ld_audit(self._defs()["PGS000013"], window_kb=250)
        assert ld.n_clusters_multi == 1
        members = [set(c["rsids"]) for c in ld.clusters if len(c["rsids"]) > 1]
        assert {"rs7903146", "rs12255372"} in members

    def test_effective_n_falls_when_ld_accounted_for(self):
        import prs_abstain as pa
        ld = pa.ld_audit(self._defs()["PGS000013"], window_kb=250)
        assert ld.effective_n_ld < ld.effective_n_independent
        assert ld.effective_n_ld == pytest.approx(3.55, abs=0.1)

    def test_clustered_weight_share_reported(self):
        import prs_abstain as pa
        ld = pa.ld_audit(self._defs()["PGS000013"], window_kb=250)
        assert ld.clustered_weight_share == pytest.approx(0.48, abs=0.02)

    def test_duplicate_positions_flagged_as_data_error(self):
        import prs_abstain as pa
        ld = pa.ld_audit(self._defs()["PGS000001"], window_kb=250)
        assert ld.duplicate_positions
        assert any("rs11552449" in d["rsids"] for d in ld.duplicate_positions)

    def test_clean_score_has_no_duplicates(self):
        import prs_abstain as pa
        assert not pa.ld_audit(self._defs()["PGS000013"], window_kb=250).duplicate_positions

    def test_duplicate_positions_block_the_score(self):
        import prs_abstain as pa
        gt = pa.load_genotype(EXAMPLES / "demo_genotype.txt")
        defs = self._defs()
        au = pa.audit_score(defs["PGS000001"], gt)
        ld = pa.ld_audit(defs["PGS000001"], window_kb=250)
        assert pa.integrity_verdict(au, ld=ld).passed is False

    def test_ld_warning_uses_corrected_effective_n(self):
        import prs_abstain as pa
        gt = pa.load_genotype(EXAMPLES / "demo_genotype.txt")
        defs = self._defs()
        au = pa.audit_score(defs["PGS000013"], gt)
        ld = pa.ld_audit(defs["PGS000013"], window_kb=250)
        v = pa.integrity_verdict(au, ld=ld, min_effective_n=10.0)
        assert any("3.5" in w or "3.6" in w for w in v.warnings)

    def test_window_size_changes_clustering(self):
        import prs_abstain as pa
        d = self._defs()["PGS000013"]
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
        """PDF is a bonus artefact; its absence must not break the run."""
        import prs_abstain as pa
        assert hasattr(pa, "_write_pdfs")

    def test_no_pdf_flag_skips_generation(self, tmp_path):
        r = run_cli(["--demo", "--output", str(tmp_path), "--no-pdf"])
        assert r.returncode == 0, r.stderr
        assert not (tmp_path / "report_clinician.pdf").exists()
