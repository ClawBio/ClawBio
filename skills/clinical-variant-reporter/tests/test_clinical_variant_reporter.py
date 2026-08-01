"""Tests for clinical-variant-reporter ACMG classification skill.

Validates the ACMG/AMP 2015 combining rules, individual criteria evaluation,
secondary findings screening, and demo mode end-to-end against a curated panel
of GIAB HG001 / ClinVar variants with known expected classifications.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

SKILL_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(SKILL_DIR))

from acmg_engine import (
    ACMG_SF_V32_GENES,
    ClassifiedVariant,
    EvidenceCriterion,
    VariantEvidence,
    classify,
    classify_variant,
    evaluate_criteria,
    is_secondary_finding_gene,
)
from clinical_variant_reporter import (
    VEP_REST_URL,
    annotate_variants_vep,
    build_evidence_from_cache,
    generate_report,
    load_demo_evidence_cache,
    parse_vcf,
    run_classification,
)


# ---------------------------------------------------------------------------
# Unit tests — combining rules
# ---------------------------------------------------------------------------
class TestCombiningRules:
    """Verify the ACMG combining rules produce correct classifications."""

    @staticmethod
    def _make_criteria(codes: list[tuple[str, str, str]]) -> list[EvidenceCriterion]:
        return [
            EvidenceCriterion(
                code=code, triggered=True, strength=strength,
                direction=direction, source="test", detail="test",
            )
            for code, strength, direction in codes
        ]

    def test_ba1_standalone_benign(self):
        criteria = self._make_criteria([("BA1", "stand_alone", "benign")])
        assert classify(criteria) == "Benign"

    def test_ba1_overrides_pathogenic(self):
        criteria = self._make_criteria([
            ("BA1", "stand_alone", "benign"),
            ("PVS1", "very_strong", "pathogenic"),
            ("PS1", "strong", "pathogenic"),
        ])
        assert classify(criteria) == "Benign"

    def test_two_strong_benign(self):
        criteria = self._make_criteria([
            ("BS1", "strong", "benign"),
            ("BS2", "strong", "benign"),
        ])
        assert classify(criteria) == "Benign"

    def test_pvs1_plus_ps_is_pathogenic(self):
        criteria = self._make_criteria([
            ("PVS1", "very_strong", "pathogenic"),
            ("PS1", "strong", "pathogenic"),
        ])
        assert classify(criteria) == "Pathogenic"

    def test_pvs1_plus_two_pm_is_pathogenic(self):
        criteria = self._make_criteria([
            ("PVS1", "very_strong", "pathogenic"),
            ("PM1", "moderate", "pathogenic"),
            ("PM2", "moderate", "pathogenic"),
        ])
        assert classify(criteria) == "Pathogenic"

    def test_pvs1_plus_pm_plus_pp_is_pathogenic(self):
        criteria = self._make_criteria([
            ("PVS1", "very_strong", "pathogenic"),
            ("PM2", "moderate", "pathogenic"),
            ("PP3", "supporting", "pathogenic"),
        ])
        assert classify(criteria) == "Pathogenic"

    def test_two_strong_pathogenic(self):
        criteria = self._make_criteria([
            ("PS1", "strong", "pathogenic"),
            ("PS4", "strong", "pathogenic"),
        ])
        assert classify(criteria) == "Pathogenic"

    def test_pvs1_plus_one_pm_is_likely_pathogenic(self):
        criteria = self._make_criteria([
            ("PVS1", "very_strong", "pathogenic"),
            ("PM2", "moderate", "pathogenic"),
        ])
        assert classify(criteria) == "Likely Pathogenic"

    def test_ps_plus_pm_is_likely_pathogenic(self):
        criteria = self._make_criteria([
            ("PS1", "strong", "pathogenic"),
            ("PM2", "moderate", "pathogenic"),
        ])
        assert classify(criteria) == "Likely Pathogenic"

    def test_three_moderate_is_likely_pathogenic(self):
        criteria = self._make_criteria([
            ("PM1", "moderate", "pathogenic"),
            ("PM2", "moderate", "pathogenic"),
            ("PM4", "moderate", "pathogenic"),
        ])
        assert classify(criteria) == "Likely Pathogenic"

    def test_bs_plus_bp_is_likely_benign(self):
        criteria = self._make_criteria([
            ("BS1", "strong", "benign"),
            ("BP4", "supporting", "benign"),
        ])
        assert classify(criteria) == "Likely Benign"

    def test_two_supporting_benign_is_likely_benign(self):
        criteria = self._make_criteria([
            ("BP4", "supporting", "benign"),
            ("BP7", "supporting", "benign"),
        ])
        assert classify(criteria) == "Likely Benign"

    def test_single_moderate_is_vus(self):
        criteria = self._make_criteria([
            ("PM2", "moderate", "pathogenic"),
        ])
        assert classify(criteria) == "Uncertain Significance"

    def test_conflicting_evidence_is_vus(self):
        criteria = self._make_criteria([
            ("PVS1", "very_strong", "pathogenic"),
            ("BS1", "strong", "benign"),
        ])
        assert classify(criteria) == "Uncertain Significance"

    def test_empty_criteria_is_vus(self):
        assert classify([]) == "Uncertain Significance"


# ---------------------------------------------------------------------------
# Unit tests — criteria evaluation
# ---------------------------------------------------------------------------
class TestCriteriaEvaluation:
    def test_ba1_triggered_when_af_above_5_percent(self):
        ev = VariantEvidence(chrom="chr1", pos=100, ref="A", alt="G", gnomad_af=0.12)
        criteria = evaluate_criteria(ev)
        ba1 = next(c for c in criteria if c.code == "BA1")
        assert ba1.triggered is True

    def test_ba1_not_triggered_when_af_below_5_percent(self):
        ev = VariantEvidence(chrom="chr1", pos=100, ref="A", alt="G", gnomad_af=0.03)
        criteria = evaluate_criteria(ev)
        ba1 = next(c for c in criteria if c.code == "BA1")
        assert ba1.triggered is False

    def test_pvs1_triggered_for_frameshift(self):
        ev = VariantEvidence(
            chrom="chr17", pos=100, ref="AG", alt="A",
            consequence="frameshift_variant", is_lof=True,
        )
        criteria = evaluate_criteria(ev)
        pvs1 = next(c for c in criteria if c.code == "PVS1")
        assert pvs1.triggered is True

    def test_pm2_triggered_for_absent_gnomad(self):
        ev = VariantEvidence(chrom="chr1", pos=100, ref="A", alt="G", gnomad_af=None)
        criteria = evaluate_criteria(ev)
        pm2 = next(c for c in criteria if c.code == "PM2")
        assert pm2.triggered is True

    def test_pp3_triggered_for_deleterious_missense(self):
        ev = VariantEvidence(
            chrom="chr1", pos=100, ref="A", alt="G",
            is_missense=True, cadd_phred=30.0,
            sift_prediction="deleterious",
            polyphen_prediction="probably_damaging",
        )
        criteria = evaluate_criteria(ev)
        pp3 = next(c for c in criteria if c.code == "PP3")
        assert pp3.triggered is True

    def test_bp7_triggered_for_synonymous_no_splice(self):
        ev = VariantEvidence(
            chrom="chr1", pos=100, ref="A", alt="G",
            consequence="synonymous_variant", is_synonymous=True,
            spliceai_max_delta=0.01,
        )
        criteria = evaluate_criteria(ev)
        bp7 = next(c for c in criteria if c.code == "BP7")
        assert bp7.triggered is True

    def test_pp5_triggered_for_pathogenic_clinvar(self):
        ev = VariantEvidence(
            chrom="chr1", pos=100, ref="A", alt="G",
            clinvar_significance="Pathogenic", clinvar_review_stars=3,
        )
        criteria = evaluate_criteria(ev)
        pp5 = next(c for c in criteria if c.code == "PP5")
        assert pp5.triggered is True

    def test_bp6_triggered_for_benign_clinvar(self):
        ev = VariantEvidence(
            chrom="chr1", pos=100, ref="A", alt="G",
            clinvar_significance="Benign", clinvar_review_stars=2,
        )
        criteria = evaluate_criteria(ev)
        bp6 = next(c for c in criteria if c.code == "BP6")
        assert bp6.triggered is True


# ---------------------------------------------------------------------------
# Unit tests — SF v3.2 screening
# ---------------------------------------------------------------------------
class TestSecondaryFindings:
    def test_sf_gene_count(self):
        assert len(ACMG_SF_V32_GENES) == 81

    def test_brca1_is_sf_gene(self):
        assert is_secondary_finding_gene("BRCA1") is True

    def test_calm1_is_sf_gene(self):
        assert is_secondary_finding_gene("CALM1") is True

    def test_cdh1_is_sf_gene(self):
        assert is_secondary_finding_gene("CDH1") is True

    def test_dpyd_is_not_sf_gene(self):
        assert is_secondary_finding_gene("DPYD") is False

    def test_random_gene_is_not_sf_gene(self):
        assert is_secondary_finding_gene("FAKE_GENE") is False


# ---------------------------------------------------------------------------
# Integration tests — demo mode
# ---------------------------------------------------------------------------
class TestDemoMode:
    @pytest.fixture
    def demo_vcf_path(self):
        return SKILL_DIR / "example_data" / "giab_acmg_panel.vcf"

    @pytest.fixture
    def demo_cache(self):
        return load_demo_evidence_cache()

    def test_demo_vcf_parseable(self, demo_vcf_path):
        records = parse_vcf(demo_vcf_path)
        assert len(records) == 20

    def test_demo_cache_has_all_variants(self, demo_vcf_path, demo_cache):
        records = parse_vcf(demo_vcf_path)
        for rec in records:
            key = f"{rec.chrom}:{rec.pos}:{rec.ref}:{rec.alt}"
            assert key in demo_cache, f"Missing cache entry for {key}"

    def test_demo_classification_counts(self, demo_vcf_path):
        records = parse_vcf(demo_vcf_path)
        classified, _ = run_classification(records, demo=True)
        assert len(classified) == 20

        counts: dict[str, int] = {}
        for cv in classified:
            counts[cv.classification] = counts.get(cv.classification, 0) + 1

        assert counts.get("Pathogenic", 0) == 4
        assert counts.get("Likely Pathogenic", 0) == 3
        assert counts.get("Uncertain Significance", 0) == 4
        assert counts.get("Benign", 0) == 3
        assert counts.get("Likely Benign", 0) == 6

    def test_demo_expected_classifications(self, demo_vcf_path):
        """Validate each demo variant against its EXPECTED INFO field."""
        records = parse_vcf(demo_vcf_path)
        classified, _ = run_classification(records, demo=True)

        expected_map = {
            "Pathogenic": "Pathogenic",
            "Likely_Pathogenic": "Likely Pathogenic",
            "VUS": "Uncertain Significance",
            "Likely_Benign": "Likely Benign",
            "Benign": "Benign",
        }

        for rec, cv in zip(records, classified):
            expected_raw = rec.info.get("EXPECTED", "")
            expected_class = expected_map.get(expected_raw, expected_raw)
            assert cv.classification == expected_class, (
                f"Variant {rec.chrom}:{rec.pos} {rec.info.get('GENE', '')} "
                f"expected {expected_class} but got {cv.classification} "
                f"(triggered: {cv.triggered_codes})"
            )

    def test_demo_secondary_findings_screening(self, demo_vcf_path):
        records = parse_vcf(demo_vcf_path)
        classified, _ = run_classification(records, demo=True)

        sf_variants = [cv for cv in classified if cv.is_secondary_finding]
        non_sf = [cv for cv in classified if not cv.is_secondary_finding]

        assert len(sf_variants) >= 17
        assert any(cv.evidence.gene == "DPYD" for cv in non_sf)

    def test_demo_report_generation(self, demo_vcf_path, tmp_path):
        records = parse_vcf(demo_vcf_path)
        classified, _ = run_classification(records, demo=True)
        generate_report(classified, tmp_path, demo=True)

        assert (tmp_path / "report.md").exists()
        assert (tmp_path / "result.json").exists()
        assert (tmp_path / "tables" / "acmg_classifications.tsv").exists()
        assert (tmp_path / "tables" / "secondary_findings.tsv").exists()
        assert (tmp_path / "reproducibility" / "commands.sh").exists()
        assert (tmp_path / "reproducibility" / "database_versions.json").exists()

        report_text = (tmp_path / "report.md").read_text()
        assert "ACMG" in report_text
        assert "Pathogenic" in report_text
        assert "ClawBio is a research" in report_text

        result = json.loads((tmp_path / "result.json").read_text())
        assert result["total_variants"] == 20
        assert result["framework"] == "ACMG/AMP 2015 (Richards et al., PMID 25741868)"

    def test_demo_transcripts_are_versioned(self, demo_vcf_path):
        """HGVS v21.1 requires versioned transcript accessions (e.g. ENST00000357654.9)."""
        records = parse_vcf(demo_vcf_path)
        classified, _ = run_classification(records, demo=True)
        import re
        versioned_pattern = re.compile(r"^ENST\d+\.\d+$")
        for cv in classified:
            transcript = cv.evidence.transcript
            if transcript:
                assert versioned_pattern.match(transcript), (
                    f"Unversioned transcript {transcript} for {cv.evidence.gene} "
                    f"at {cv.evidence.chrom}:{cv.evidence.pos}"
                )

    def test_gene_filter(self, demo_vcf_path):
        records = parse_vcf(demo_vcf_path)
        classified, _ = run_classification(records, demo=True, gene_filter={"BRCA1", "TP53"})
        genes = {cv.evidence.gene for cv in classified}
        assert genes <= {"BRCA1", "TP53"}
        assert len(classified) >= 2


# ---------------------------------------------------------------------------
# Unit tests — VEP live-path parameters
# ---------------------------------------------------------------------------
class TestVepLivePath:
    """Verify that annotate_variants_vep sends transcript_version=1 to Ensembl."""

    def test_transcript_version_param_sent(self, monkeypatch):
        from unittest.mock import MagicMock
        from clinical_variant_reporter import VcfRecord

        captured_kwargs: list[dict] = []

        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.json.return_value = [
            {
                "input": "1 100 100 A/G 1",
                "most_severe_consequence": "missense_variant",
                "transcript_consequences": [
                    {
                        "gene_symbol": "FAKEGENE",
                        "impact": "MODERATE",
                        "consequence_terms": ["missense_variant"],
                        "transcript_id": "ENST00000000001.3",
                    }
                ],
            }
        ]

        def mock_post(*args, **kwargs):
            captured_kwargs.append(kwargs)
            return mock_response

        import requests
        monkeypatch.setattr(requests, "post", mock_post)
        monkeypatch.setattr(
            "clinical_variant_reporter.VEP_RATE_LIMIT_SECONDS", 0,
        )

        records = [VcfRecord(chrom="1", pos=100, id=".", ref="A", alt="G",
                             qual=".", filt="PASS", info={})]
        annotate_variants_vep(records)

        assert len(captured_kwargs) == 1
        params = captured_kwargs[0].get("params", {})
        assert params.get("transcript_version") == 1, (
            f"Expected transcript_version=1 in VEP params, got: {params}"
        )


# ---------------------------------------------------------------------------
# Regression — VEP clin_sig_allele may be a string, not a dict
# ---------------------------------------------------------------------------
class TestClinSigAlleleTypeGuard:
    """VEP returns colocated_variants[].clin_sig_allele as a dict for some
    variants and a string (e.g. "T:pathogenic") for others. Extraction must not
    crash on the string form (previously AttributeError: 'str' has no 'get')."""

    def _record(self):
        from clinical_variant_reporter import VcfRecord
        return VcfRecord(chrom="13", pos=32339267, id="rs886040553",
                         ref="A", alt="T", qual=".", filt="PASS", info={})

    def _vep(self, clin_sig_allele):
        return {
            "most_severe_consequence": "stop_gained",
            "transcript_consequences": [
                {"gene_symbol": "BRCA2", "impact": "HIGH",
                 "consequence_terms": ["stop_gained"],
                 "transcript_id": "ENST00000544455.6"}
            ],
            "colocated_variants": [
                {"clin_sig": ["pathogenic"], "clin_sig_allele": clin_sig_allele}
            ],
        }

    def test_string_clin_sig_allele_does_not_crash(self):
        from clinical_variant_reporter import _extract_evidence_from_vep
        ev = _extract_evidence_from_vep(self._vep("T:pathogenic"), self._record())
        assert ev.gene == "BRCA2"
        assert ev.clinvar_review_stars == 0  # defaults gracefully when not a dict

    def test_dict_clin_sig_allele_still_reads_stars(self):
        from clinical_variant_reporter import _extract_evidence_from_vep
        ev = _extract_evidence_from_vep(
            self._vep({"review_status_stars": 3}), self._record())
        assert ev.clinvar_review_stars == 3


# ---------------------------------------------------------------------------
# Regression — Data Sources table hardcoded stale ClinVar/gnomAD versions
# ---------------------------------------------------------------------------
class TestDataSourcesVersioning:
    """The Data Sources section must reflect the ClinVar/dbSNP version Ensembl
    actually served for the run, keep gnomAD's own hardcoded version (Ensembl's
    variation endpoint does not list gnomAD), and label demo mode as a cache
    rather than implying a live query (issue #303)."""

    def _report_text(self, tmp_path, demo, source_versions=None):
        from clinical_variant_reporter import generate_report
        generate_report([], tmp_path, demo=demo, source_versions=source_versions)
        return (tmp_path / "report.md").read_text()

    def test_live_mode_reports_real_clinvar_version_and_keeps_gnomad(self, tmp_path):
        # source_versions is threaded in from annotation time (see
        # TestAnnotateCapturesSourceVersions below), not fetched by generate_report.
        text = self._report_text(
            tmp_path, demo=False,
            source_versions={"clinvar": "09/2025", "dbsnp": "156", "omim": "09/2025"},
        )
        assert "ClinVar | 09/2025 (via Ensembl VEP colocated variants)" in text
        assert "gnomAD | v4.1 (via Ensembl VEP colocated variants)" in text
        assert "2025-03-01 release" not in text

    def test_live_mode_suppresses_claim_when_nothing_annotated(self, tmp_path):
        text = self._report_text(tmp_path, demo=False, source_versions=None)
        assert "unavailable (no variant was successfully annotated this run)" in text
        assert "v4.1" not in text
        assert "09/2025" not in text

    def test_live_mode_flags_missing_clinvar_entry_even_if_annotation_ran(self, tmp_path):
        # Some annotation succeeded (source_versions is not None) but Ensembl's
        # response happened to carry no ClinVar entry this run.
        text = self._report_text(tmp_path, demo=False, source_versions={"dbsnp": "156"})
        assert "ClinVar | unavailable (Ensembl did not report a ClinVar version this run)" in text
        assert "gnomAD | v4.1 (via Ensembl VEP colocated variants)" in text

    def test_demo_mode_labels_as_demo_cache(self, tmp_path, monkeypatch):
        import requests

        def fail_if_called(*args, **kwargs):
            raise AssertionError("demo mode must not touch the network")

        monkeypatch.setattr(requests, "get", fail_if_called)

        text = self._report_text(tmp_path, demo=True)
        assert "ClinVar | 2025-03-01 release (demo cache)" in text
        assert "gnomAD | v4.1 (demo cache)" in text

    def test_demo_pipeline_makes_zero_network_calls(self, tmp_path, monkeypatch):
        """End-to-end demo mode (parse -> classify -> report) must stay fully
        offline: patch both requests.get and requests.post to fail loudly."""
        import requests
        from clinical_variant_reporter import (
            generate_report, parse_vcf, run_classification,
        )

        def fail_if_called(*args, **kwargs):
            raise AssertionError("demo pipeline must not touch the network")

        monkeypatch.setattr(requests, "get", fail_if_called)
        monkeypatch.setattr(requests, "post", fail_if_called)

        demo_vcf = SKILL_DIR / "example_data" / "giab_acmg_panel.vcf"
        records = parse_vcf(demo_vcf)
        classified, source_versions = run_classification(records, demo=True)
        assert source_versions is None
        generate_report(classified, tmp_path, demo=True, source_versions=source_versions)
        assert (tmp_path / "report.md").exists()


# ---------------------------------------------------------------------------
# Regression — annotation-time provenance: version capture must be tied to
# whether VEP annotation actually succeeded, not to report-generation time
# ---------------------------------------------------------------------------
class TestAnnotateCapturesSourceVersions:
    """annotate_variants_vep must capture Ensembl's data-source versions once,
    only when at least one VEP batch actually succeeded, so a report can never
    claim provenance for annotation that never happened (issue #303 review)."""

    def _vcf_record(self):
        from clinical_variant_reporter import VcfRecord
        return VcfRecord(chrom="1", pos=100, id=".", ref="A", alt="G",
                          qual=".", filt="PASS", info={})

    def _mock_vep_success(self):
        from unittest.mock import MagicMock
        resp = MagicMock()
        resp.status_code = 200
        resp.json.return_value = [{
            "input": "1 100 100 A/G 1",
            "most_severe_consequence": "missense_variant",
            "transcript_consequences": [{
                "gene_symbol": "FAKEGENE", "impact": "MODERATE",
                "consequence_terms": ["missense_variant"],
                "transcript_id": "ENST00000000001.3",
            }],
        }]
        return resp

    def test_source_versions_captured_when_batch_succeeds(self, monkeypatch):
        from unittest.mock import MagicMock
        import requests
        from clinical_variant_reporter import annotate_variants_vep

        monkeypatch.setattr(requests, "post", lambda *a, **k: self._mock_vep_success())
        monkeypatch.setattr("clinical_variant_reporter.VEP_RATE_LIMIT_SECONDS", 0)

        info_resp = MagicMock()
        info_resp.json.return_value = [
            {"name": "ClinVar", "version": "09/2025"},
            {"name": "dbSNP", "version": "156"},
            {"name": "OMIM", "version": "09/2025"},
        ]
        monkeypatch.setattr(requests, "get", lambda *a, **k: info_resp)

        evidence_list, source_versions = annotate_variants_vep([self._vcf_record()])
        assert len(evidence_list) == 1
        assert source_versions == {"clinvar": "09/2025", "dbsnp": "156", "omim": "09/2025"}

    def test_source_versions_none_when_every_batch_fails(self, monkeypatch):
        import requests
        from clinical_variant_reporter import annotate_variants_vep

        def raise_timeout(*args, **kwargs):
            raise requests.exceptions.Timeout("no response")

        def fail_if_called(*args, **kwargs):
            raise AssertionError("version fetch must not run when nothing was annotated")

        monkeypatch.setattr(requests, "post", raise_timeout)
        monkeypatch.setattr(requests, "get", fail_if_called)
        monkeypatch.setattr("clinical_variant_reporter.VEP_RATE_LIMIT_SECONDS", 0)

        evidence_list, source_versions = annotate_variants_vep([self._vcf_record()])
        assert len(evidence_list) == 1  # unannotated placeholder, still returned
        assert source_versions is None


# ---------------------------------------------------------------------------
# Regression — _fetch_ensembl_data_versions error handling and parsing
# ---------------------------------------------------------------------------
class TestFetchEnsemblDataVersions:
    """The live-mode version fetch must name the exception and warn on stderr
    (matching the VEP handler's style), and must not crash on a 4xx/5xx
    response, a non-JSON body, or a payload missing the ClinVar entry."""

    def test_http_error_returns_none_and_warns(self, monkeypatch, capsys):
        import requests
        from clinical_variant_reporter import _fetch_ensembl_data_versions

        def raise_http_error(*args, **kwargs):
            raise requests.exceptions.HTTPError("500 Server Error")

        monkeypatch.setattr(requests, "get", raise_http_error)

        assert _fetch_ensembl_data_versions() is None
        assert "WARNING" in capsys.readouterr().err

    def test_non_json_body_returns_none(self, monkeypatch):
        from unittest.mock import MagicMock
        import requests
        from clinical_variant_reporter import _fetch_ensembl_data_versions

        resp = MagicMock()
        resp.json.side_effect = ValueError("not JSON")
        monkeypatch.setattr(requests, "get", lambda *a, **k: resp)

        assert _fetch_ensembl_data_versions() is None

    def test_payload_missing_clinvar_entry(self, monkeypatch):
        from unittest.mock import MagicMock
        import requests
        from clinical_variant_reporter import _fetch_ensembl_data_versions

        resp = MagicMock()
        resp.json.return_value = [{"name": "dbSNP", "version": "156"}]
        monkeypatch.setattr(requests, "get", lambda *a, **k: resp)

        versions = _fetch_ensembl_data_versions()
        assert versions == {"dbsnp": "156"}
        assert "clinvar" not in versions
