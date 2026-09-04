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
    interpret_clinvar_significance,
    is_secondary_finding_gene,
    normalize_clinvar_terms,
)
from clinical_variant_reporter import (
    ENSEMBL_INFO_VARIATION_PATH,
    VEP_REST_HOST,
    VEP_REST_HOST_GRCH37,
    VEP_REST_PATH,
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

    @staticmethod
    def _mock_vep_response():
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

    def test_vep_post_uses_grch38_host_by_default(self, monkeypatch):
        """rest.ensembl.org always serves GRCh38 regardless of the assembly
        query param, so the default (no --assembly) run must hit that host."""
        from clinical_variant_reporter import VcfRecord

        captured_urls: list[str] = []

        def mock_post(*args, **kwargs):
            captured_urls.append(args[0] if args else kwargs.get("url"))
            return self._mock_vep_response()

        import requests
        monkeypatch.setattr(requests, "post", mock_post)
        monkeypatch.setattr("clinical_variant_reporter.VEP_RATE_LIMIT_SECONDS", 0)

        records = [VcfRecord(chrom="1", pos=100, id=".", ref="A", alt="G",
                             qual=".", filt="PASS", info={})]
        annotate_variants_vep(records, assembly="GRCh38")

        assert captured_urls == [VEP_REST_HOST + VEP_REST_PATH]

    def test_vep_post_switches_host_for_grch37(self, monkeypatch):
        """--assembly GRCh37 must POST to grch37.rest.ensembl.org, not the
        main host, which silently serves GRCh38 regardless of the assembly
        query param (review on #327: 'derive the host from assembly for
        both the VEP POST and the version GET')."""
        from clinical_variant_reporter import VcfRecord

        captured_urls: list[str] = []

        def mock_post(*args, **kwargs):
            captured_urls.append(args[0] if args else kwargs.get("url"))
            return self._mock_vep_response()

        import requests
        monkeypatch.setattr(requests, "post", mock_post)
        monkeypatch.setattr("clinical_variant_reporter.VEP_RATE_LIMIT_SECONDS", 0)

        records = [VcfRecord(chrom="1", pos=100, id=".", ref="A", alt="G",
                             qual=".", filt="PASS", info={})]
        annotate_variants_vep(records, assembly="GRCh37")

        assert captured_urls == [VEP_REST_HOST_GRCH37 + VEP_REST_PATH]
        assert captured_urls != [VEP_REST_HOST + VEP_REST_PATH]


# ---------------------------------------------------------------------------
# Regression — VEP clin_sig_allele is an ALT-linked string
# ---------------------------------------------------------------------------
class TestClinSigAlleleTypeGuard:
    """VEP returns colocated_variants[].clin_sig_allele as an ALT-linked
    string such as ``T:pathogenic``. Unsupported shapes must fail closed rather
    than being mistaken for review metadata."""

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
        assert ev.clinvar_significance == "pathogenic"
        assert ev.clinvar_review_stars == 0  # REST does not link stars to this ALT

    def test_dict_clin_sig_allele_fails_closed_without_an_alt_assertion(self):
        from clinical_variant_reporter import _extract_evidence_from_vep
        ev = _extract_evidence_from_vep(
            self._vep({"review_status_stars": 3}), self._record())
        assert ev.clinvar_significance == ""
        assert ev.clinvar_review_stars == 0


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
        assert "ClinVar | 09/2025 (via Ensembl /info/variation/homo_sapiens)" in text
        assert "gnomAD | v4.1 (hardcoded; Ensembl serves no gnomAD version endpoint)" in text
        assert "2025-03-01 release" not in text

    def test_live_mode_suppresses_claim_when_nothing_annotated(self, tmp_path):
        # source_versions is None: annotate_variants_vep never had a
        # successful VEP batch, so no source, including gnomAD's own
        # hardcoded constant, may be claimed for this run.
        text = self._report_text(tmp_path, demo=False, source_versions=None)
        assert "unavailable (no variant was successfully annotated this run)" in text
        assert "v4.1" not in text
        assert "09/2025" not in text

    def test_live_mode_flags_missing_clinvar_entry_even_if_annotation_ran(self, tmp_path):
        # Some annotation succeeded (source_versions is not None) but Ensembl's
        # response happened to carry no ClinVar entry this run.
        text = self._report_text(tmp_path, demo=False, source_versions={"dbsnp": "156"})
        assert "ClinVar | unavailable (Ensembl did not report a ClinVar version this run)" in text
        assert "gnomAD | v4.1 (hardcoded; Ensembl serves no gnomAD version endpoint)" in text

    def test_live_mode_flags_version_lookup_failure_distinctly_from_no_annotation(self, tmp_path):
        # source_versions == {} (empty dict, not None): at least one VEP batch
        # succeeded but the separate Ensembl version lookup itself failed.
        # Before this fix both states printed the identical "no variant was
        # successfully annotated this run" sentence, which is false when
        # annotation actually ran (review on #327, blocking item 1).
        text = self._report_text(tmp_path, demo=False, source_versions={})
        assert "ClinVar | unavailable (annotation succeeded, but the Ensembl version lookup failed)" in text
        assert "no variant was successfully annotated this run" not in text
        # gnomAD is a hardcoded constant independent of this lookup, so
        # annotation having succeeded still lets it be claimed.
        assert "gnomAD | v4.1 (hardcoded; Ensembl serves no gnomAD version endpoint)" in text

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
# Regression — result.json / database_versions.json must carry the
# structured ClinVar/dbSNP/OMIM/gnomAD versions, not only prose (review on
# #327 smaller items: "Keep the audit JSON machine-readable")
# ---------------------------------------------------------------------------
class TestStructuredDataSourceVersions:
    def test_live_mode_writes_structured_versions_to_result_and_db_json(self, tmp_path):
        generate_report(
            [], tmp_path, demo=False,
            source_versions={"clinvar": "09/2025", "dbsnp": "156", "omim": "09/2025"},
        )
        result = json.loads((tmp_path / "result.json").read_text())
        assert result["data_source_versions"] == {
            "clinvar": "09/2025", "dbsnp": "156", "omim": "09/2025", "gnomad": "v4.1",
        }
        db_versions = json.loads(
            (tmp_path / "reproducibility" / "database_versions.json").read_text()
        )
        assert db_versions["data_source_versions"] == {
            "clinvar": "09/2025", "dbsnp": "156", "omim": "09/2025", "gnomad": "v4.1",
        }

    def test_live_mode_nulls_missing_sources_rather_than_prose(self, tmp_path):
        # source_versions == {}: annotation succeeded, lookup failed. The
        # structured form must use null, not the "unavailable (...)" prose
        # that the Markdown table shows, so it stays machine-parseable.
        generate_report([], tmp_path, demo=False, source_versions={})
        result = json.loads((tmp_path / "result.json").read_text())
        assert result["data_source_versions"] == {
            "clinvar": None, "dbsnp": None, "omim": None, "gnomad": "v4.1",
        }

    def test_live_mode_nulls_gnomad_too_when_nothing_annotated(self, tmp_path):
        # source_versions is None: no VEP batch succeeded this run, so no
        # source-version claim can be made at all -- including gnomAD's
        # hardcoded constant. report.md already says "unavailable (no
        # variant was successfully annotated this run)" for both ClinVar and
        # gnomAD in this case (_data_source_versions_for_report's `note`
        # branch); result.json previously still asserted the gnomAD constant
        # here, so the two artefacts of one run disagreed on whether a
        # gnomAD claim may be made (review on #327, blocking item 2).
        generate_report([], tmp_path, demo=False, source_versions=None)
        result = json.loads((tmp_path / "result.json").read_text())
        assert result["data_source_versions"] == {
            "clinvar": None, "dbsnp": None, "omim": None, "gnomad": None,
        }

    def test_annotation_backend_stamps_the_assembly_actually_used(self, tmp_path):
        # review #327: annotation_backend hardcoded "GRCh38" regardless of
        # the assembly argument, so a GRCh37 run's own database_versions.json
        # contradicted its GRCh37 data_source_versions (a ClinVar version
        # that host never served, stamped under a GRCh38 backend label).
        generate_report([], tmp_path, demo=False, assembly="GRCh37")
        db_versions = json.loads(
            (tmp_path / "reproducibility" / "database_versions.json").read_text()
        )
        assert db_versions["annotation_backend"] == "Ensembl VEP REST (GRCh37)"

    def test_demo_mode_structured_versions(self, tmp_path):
        generate_report([], tmp_path, demo=True)
        result = json.loads((tmp_path / "result.json").read_text())
        assert result["data_source_versions"]["clinvar"] == "2025-03-01 release"
        assert result["data_source_versions"]["gnomad"] == "v4.1"


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

        captured_args: list[tuple] = []
        captured_kwargs: list[dict] = []

        info_resp = MagicMock()
        info_resp.json.return_value = [
            {"name": "ClinVar", "version": "09/2025"},
            {"name": "dbSNP", "version": "156"},
            {"name": "OMIM", "version": "09/2025"},
        ]

        def mock_get(*args, **kwargs):
            captured_args.append(args)
            captured_kwargs.append(kwargs)
            return info_resp

        monkeypatch.setattr(requests, "get", mock_get)

        evidence_list, source_versions = annotate_variants_vep([self._vcf_record()])
        assert len(evidence_list) == 1
        assert source_versions == {"clinvar": "09/2025", "dbsnp": "156", "omim": "09/2025"}

        # The suite must not pass identically if the code reverts to the old
        # /info/rest endpoint (review on #327, blocking item 3): assert the
        # actual URL called, not just that requests.get was called at all.
        # Pinned against the LITERAL path, not just ENSEMBL_INFO_VARIATION_PATH --
        # a constant-relative-only assertion stays green even if the constant
        # itself regresses to the wrong endpoint (mutation-tested on review:
        # changing the constant to "/info/software" left 72 tests unchanged).
        assert len(captured_args) == 1
        called_url = captured_args[0][0] if captured_args[0] else captured_kwargs[0].get("url")
        assert called_url == "https://rest.ensembl.org/info/variation/homo_sapiens"
        assert called_url == VEP_REST_HOST + ENSEMBL_INFO_VARIATION_PATH
        assert captured_kwargs[0].get("timeout") == 10

    def test_batch_failure_raises_before_version_fetch(self, monkeypatch):
        import requests
        from clinical_variant_reporter import annotate_variants_vep

        def raise_timeout(*args, **kwargs):
            raise requests.exceptions.Timeout("no response")

        def fail_if_called(*args, **kwargs):
            raise AssertionError("version fetch must not run when nothing was annotated")

        monkeypatch.setattr(requests, "post", raise_timeout)
        monkeypatch.setattr(requests, "get", fail_if_called)
        monkeypatch.setattr("clinical_variant_reporter.VEP_RATE_LIMIT_SECONDS", 0)

        with pytest.raises(RuntimeError, match="VEP annotation failed for a batch of 1 variants"):
            annotate_variants_vep([self._vcf_record()])

    def test_complete_annotation_success_still_returns_evidence_and_versions(self, monkeypatch):
        from unittest.mock import MagicMock
        import requests
        from clinical_variant_reporter import VcfRecord, annotate_variants_vep

        records = [
            self._vcf_record(),
            VcfRecord(chrom="2", pos=200, id=".", ref="C", alt="T", qual=".", filt="PASS", info={}),
        ]

        resp = MagicMock()
        resp.raise_for_status.return_value = None
        resp.json.return_value = [
            {
                "input": "1 100 100 A/G 1",
                "most_severe_consequence": "missense_variant",
                "transcript_consequences": [{"gene_symbol": "GENE1", "consequence_terms": ["missense_variant"]}],
            },
            {
                "input": "2 200 200 C/T 1",
                "most_severe_consequence": "synonymous_variant",
                "transcript_consequences": [{"gene_symbol": "GENE2", "consequence_terms": ["synonymous_variant"]}],
            },
        ]
        info_resp = MagicMock()
        info_resp.json.return_value = [{"name": "ClinVar", "version": "09/2025"}]

        monkeypatch.setattr(requests, "post", lambda *a, **k: resp)
        monkeypatch.setattr(requests, "get", lambda *a, **k: info_resp)
        monkeypatch.setattr("clinical_variant_reporter.VEP_RATE_LIMIT_SECONDS", 0)

        evidence_list, source_versions = annotate_variants_vep(records)

        assert [ev.gene for ev in evidence_list] == ["GENE1", "GENE2"]
        assert source_versions == {"clinvar": "09/2025"}

    def test_partial_batch_failure_raises_instead_of_returning_partial_results(self, monkeypatch):
        import requests
        from clinical_variant_reporter import VcfRecord, annotate_variants_vep

        records = [
            self._vcf_record(),
            VcfRecord(chrom="2", pos=200, id=".", ref="C", alt="T", qual=".", filt="PASS", info={}),
        ]
        calls = {"count": 0}

        def post_once_then_fail(*args, **kwargs):
            calls["count"] += 1
            if calls["count"] == 1:
                return self._mock_vep_success()
            raise requests.exceptions.ConnectionError("transient VEP outage")

        def fail_if_called(*args, **kwargs):
            raise AssertionError("version fetch must not run after a failed annotation batch")

        monkeypatch.setattr("clinical_variant_reporter.VEP_BATCH_SIZE", 1)
        monkeypatch.setattr("clinical_variant_reporter.VEP_RATE_LIMIT_SECONDS", 0)
        monkeypatch.setattr(requests, "post", post_once_then_fail)
        monkeypatch.setattr(requests, "get", fail_if_called)

        with pytest.raises(RuntimeError, match="batch starting at index 1"):
            annotate_variants_vep(records)

    def test_malformed_vep_json_raises_instead_of_classifying_unannotated(self, monkeypatch):
        from unittest.mock import MagicMock
        import requests
        from clinical_variant_reporter import annotate_variants_vep

        resp = MagicMock()
        resp.raise_for_status.return_value = None
        resp.json.return_value = {"error": "not the expected list payload"}

        monkeypatch.setattr(requests, "post", lambda *a, **k: resp)
        monkeypatch.setattr(requests, "get", lambda *a, **k: None)
        monkeypatch.setattr("clinical_variant_reporter.VEP_RATE_LIMIT_SECONDS", 0)

        with pytest.raises(RuntimeError, match="Malformed VEP response"):
            annotate_variants_vep([self._vcf_record()])

    def test_cli_propagates_annotation_failure_without_writing_report(self, monkeypatch, tmp_path):
        import requests
        import clinical_variant_reporter

        vcf_path = tmp_path / "input.vcf"
        vcf_path.write_text(
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "1\t100\t.\tA\tG\t.\tPASS\t.\n",
            encoding="utf-8",
        )
        output_dir = tmp_path / "out"

        def raise_timeout(*args, **kwargs):
            raise requests.exceptions.Timeout("no response")

        monkeypatch.setattr(requests, "post", raise_timeout)
        monkeypatch.setattr(requests, "get", lambda *a, **k: None)
        monkeypatch.setattr("clinical_variant_reporter.VEP_RATE_LIMIT_SECONDS", 0)
        monkeypatch.setattr(
            sys,
            "argv",
            ["clinical_variant_reporter.py", "--input", str(vcf_path), "--output", str(output_dir)],
        )

        with pytest.raises(RuntimeError, match="VEP annotation failed"):
            clinical_variant_reporter.main()

        assert not (output_dir / "result.json").exists()

    def test_source_versions_empty_dict_when_annotation_succeeds_but_lookup_fails(self, monkeypatch):
        """The bug behind blocking item 1: before this fix, a failed version
        lookup after successful annotation was indistinguishable from no
        annotation at all (both produced source_versions is None). Now the
        two states must differ: None means nothing was annotated, {} means
        annotation succeeded but the lookup itself failed."""
        import requests
        from clinical_variant_reporter import annotate_variants_vep

        monkeypatch.setattr(requests, "post", lambda *a, **k: self._mock_vep_success())
        monkeypatch.setattr("clinical_variant_reporter.VEP_RATE_LIMIT_SECONDS", 0)

        def raise_connection_error(*args, **kwargs):
            raise requests.exceptions.ConnectionError("version host unreachable")

        monkeypatch.setattr(requests, "get", raise_connection_error)

        evidence_list, source_versions = annotate_variants_vep([self._vcf_record()])
        assert len(evidence_list) == 1
        assert source_versions == {}
        assert source_versions is not None

    def test_version_get_uses_grch38_host_by_default(self, monkeypatch):
        from unittest.mock import MagicMock
        import requests
        from clinical_variant_reporter import annotate_variants_vep

        monkeypatch.setattr(requests, "post", lambda *a, **k: self._mock_vep_success())
        monkeypatch.setattr("clinical_variant_reporter.VEP_RATE_LIMIT_SECONDS", 0)

        captured_urls: list[str] = []
        info_resp = MagicMock()
        info_resp.json.return_value = [{"name": "ClinVar", "version": "09/2025"}]

        def mock_get(*args, **kwargs):
            captured_urls.append(args[0] if args else kwargs.get("url"))
            return info_resp

        monkeypatch.setattr(requests, "get", mock_get)

        annotate_variants_vep([self._vcf_record()], assembly="GRCh38")
        assert captured_urls == [VEP_REST_HOST + ENSEMBL_INFO_VARIATION_PATH]

    def test_version_get_switches_host_for_grch37(self, monkeypatch):
        """grch37.rest.ensembl.org and rest.ensembl.org disagree on ClinVar's
        version (measured in the review: 06/2023 vs 09/2025, 28 months
        apart), so a --assembly GRCh37 run must fetch from the GRCh37 host
        or it certifies a version that was never served for GRCh37 data
        (review on #327, blocking item 2)."""
        from unittest.mock import MagicMock
        import requests
        from clinical_variant_reporter import annotate_variants_vep

        monkeypatch.setattr(requests, "post", lambda *a, **k: self._mock_vep_success())
        monkeypatch.setattr("clinical_variant_reporter.VEP_RATE_LIMIT_SECONDS", 0)

        captured_urls: list[str] = []
        info_resp = MagicMock()
        info_resp.json.return_value = [{"name": "ClinVar", "version": "06/2023"}]

        def mock_get(*args, **kwargs):
            captured_urls.append(args[0] if args else kwargs.get("url"))
            return info_resp

        monkeypatch.setattr(requests, "get", mock_get)

        evidence_list, source_versions = annotate_variants_vep(
            [self._vcf_record()], assembly="GRCh37",
        )
        assert captured_urls == [VEP_REST_HOST_GRCH37 + ENSEMBL_INFO_VARIATION_PATH]
        assert captured_urls != [VEP_REST_HOST + ENSEMBL_INFO_VARIATION_PATH]
        assert source_versions == {"clinvar": "06/2023"}


# ---------------------------------------------------------------------------
# Regression — _fetch_ensembl_data_versions error handling and parsing
# ---------------------------------------------------------------------------
class TestFetchEnsemblDataVersions:
    """The live-mode version fetch must name the exception and warn on stderr
    (matching the VEP handler's style), and must not crash on a 4xx/5xx
    response, a non-JSON body, or a payload missing the ClinVar entry."""

    def test_http_error_returns_none_and_warns(self, monkeypatch, capsys):
        """requests.get itself raising (network error, DNS failure, etc.),
        as opposed to a real HTTP response whose raise_for_status() raises
        (covered separately below)."""
        import requests
        from clinical_variant_reporter import _fetch_ensembl_data_versions

        captured_kwargs: list[dict] = []

        def raise_http_error(*args, **kwargs):
            captured_kwargs.append(kwargs)
            raise requests.exceptions.HTTPError("500 Server Error")

        monkeypatch.setattr(requests, "get", raise_http_error)

        assert _fetch_ensembl_data_versions() is None
        assert "WARNING" in capsys.readouterr().err
        assert captured_kwargs[0].get("timeout") == 10

    def test_raise_for_status_real_failure_returns_none_and_warns(self, monkeypatch, capsys):
        """The gap the review named explicitly: a MagicMock's
        raise_for_status() is a silent no-op by default, so a prior version
        of this suite never actually exercised the resp.raise_for_status()
        call in _fetch_ensembl_data_versions. This test uses a response
        object whose raise_for_status() genuinely raises for a 4xx/5xx
        status, driving the real code path rather than short-circuiting at
        requests.get() itself."""
        import requests
        from clinical_variant_reporter import _fetch_ensembl_data_versions

        class RealFailureResponse:
            status_code = 429

            def raise_for_status(self):
                raise requests.exceptions.HTTPError("429 Too Many Requests")

            def json(self):
                raise AssertionError("json() must not be called after raise_for_status() raises")

        monkeypatch.setattr(requests, "get", lambda *a, **k: RealFailureResponse())

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

        captured_args: list[tuple] = []
        captured_kwargs: list[dict] = []

        resp = MagicMock()
        resp.json.return_value = [{"name": "dbSNP", "version": "156"}]

        def mock_get(*args, **kwargs):
            captured_args.append(args)
            captured_kwargs.append(kwargs)
            return resp

        monkeypatch.setattr(requests, "get", mock_get)

        versions = _fetch_ensembl_data_versions()
        assert versions == {"dbsnp": "156"}
        assert "clinvar" not in versions

        called_url = captured_args[0][0] if captured_args[0] else captured_kwargs[0].get("url")
        assert called_url == VEP_REST_HOST + ENSEMBL_INFO_VARIATION_PATH
        assert captured_kwargs[0].get("timeout") == 10
        assert captured_kwargs[0].get("headers", {}).get("Accept") == "application/json"


# ---------------------------------------------------------------------------
# Regression — VEP clin_sig arrives as a list; ClinVar terms are parsed, not
# joined and substring-matched (issue #328)
# ---------------------------------------------------------------------------
def _legacy_directions(significance: str) -> tuple[bool, bool]:
    """The substring rule this fix replaces (acmg_engine.py before #328).

    Kept as the regression oracle for scalar inputs: on every string the old
    rule read correctly, the structured parser must agree with it. The two
    deliberate divergences (mixed pathogenic + benign strings; benign plus
    conflicting_data_from_submitters, which the old BP6 never checked for)
    have their own test.
    """
    lowered = significance.lower()
    pathogenic = "pathogenic" in lowered and "conflicting" not in lowered
    benign = ("benign" in lowered or "likely benign" in lowered) and "pathogenic" not in lowered
    return pathogenic, benign


def _clinvar_evidence(significance, stars: int = 3, **overrides) -> VariantEvidence:
    fields = dict(
        chrom="17", pos=43093464, ref="A", alt="G", gene="BRCA1",
        consequence="missense_variant", is_missense=True,
        clinvar_significance=significance, clinvar_review_stars=stars,
    )
    fields.update(overrides)
    return VariantEvidence(**fields)


def _criterion(ev: VariantEvidence, code: str) -> EvidenceCriterion:
    return next(c for c in evaluate_criteria(ev) if c.code == code)


class TestClinVarTermNormalisation:
    """Both shapes ClinVar data reaches the engine in must reduce to the same
    ordered, de-duplicated tuple of snake_case terms: VEP's JSON list, and the
    scalar strings ClinVar itself emits, which join several terms with "/"
    (compound), "|" (VCF CLNSIG multi-value), ";" (variant_summary) or ", "."""

    @pytest.mark.parametrize("raw, expected", [
        ("Pathogenic", ("pathogenic",)),
        ("Likely pathogenic", ("likely_pathogenic",)),
        (["pathogenic", "pathogenic"], ("pathogenic",)),
        (["Pathogenic", "pathogenic"], ("pathogenic",)),
        ("Pathogenic/Likely pathogenic", ("pathogenic", "likely_pathogenic")),
        ("Pathogenic/Likely_pathogenic", ("pathogenic", "likely_pathogenic")),
        ("Pathogenic|risk_factor", ("pathogenic", "risk_factor")),
        ("Benign/Likely benign", ("benign", "likely_benign")),
        ("Uncertain significance; risk factor", ("uncertain_significance", "risk_factor")),
        ("Conflicting interpretations of pathogenicity, Pathogenic",
         ("conflicting_interpretations_of_pathogenicity", "pathogenic")),
        # clin_sig_allele spelling of "Pathogenic, low penetrance" (VEP 116, rs1800562)
        ("pathogenic/pathogenic,_low_penetrance", ("pathogenic", "low_penetrance")),
        (["Pathogenic/Likely_pathogenic", "risk_factor"],
         ("pathogenic", "likely_pathogenic", "risk_factor")),
        (["uncertain_significance", "pathogenic"], ("uncertain_significance", "pathogenic")),
        (("benign",), ("benign",)),
        ("", ()),
        (None, ()),
        ([], ()),
        (["", None, "  "], ()),
    ])
    def test_terms_are_split_normalised_and_deduplicated(self, raw, expected):
        assert normalize_clinvar_terms(raw) == expected

    def test_unsupported_shape_fails_loudly(self):
        with pytest.raises(TypeError):
            normalize_clinvar_terms({"T": "pathogenic"})


class TestClinVarAssertionSemantics:
    """Explicit semantics for the three ClinVar-backed rules (PS1, PP5, BP6):

    - a term from the pathogenic family supports pathogenic, one from the
      benign family supports benign;
    - any term ClinVar or VEP uses for "submitters disagree" (every spelling
      in the wild) is a conflict;
    - pathogenic-family and benign-family terms co-occurring is a conflict
      even when no such token is present (VEP's clin_sig list aggregates
      every ClinVar record at the site);
    - a conflict supports neither direction;
    - anything else (uncertain, drug response, risk factor, unknown) is
      neutral: it neither supports nor blocks.
    """

    # (input, supports_pathogenic, supports_benign, is_conflicting)
    CASES = [
        # single terms
        ("Pathogenic", True, False, False),
        ("Likely pathogenic", True, False, False),
        ("Benign", False, True, False),
        ("Likely benign", False, True, False),
        ("Uncertain significance", False, False, False),
        ("drug_response", False, False, False),
        ("not_provided", False, False, False),
        ("novel_value_from_a_future_release", False, False, False),
        ("", False, False, False),
        (None, False, False, False),
        # duplicates
        (["pathogenic", "pathogenic"], True, False, False),
        (["benign", "Benign"], False, True, False),
        # ClinVar compound / multi-value scalar forms
        ("Pathogenic/Likely pathogenic", True, False, False),
        ("Pathogenic/Likely_pathogenic", True, False, False),
        ("Pathogenic|risk_factor", True, False, False),
        ("Benign/Likely benign", False, True, False),
        ("Pathogenic, low penetrance", True, False, False),
        # VEP list vocabulary as served by rest.ensembl.org release 116
        (["pathogenic_low_penetrance"], True, False, False),
        (["likely_pathogenic_low_penetrance"], True, False, False),
        ("Likely pathogenic, low penetrance", True, False, False),
        (["uncertain_significance", "pathogenic"], True, False, False),
        (["uncertain_significance", "not_provided", "pathogenic", "other",
          "risk_factor", "pathogenic_low_penetrance"], True, False, False),
        (["uncertain_significance", "benign", "likely_benign"], False, True, False),
        # explicit conflicting tokens, every spelling in the wild
        (["pathogenic", "conflicting_interpretations_of_pathogenicity"], False, False, True),
        (["pathogenic", "conflicting_classifications_of_pathogenicity"], False, False, True),
        (["benign", "conflicting_data_from_submitters"], False, False, True),
        ("Conflicting_classifications_of_pathogenicity", False, False, True),
        ("Conflicting interpretations of pathogenicity, Pathogenic", False, False, True),
        ("Conflicting_interpretations_of_pathogenicity|Pathogenic", False, False, True),
        # both families present is a conflict whatever the tokens say
        (["benign", "pathogenic"], False, False, True),
        (["likely_benign", "likely_pathogenic"], False, False, True),
        ("Benign; Pathogenic", False, False, True),
        ("Pathogenic/Likely_pathogenic|Benign", False, False, True),
        # rs1042522 aggregate at release 116: C is pathogenic, T is benign
        (["uncertain_significance", "benign", "likely_benign", "pathogenic"], False, False, True),
    ]

    @pytest.mark.parametrize("raw, pathogenic, benign, conflicting", CASES)
    def test_directions(self, raw, pathogenic, benign, conflicting):
        assertion = interpret_clinvar_significance(raw)
        observed = (assertion.supports_pathogenic, assertion.supports_benign, assertion.is_conflicting)
        assert observed == (pathogenic, benign, conflicting)

    def test_a_conflict_never_supports_either_direction(self):
        for raw, *_ in self.CASES:
            assertion = interpret_clinvar_significance(raw)
            if assertion.is_conflicting:
                assert not assertion.supports_pathogenic
                assert not assertion.supports_benign

    def test_explicit_token_is_load_bearing_without_a_benign_term(self):
        """No benign term here, so only the explicit-token clause can block.
        Deleting that clause turns this test red."""
        assertion = interpret_clinvar_significance(
            ["pathogenic", "conflicting_classifications_of_pathogenicity"])
        assert assertion.is_conflicting
        assert assertion.supports_pathogenic is False
        assert assertion.conflicting_terms == ("conflicting_classifications_of_pathogenicity",)

    def test_cooccurrence_is_load_bearing_without_an_explicit_token(self):
        """No conflicting token here, so only the co-occurrence clause can block.
        Deleting that clause turns this test red."""
        assertion = interpret_clinvar_significance(["benign", "pathogenic"])
        assert assertion.is_conflicting
        assert assertion.conflicting_terms == ()
        assert assertion.pathogenic_terms == ("pathogenic",)
        assert assertion.benign_terms == ("benign",)

    def test_unseen_conflicting_spelling_is_still_a_conflict(self):
        """ClinVar renamed interpretations->classifications in 2024; match the
        family stem so the next rename does not reopen this bug."""
        assertion = interpret_clinvar_significance(["pathogenic", "conflicting_something_new"])
        assert assertion.is_conflicting

    def test_conflict_reason_is_human_readable(self):
        explicit = interpret_clinvar_significance(
            ["pathogenic", "conflicting_classifications_of_pathogenicity"])
        assert "conflicting_classifications_of_pathogenicity" in explicit.conflict_reason
        mixed = interpret_clinvar_significance(["benign", "pathogenic"])
        assert "pathogenic" in mixed.conflict_reason and "benign" in mixed.conflict_reason
        assert interpret_clinvar_significance("Pathogenic").conflict_reason == ""

    SCALAR_PARITY = [
        "Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic",
        "Pathogenic/Likely_pathogenic", "Pathogenic|risk_factor",
        "Pathogenic, low penetrance", "Likely pathogenic, low penetrance",
        "Benign", "Likely benign",
        "Benign/Likely benign", "Uncertain significance", "risk factor",
        "drug response", "not provided", "",
        "Conflicting interpretations of pathogenicity, Pathogenic",
        "Conflicting_classifications_of_pathogenicity",
        "Conflicting_interpretations_of_pathogenicity|Pathogenic",
        "Benign|Conflicting_interpretations_of_pathogenicity",
        "Benign|Conflicting_classifications_of_pathogenicity",
        "conflicting data from submitters",
    ]

    @pytest.mark.parametrize("raw", SCALAR_PARITY)
    def test_parity_with_the_replaced_rule_on_scalar_inputs(self, raw):
        assertion = interpret_clinvar_significance(raw)
        assert (assertion.supports_pathogenic, assertion.supports_benign) == _legacy_directions(raw)

    @pytest.mark.parametrize("raw, legacy", [
        # The old rule read "pathogenic" anywhere in a mixed string as
        # pathogenic support.
        ("Benign; Pathogenic", (True, False)),
        ("Pathogenic/Benign", (True, False)),
        ("Likely benign|Likely pathogenic", (True, False)),
        # The old BP6 had no conflicting check at all; the other two
        # conflicting tokens only blocked it because "pathogenicity" contains
        # "pathogenic". This token does not, so BP6 used to fire on it.
        ("Benign|conflicting_data_from_submitters", (False, True)),
        ("Likely_benign;conflicting_data_from_submitters", (False, True)),
    ])
    def test_deliberate_divergences_from_the_replaced_rule(self, raw, legacy):
        """The two intended behaviour changes on scalar input, pinned so the
        parity contract above is honest about what it excludes."""
        assert _legacy_directions(raw) == legacy
        assertion = interpret_clinvar_significance(raw)
        assert assertion.is_conflicting
        assert (assertion.supports_pathogenic, assertion.supports_benign) == (False, False)


class TestClinVarCriteriaWithStructuredTerms:
    def test_list_valued_pathogenic_fires_ps1_and_pp5(self):
        ev = _clinvar_evidence(["uncertain_significance", "pathogenic"], stars=3)
        assert _criterion(ev, "PS1").triggered is True
        assert _criterion(ev, "PP5").triggered is True
        assert _criterion(ev, "BP6").triggered is False

    def test_compound_scalars_keep_firing(self):
        assert _criterion(_clinvar_evidence("Pathogenic/Likely pathogenic", stars=3), "PP5").triggered
        assert _criterion(_clinvar_evidence("Pathogenic|risk_factor", stars=3), "PS1").triggered
        assert _criterion(_clinvar_evidence("Benign/Likely benign", stars=2), "BP6").triggered

    @pytest.mark.parametrize("raw", [
        ["pathogenic", "conflicting_interpretations_of_pathogenicity"],
        ["pathogenic", "conflicting_classifications_of_pathogenicity"],
        ["benign", "conflicting_classifications_of_pathogenicity"],
        ["benign", "pathogenic"],
        "Conflicting interpretations of pathogenicity, Pathogenic",
    ])
    def test_conflict_withholds_ps1_pp5_bp6_even_with_stars(self, raw):
        ev = _clinvar_evidence(raw, stars=4)
        for code in ("PS1", "PP5", "BP6"):
            crit = _criterion(ev, code)
            assert crit.triggered is False, code
            assert "conflict" in crit.detail.lower(), code

    def test_star_gate_is_unchanged(self):
        ev = _clinvar_evidence(["pathogenic"], stars=1)
        assert _criterion(ev, "PS1").triggered is False
        assert _criterion(ev, "PP5").triggered is False

    def test_criterion_source_renders_terms_not_python_repr(self):
        ev = _clinvar_evidence(["uncertain_significance", "pathogenic"])
        source = _criterion(ev, "PS1").source
        assert "['" not in source
        assert "ClinVar=uncertain_significance; pathogenic, stars=3" == source

    def test_scalar_source_is_byte_identical_to_before(self):
        ev = _clinvar_evidence("Pathogenic")
        assert _criterion(ev, "PP5").source == "ClinVar=Pathogenic, stars=3"

    def test_conflicting_list_lands_on_vus_not_likely_pathogenic(self):
        in_silico = dict(gnomad_af=None, cadd_phred=28.0, sift_prediction="deleterious")
        conflicted = _clinvar_evidence(
            ["pathogenic", "conflicting_classifications_of_pathogenicity"], stars=3, **in_silico)
        assert classify_variant(conflicted).classification == "Uncertain Significance"
        # Same evidence without the conflict token: PS1 + PM2 + PP3 + PP5 -> LP.
        clean = _clinvar_evidence(["uncertain_significance", "pathogenic"], stars=3, **in_silico)
        assert classify_variant(clean).classification == "Likely Pathogenic"


class TestListValuedClinSigEndToEnd:
    """Issue #328: VEP REST returns colocated_variants[].clin_sig as a JSON
    list. Before the fix _extract_evidence_from_vep stored it unchanged and
    every ClinVar rule called .lower() on it, so live mode raised
    AttributeError on the first ClinVar-annotated variant and wrote nothing."""

    def _record(self):
        from clinical_variant_reporter import VcfRecord
        return VcfRecord(chrom="13", pos=32339267, id="rs886040553",
                         ref="A", alt="T", qual=".", filt="PASS", info={})

    def _vep(self, clin_sig, clin_sig_allele, consequence="missense_variant", impact="MODERATE"):
        return {
            "most_severe_consequence": consequence,
            "transcript_consequences": [
                {"gene_symbol": "BRCA2", "impact": impact,
                 "consequence_terms": [consequence],
                 "transcript_id": "ENST00000544455.6",
                 "cadd_phred": 28.0, "sift_prediction": "deleterious"}
            ],
            "colocated_variants": [
                {"id": "rs886040553", "clin_sig": clin_sig, "clin_sig_allele": clin_sig_allele}
            ],
        }

    def test_release_116_payload_classifies_without_crashing(self):
        """Shape captured live from rest.ensembl.org (release 116) for
        13:32339267 A>T: a list clin_sig and a string clin_sig_allele."""
        from clinical_variant_reporter import _extract_evidence_from_vep
        payload = self._vep(
            ["uncertain_significance", "pathogenic"],
            "G:conflicting_classifications_of_pathogenicity;T:pathogenic;G:uncertain_significance",
        )
        ev = _extract_evidence_from_vep(payload, self._record())
        assert ev.clinvar_significance == "pathogenic"
        assert ev.clinvar_review_stars == 0  # REST payload carries no trustworthy stars
        cv = classify_variant(ev)  # raised AttributeError before the fix
        assert cv.classification == "Uncertain Significance"
        assert not {"PS1", "PP5", "BP6"} & set(cv.triggered_codes)

    def test_unlinked_star_dict_does_not_reenable_site_level_evidence(self):
        from clinical_variant_reporter import _extract_evidence_from_vep
        payload = self._vep(["benign", "pathogenic"], {"review_status_stars": 3})
        cv = classify_variant(_extract_evidence_from_vep(payload, self._record()))
        assert cv.evidence.clinvar_significance == ""
        assert cv.evidence.clinvar_review_stars == 0
        assert not {"PS1", "PP5", "BP6"} & set(cv.triggered_codes)
        assert cv.classification == "Uncertain Significance"

    def test_allele_specific_terms_without_review_stars_do_not_fire_clinvar_rules(self):
        from clinical_variant_reporter import _extract_evidence_from_vep
        payload = self._vep(
            ["uncertain_significance", "pathogenic"],
            "G:benign;T:pathogenic",
        )
        cv = classify_variant(_extract_evidence_from_vep(payload, self._record()))
        assert cv.evidence.clinvar_significance == "pathogenic"
        assert cv.evidence.clinvar_review_stars == 0
        assert not {"PS1", "PP5", "BP6"} & set(cv.triggered_codes)
        assert cv.classification == "Uncertain Significance"

    def test_evidence_exposes_display_text_for_both_shapes(self):
        assert _clinvar_evidence(["uncertain_significance", "pathogenic"]).clinvar_significance_text \
            == "uncertain_significance; pathogenic"
        assert _clinvar_evidence("Uncertain significance").clinvar_significance_text \
            == "Uncertain significance"
        assert _clinvar_evidence("").clinvar_significance_text == ""
        assert _clinvar_evidence(None).clinvar_significance_text == ""

    def test_report_and_table_render_terms_not_list_repr(self, tmp_path):
        cv = classify_variant(_clinvar_evidence(
            ["uncertain_significance", "pathogenic"], stars=3,
            gnomad_af=None, cadd_phred=28.0, sift_prediction="deleterious",
        ))
        assert cv.classification == "Likely Pathogenic"  # so it lands in the actionable section

        generate_report([cv], tmp_path, demo=False, source_versions={})
        report = (tmp_path / "report.md").read_text()
        table = (tmp_path / "tables" / "acmg_classifications.tsv").read_text()

        assert "['" not in report and "['" not in table
        assert "**ClinVar**: uncertain_significance; pathogenic (stars: 3)" in report
        row = next(line for line in table.splitlines() if line.startswith("17\t43093464"))
        assert "\tuncertain_significance; pathogenic\t3\t" in row


class TestAlleleSpecificClinVarSignificance:
    """Live VEP ``clin_sig`` is site-level; only the queried ALT is evidence.

    The API's ``clin_sig_allele`` payload is a semicolon-delimited set of
    ``ALT:terms`` entries.  A multi-allelic locus must never transfer a
    ClinVar assertion from a different alternate allele.
    """

    def _record(self, alt="T"):
        from clinical_variant_reporter import VcfRecord
        return VcfRecord(chrom="13", pos=32339267, id="rs886040553",
                         ref="A", alt=alt, qual=".", filt="PASS", info={})

    def _vep(self, clin_sig_allele):
        return {
            "most_severe_consequence": "missense_variant",
            "transcript_consequences": [{
                "gene_symbol": "BRCA2", "impact": "MODERATE",
                "consequence_terms": ["missense_variant"],
            }],
            "colocated_variants": [{
                # This aggregate deliberately contains a benign assertion for G
                # and a pathogenic assertion for T.  It is not allele-specific.
                "clin_sig": ["benign", "pathogenic"],
                "clin_sig_allele": clin_sig_allele,
            }],
        }

    def test_uses_only_queried_alt_and_accepts_vep_ampersand_terms(self):
        from clinical_variant_reporter import _extract_evidence_from_vep

        ev = _extract_evidence_from_vep(
            self._vep("G:benign;T:pathogenic&likely_pathogenic"), self._record("T"))

        assert ev.clinvar_assertion.terms == ("pathogenic", "likely_pathogenic")
        assert ev.clinvar_assertion.is_conflicting is False
        assert ev.clinvar_review_stars == 0

    def test_does_not_copy_a_different_alt_from_a_multiallelic_site(self):
        from clinical_variant_reporter import _extract_evidence_from_vep

        ev = _extract_evidence_from_vep(self._vep("G:benign"), self._record("T"))

        assert ev.clinvar_significance == ""
        assert ev.clinvar_review_stars == 0
        assert ev.clinvar_assertion.terms == ()

    @pytest.mark.parametrize("clin_sig_allele", [
        "T",                         # no allele-to-assertion separator
        "G:benign;T:",               # no assertion for the queried allele
        "G:benign;T:pathogenic:bad", # ambiguous assertion separator
        {"review_status_stars": 3},  # no allele-specific ClinVar assertion
    ])
    def test_missing_or_unsafe_allele_payload_fails_closed(self, clin_sig_allele):
        from clinical_variant_reporter import _extract_evidence_from_vep

        ev = _extract_evidence_from_vep(self._vep(clin_sig_allele), self._record("T"))

        assert ev.clinvar_significance == ""
        assert ev.clinvar_review_stars == 0
