from __future__ import annotations

import csv
import importlib.util
import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
import requests

SKILL_DIR = Path(__file__).resolve().parents[1]
PROJECT_ROOT = SKILL_DIR.parent.parent

sys.path.insert(0, str(SKILL_DIR))
sys.path.insert(0, str(PROJECT_ROOT))

import illumina_bridge  # noqa: E402
from illumina_bundle import (  # noqa: E402
    discover_bundle_artifacts,
    is_recognizable_illumina_bundle,
    parse_qc_metrics,
    parse_sample_sheet,
    summarize_sample_sheet,
)
from illumina_providers import ICAMetadataProvider  # noqa: E402

ORCHESTRATOR_DIR = PROJECT_ROOT / "skills" / "bio-orchestrator"
sys.path.insert(0, str(ORCHESTRATOR_DIR))
import orchestrator  # noqa: E402

_RUNNER_SPEC = importlib.util.spec_from_file_location("clawbio_runner", PROJECT_ROOT / "clawbio.py")
clawbio_runner = importlib.util.module_from_spec(_RUNNER_SPEC)
assert _RUNNER_SPEC.loader is not None
_RUNNER_SPEC.loader.exec_module(clawbio_runner)


DEMO_BUNDLE = SKILL_DIR / "demo_bundle"
ICA_PROJECT_ID = "00000000-0000-4000-8000-000000000074"
ICA_ANALYSIS_ID = "00000000-0000-4000-8000-000000000075"


@pytest.fixture(autouse=True)
def isolated_ica_environment(monkeypatch):
    """Ambient tenant credentials must not turn offline tests into live lookups."""
    monkeypatch.delenv("ILLUMINA_ICA_API_KEY", raising=False)
    monkeypatch.delenv("ILLUMINA_ICA_BASE_URL", raising=False)


@pytest.fixture
def ica_v3_payload():
    """Relevant Project/AnalysisV3 response fields, with invented identifiers."""
    return {
        "project": {"id": ICA_PROJECT_ID, "name": "Demo ICA Project", "active": True},
        "run": {
            "id": ICA_ANALYSIS_ID,
            "reference": "synthetic-analysis-74",
            "userReference": "DRAGEN Germline Demo",
            "status": "SUCCEEDED",
            "pipeline": {"code": "DRAGEN Germline 4.3"},
        },
    }


class FakeResponse:
    def __init__(self, payload, status_code=200):
        self.payload = payload
        self.status_code = status_code

    def raise_for_status(self):
        response = requests.Response()
        response.status_code = self.status_code
        response.url = "https://tenant.illumina.com/ica/rest/synthetic-private-id"
        response.raise_for_status()

    def json(self):
        if isinstance(self.payload, Exception):
            raise self.payload
        return self.payload


class FakeSession:
    def __init__(self, responses=None):
        self.headers = {}
        self.responses = list(responses or [])
        self.calls = []

    def get(self, url, timeout, headers=None):
        self.calls.append(
            {
                "url": url,
                "timeout": timeout,
                "headers": dict(headers or {}),
            }
        )
        if not self.responses:
            raise AssertionError(f"Unexpected GET request: {url}")
        response = self.responses.pop(0)
        if isinstance(response, Exception):
            raise response
        return response

BASESPACE_SAMPLE_SHEET = """[Header],
FileFormatVersion,2
RunName,Demo_Run
InstrumentPlatform,NovaSeq
AnalysisLocation,Cloud

[BCLConvert_Data]
Sample_ID,Index,Index2
TumorA_dna,AAAACCCC,GGGGTTTT
TumorB_rna,CCCCAAAA,TTTTGGGG

[Cloud_TSO500S_Data]
Sample_ID,Sample_Type,Pair_ID,Sample_Feature,Index_ID,Index,Index2
TumorA_dna,DNA,TumorA,HRD,UDP0001,AAAACCCC,GGGGTTTT
TumorB_rna,RNA,TumorB,,UDP0002,CCCCAAAA,TTTTGGGG

[Cloud_Data]
Sample_ID,ProjectName,LibraryName,LibraryPrepKitName,IndexAdapterKitName
TumorA_dna,DemoProject,TumorA_dna_AAAACCCC_GGGGTTTT,TSO500_v2,TSO500v2_ForwardOrientation
TumorB_rna,DemoProject,TumorB_rna_CCCCAAAA_TTTTGGGG,TSO500_v2,TSO500v2_ForwardOrientation
"""

METRICS_OUTPUT_TSV = """DRAGEN TruSight Oncology 500 v2.6.2 Analysis Software - Metrics Output

[Header]
Output Date\t2025-11-25
Output Time\t10:24:49
Workflow Version\t2.6.2.4

[Run QC Metrics]
Metric (UOM)\tLSL Guideline\tUSL Guideline\tValue
PCT_PF_READS (%)\t55.0\tNA\t77.5
PCT_Q30_R1 (%)\t80.0\tNA\t92.8
PCT_Q30_R2 (%)\t80.0\tNA\t92.3

[Analysis Status]
\tTumorA_dna\tTumorB_rna
COMPLETED_ALL_STEPS\tTRUE\tTRUE
FAILED_STEPS\tNA\tNA
STEPS_NOT_EXECUTED\tNA\tNA
"""


@pytest.fixture
def copied_bundle(tmp_path):
    dest = tmp_path / "bundle"
    shutil.copytree(DEMO_BUNDLE, dest)
    return dest


def test_discover_bundle_artifacts_success(copied_bundle):
    artifacts = discover_bundle_artifacts(copied_bundle)
    assert artifacts.sample_sheet_path.name == "SampleSheet.csv"
    assert artifacts.vcf_path.name == "demo.vcf"
    assert artifacts.qc_path.name == "qc_metrics.json"


def test_discover_bundle_artifacts_explicit_override_precedence(copied_bundle, tmp_path):
    alt_vcf = tmp_path / "manual.vcf"
    alt_vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    artifacts = discover_bundle_artifacts(copied_bundle, vcf_override=alt_vcf)
    assert artifacts.vcf_path == alt_vcf.resolve()


def test_discover_bundle_missing_sample_sheet_fails(copied_bundle):
    (copied_bundle / "SampleSheet.csv").unlink()
    with pytest.raises(FileNotFoundError):
        discover_bundle_artifacts(copied_bundle)


def test_discover_bundle_missing_vcf_fails(copied_bundle):
    (copied_bundle / "demo.vcf").unlink()
    with pytest.raises(FileNotFoundError):
        discover_bundle_artifacts(copied_bundle)


def test_parse_sample_sheet_extracts_rows():
    rows = parse_sample_sheet(DEMO_BUNDLE / "SampleSheet.csv")
    assert len(rows) == 2
    assert rows[0]["sample_id"] == "DEMO_SAMPLE_01"
    assert rows[1]["sample_project"] == "ClawBioDemo"


def test_parse_sample_sheet_merges_basespace_sections(tmp_path):
    sample_sheet = tmp_path / "samplesheet.csv"
    sample_sheet.write_text(BASESPACE_SAMPLE_SHEET, encoding="utf-8")
    rows = parse_sample_sheet(sample_sheet)
    assert len(rows) == 2
    assert rows[0]["sample_id"] == "TumorA_dna"
    assert rows[0]["sample_name"] == "TumorA"
    assert rows[0]["sample_type"] == "DNA"
    assert rows[0]["sample_feature"] == "HRD"
    assert rows[0]["library_name"] == "TumorA_dna_AAAACCCC_GGGGTTTT"
    assert rows[1]["sample_type"] == "RNA"
    assert rows[1]["sample_project"] == "DemoProject"


def test_parse_qc_metrics_json_normalizes_fixture():
    qc = parse_qc_metrics(DEMO_BUNDLE / "qc_metrics.json")
    assert qc["run_id"] == "demo-run-001"
    assert qc["percent_q30"] == 92.7
    assert qc["instrument"] == "NovaSeq X Plus"


def test_parse_qc_metrics_tsv_normalizes_metrics_output(tmp_path):
    metrics_tsv = tmp_path / "MetricsOutput.tsv"
    metrics_tsv.write_text(METRICS_OUTPUT_TSV, encoding="utf-8")
    qc = parse_qc_metrics(metrics_tsv)
    assert qc["analysis_software"] == "DRAGEN TruSight Oncology 500 v2.6.2 Analysis Software - Metrics Output"
    assert qc["workflow_version"] == "2.6.2.4"
    assert qc["percent_pf_reads"] == 77.5
    assert qc["percent_q30"] == 92.55
    assert qc["completed_samples"] == 2
    assert qc["reported_sample_count"] == 2


def test_parse_qc_metrics_malformed_json_raises(tmp_path):
    bad_qc = tmp_path / "bad_qc.json"
    bad_qc.write_text("{not valid json", encoding="utf-8")
    with pytest.raises(ValueError):
        parse_qc_metrics(bad_qc)


def test_discover_bundle_prefers_primary_result_vcf(tmp_path):
    bundle = tmp_path / "bundle"
    (bundle / "Results" / "TumorA" / "TumorA_dna").mkdir(parents=True)
    (bundle / "Logs_Intermediates" / "TumorA" / "TumorA_dna").mkdir(parents=True)
    (bundle / "samplesheet.csv").write_text(BASESPACE_SAMPLE_SHEET, encoding="utf-8")
    (bundle / "MetricsOutput.tsv").write_text(METRICS_OUTPUT_TSV, encoding="utf-8")
    preferred_vcf = bundle / "Results" / "TumorA" / "TumorA_dna" / "TumorA_dna.hard-filtered.vcf"
    preferred_vcf.write_text("##fileformat=VCFv4.2\n", encoding="utf-8")
    (bundle / "Results" / "TumorA" / "TumorA_dna" / "TumorA_dna.cnv.vcf").write_text(
        "##fileformat=VCFv4.2\n",
        encoding="utf-8",
    )
    (bundle / "Logs_Intermediates" / "TumorA" / "TumorA_dna" / "TumorA_dna.hard-filtered.vcf.gz").write_text(
        "placeholder",
        encoding="utf-8",
    )

    artifacts = discover_bundle_artifacts(bundle)
    assert artifacts.vcf_path == preferred_vcf.resolve()


def test_build_summary_and_data_is_deterministic(copied_bundle):
    artifacts = discover_bundle_artifacts(copied_bundle)
    sample_rows = parse_sample_sheet(artifacts.sample_sheet_path)
    sample_summary = summarize_sample_sheet(sample_rows)
    qc_summary = parse_qc_metrics(artifacts.qc_path)
    provider = ICAMetadataProvider(api_key="")
    metadata = provider.enrich(
        bundle_dir=copied_bundle,
        project_id="ica-project-demo",
        run_id="ica-run-demo",
        allow_mock=True,
    )
    merged_rows, merge = illumina_bridge.merge_sample_metadata(sample_rows, metadata)
    hints = illumina_bridge.build_downstream_routing_hints(
        vcf_path=artifacts.vcf_path,
        sample_count=sample_summary["sample_count"],
    )
    summary1, data1 = illumina_bridge.build_summary_and_data(
        bundle=artifacts,
        sample_rows=merged_rows,
        sample_summary=sample_summary,
        qc_summary=qc_summary,
        metadata_result=metadata,
        metadata_merge=merge,
        downstream_hints=hints,
    )
    summary2, data2 = illumina_bridge.build_summary_and_data(
        bundle=artifacts,
        sample_rows=merged_rows,
        sample_summary=sample_summary,
        qc_summary=qc_summary,
        metadata_result=metadata,
        metadata_merge=merge,
        downstream_hints=hints,
    )
    assert summary1 == summary2
    assert data1 == data2


@pytest.mark.parametrize("active", [True, False])
def test_ica_v3_project_and_analysis_fields(copied_bundle, ica_v3_payload, active):
    ica_v3_payload["project"]["active"] = active
    session = FakeSession([FakeResponse(ica_v3_payload["project"]), FakeResponse(ica_v3_payload["run"])])
    result = ICAMetadataProvider(api_key="test-key", session=session).enrich(
        bundle_dir=copied_bundle, project_id=ICA_PROJECT_ID, run_id=ICA_ANALYSIS_ID,
    )
    assert result.status == "enriched"
    assert result.project == {
        "id": ICA_PROJECT_ID, "name": "Demo ICA Project", "active": active, "status": "",
    }
    assert result.run["name"] == "DRAGEN Germline Demo"
    assert result.run["status"] == "SUCCEEDED"
    assert result.run["pipeline"] == "DRAGEN Germline 4.3"
    assert [call["url"] for call in session.calls] == [
        f"https://ica.illumina.com/ica/rest/api/projects/{ICA_PROJECT_ID}",
        f"https://ica.illumina.com/ica/rest/api/projects/{ICA_PROJECT_ID}/analyses/{ICA_ANALYSIS_ID}",
    ]
    assert session.headers["Accept"] == "application/vnd.illumina.v3+json"


@pytest.mark.parametrize("analysis_status", [
    "REQUESTED", "AWAITINGINPUT", "INPROGRESS", "FAILED", "FAILEDFINAL", "ABORTED",
])
def test_ica_non_succeeded_analysis_warns_without_changing_enriched_status(
    copied_bundle, ica_v3_payload, analysis_status,
):
    ica_v3_payload["run"]["status"] = analysis_status
    session = FakeSession([FakeResponse(ica_v3_payload["project"]), FakeResponse(ica_v3_payload["run"])])
    result = ICAMetadataProvider(api_key="test-key", session=session).enrich(
        bundle_dir=copied_bundle, project_id=ICA_PROJECT_ID, run_id=ICA_ANALYSIS_ID,
    )
    assert result.status == "enriched"
    assert result.run["status"] == analysis_status
    assert any(analysis_status in warning and "not SUCCEEDED" in warning for warning in result.warnings)


def test_ica_reference_fallback_and_missing_active(copied_bundle, ica_v3_payload):
    del ica_v3_payload["run"]["userReference"]
    del ica_v3_payload["project"]["active"]
    session = FakeSession([FakeResponse(ica_v3_payload["project"]), FakeResponse(ica_v3_payload["run"])])
    result = ICAMetadataProvider(api_key="test-key", session=session).enrich(
        bundle_dir=copied_bundle, project_id=ICA_PROJECT_ID, run_id=ICA_ANALYSIS_ID,
    )
    assert result.run["name"] == "synthetic-analysis-74"
    assert result.project["active"] is None


def test_ica_existing_status_and_pipeline_name_remain_compatible(copied_bundle, ica_v3_payload):
    ica_v3_payload["project"]["status"] = "legacy-status"
    ica_v3_payload["run"]["pipeline"]["name"] = "Legacy pipeline name"
    session = FakeSession([FakeResponse(ica_v3_payload["project"]), FakeResponse(ica_v3_payload["run"])])
    result = ICAMetadataProvider(api_key="test-key", session=session).enrich(
        bundle_dir=copied_bundle, project_id=ICA_PROJECT_ID, run_id=ICA_ANALYSIS_ID,
    )
    assert result.project["status"] == "legacy-status"
    assert result.run["pipeline"] == "Legacy pipeline name"


def test_ica_v3_does_not_infer_sample_metadata(copied_bundle, ica_v3_payload):
    # An uncontracted legacy field must not silently enable sample-name matching.
    ica_v3_payload["run"]["samples"] = [{
        "sample_id": "DEMO_SAMPLE_01", "id": "unverified-sample", "status": "SUCCEEDED",
    }]
    session = FakeSession([FakeResponse(ica_v3_payload["project"]), FakeResponse(ica_v3_payload["run"])])
    result = ICAMetadataProvider(api_key="test-key", session=session).enrich(
        bundle_dir=copied_bundle, project_id=ICA_PROJECT_ID, run_id=ICA_ANALYSIS_ID,
    )
    assert result.samples == []
    assert any("sample-level" in warning and "unavailable" in warning for warning in result.warnings)
    assert not any("not SUCCEEDED" in warning for warning in result.warnings)
    assert len(session.calls) == 2


def test_ica_mock_contains_v3_metadata_without_sample_enrichment(copied_bundle):
    result = ICAMetadataProvider(api_key="").enrich(
        bundle_dir=copied_bundle, project_id=ICA_PROJECT_ID, run_id=ICA_ANALYSIS_ID, allow_mock=True,
    )
    assert result.status == "mocked-demo"
    assert result.project["active"] is True
    assert result.project["status"] == ""
    assert result.run["name"] == "DRAGEN Germline Demo"
    assert result.samples == []
    assert any("sample-level" in warning and "unavailable" in warning for warning in result.warnings)


def test_ica_mock_keeps_local_samples_without_remote_annotations(copied_bundle):
    provider = ICAMetadataProvider(api_key="")
    result = provider.enrich(
        bundle_dir=copied_bundle,
        project_id="ica-project-demo",
        run_id="ica-run-demo",
        allow_mock=True,
    )
    rows = parse_sample_sheet(copied_bundle / "SampleSheet.csv")
    merged_rows, merge = illumina_bridge.merge_sample_metadata(rows, result)
    assert result.status == "mocked-demo"
    assert result.project["name"] == "Demo ICA Project"
    assert result.run["status"] == "SUCCEEDED"
    assert merge == {"samples_in_bundle": 2, "samples_enriched": 0, "samples_unmatched": 2}
    assert merged_rows[0]["sample_id"] == "DEMO_SAMPLE_01"
    assert merged_rows[0]["ica_sample_id"] == ""


def test_ica_provider_missing_api_key_yields_warning(copied_bundle):
    session = FakeSession()
    provider = ICAMetadataProvider(api_key="", session=session)
    assert provider.api_key is None
    result = provider.enrich(
        bundle_dir=copied_bundle,
        project_id="ica-project-demo",
        run_id="ica-run-demo",
        allow_mock=False,
    )
    assert result.status == "warning"
    assert "ILLUMINA_ICA_API_KEY" in result.warnings[0]
    assert session.calls == []


@pytest.mark.parametrize("project_id, analysis_id", [
    (None, ICA_ANALYSIS_ID), (ICA_PROJECT_ID, None), (None, None),
])
def test_ica_missing_ids_explain_flags_without_requests(copied_bundle, project_id, analysis_id):
    session = FakeSession()
    result = ICAMetadataProvider(api_key="test-key", session=session).enrich(
        bundle_dir=copied_bundle, project_id=project_id, run_id=analysis_id,
    )
    assert result.status == "skipped"
    assert "--ica-project-id" in result.warnings[0]
    assert "--ica-run-id" in result.warnings[0]
    assert "analysis ID" in result.warnings[0]
    assert session.calls == []


@pytest.mark.parametrize("lookup", ["project", "analysis"])
@pytest.mark.parametrize("status_code, hint", [
    (401, "ILLUMINA_ICA_API_KEY"),
    (403, "permissions"),
    (404, "--ica-"),
    (429, "retry"),
    (500, "retry"),
    (503, "retry"),
    (418, "access"),
])
def test_ica_http_failure_identifies_lookup_and_keeps_local_import(
    tmp_path, monkeypatch, ica_v3_payload, lookup, status_code, hint,
):
    responses = [FakeResponse({}, status_code=status_code)]
    if lookup == "analysis":
        responses.insert(0, FakeResponse(ica_v3_payload["project"]))
    session = FakeSession(responses)
    monkeypatch.setenv("ILLUMINA_ICA_API_KEY", "test-key")
    monkeypatch.setattr(requests, "Session", lambda: session)
    output_dir = tmp_path / "import"
    result = illumina_bridge.import_bundle(
        bundle_dir=DEMO_BUNDLE, output_dir=output_dir, metadata_provider_name="ica",
        ica_project_id=ICA_PROJECT_ID, ica_run_id=ICA_ANALYSIS_ID,
    )
    metadata = result["data"]["metadata_enrichment"]
    warning = " ".join(metadata["warnings"])
    assert metadata["status"] == "warning"
    assert f"ICA {lookup} lookup" in warning
    assert f"HTTP {status_code}" in warning
    assert hint.lower() in warning.lower()
    if status_code == 404:
        assert ("--ica-project-id" if lookup == "project" else "--ica-run-id") in warning
    assert "synthetic-private-id" not in warning
    assert "https://" not in warning
    assert len(session.calls) == (1 if lookup == "project" else 2)
    assert result["summary"]["sample_count"] == 2
    assert result["summary"]["metadata_status"] == "warning"
    assert metadata["merge"]["samples_enriched"] == 0
    assert (output_dir / "report.md").exists()
    assert (output_dir / "result.json").exists()
    assert (output_dir / "tables" / "sample_manifest.csv").exists()
    assert (output_dir / "reproducibility" / "commands.sh").exists()
    assert warning in (output_dir / "report.md").read_text()


@pytest.mark.parametrize("lookup", ["project", "analysis"])
@pytest.mark.parametrize("error, hint", [
    (requests.Timeout("synthetic-private-id"), "timed out"),
    (requests.ConnectionError("synthetic-private-id"), "network"),
    (requests.RequestException("synthetic-private-id"), "request failed"),
    (ValueError("synthetic-private-id"), "JSON"),
])
def test_ica_transport_failure_has_actionable_warning(
    copied_bundle, ica_v3_payload, lookup, error, hint,
):
    responses = [FakeResponse(error) if isinstance(error, ValueError) else error]
    if lookup == "analysis":
        responses.insert(0, FakeResponse(ica_v3_payload["project"]))
    session = FakeSession(responses)
    provider = ICAMetadataProvider(api_key="test-key", session=session)
    result = provider.enrich(
        bundle_dir=copied_bundle,
        project_id=ICA_PROJECT_ID,
        run_id=ICA_ANALYSIS_ID,
    )
    assert result.status == "warning"
    warning = " ".join(result.warnings)
    assert f"ICA {lookup} lookup" in warning
    assert hint.lower() in warning.lower()
    assert "synthetic-private-id" not in warning
    assert len(session.calls) == (1 if lookup == "project" else 2)


@pytest.mark.parametrize("payload", [None, [], "not an object"])
def test_ica_non_object_response_warns_instead_of_crashing(copied_bundle, payload):
    provider = ICAMetadataProvider(api_key="test-key", session=FakeSession([FakeResponse(payload)]))
    result = provider.enrich(bundle_dir=copied_bundle, project_id=ICA_PROJECT_ID, run_id=ICA_ANALYSIS_ID)
    assert result.status == "warning"
    assert "ICA project lookup" in result.warnings[0]
    assert "JSON" in result.warnings[0]


def test_ica_provider_does_not_store_api_key_in_session_headers():
    session = FakeSession()
    provider = ICAMetadataProvider(api_key="super-secret-key", session=session)
    assert provider.api_key == "super-secret-key"
    assert session.headers["Accept"] == "application/vnd.illumina.v3+json"
    assert "X-API-Key" not in session.headers


def test_ica_provider_fetch_json_sends_api_key_per_request():
    session = FakeSession([FakeResponse({"id": "project-1"})])
    provider = ICAMetadataProvider(api_key="super-secret-key", session=session)
    payload = provider._fetch_json("/api/projects/project-1")
    assert payload["id"] == "project-1"
    assert session.calls[0]["headers"]["X-API-Key"] == "super-secret-key"


def test_ica_provider_repr_redacts_api_key():
    provider = ICAMetadataProvider(api_key="super-secret-key")
    rendered = repr(provider)
    assert "super-secret-key" not in rendered
    assert "has_api_key=True" in rendered


def test_ica_provider_invalid_base_url_falls_back_with_warning(copied_bundle, ica_v3_payload):
    session = FakeSession(
        [
            FakeResponse(ica_v3_payload["project"]),
            FakeResponse(ica_v3_payload["run"]),
        ]
    )
    provider = ICAMetadataProvider(
        api_key="test-key",
        base_url="http://example.com/ica/rest",
        session=session,
    )
    assert provider.base_url == "https://ica.illumina.com/ica/rest"
    result = provider.enrich(
        bundle_dir=copied_bundle,
        project_id="project-1",
        run_id="run-1",
    )
    assert result.status == "enriched"
    assert any("ILLUMINA_ICA_BASE_URL" in warning for warning in result.warnings)


def test_ica_provider_accepts_trusted_illumina_base_url():
    session = FakeSession()
    provider = ICAMetadataProvider(
        api_key="test-key",
        base_url="https://tenant.illumina.com/ica/rest",
        session=session,
    )
    assert provider.base_url == "https://tenant.illumina.com/ica/rest"
    assert provider._initialization_warnings == []


@pytest.mark.parametrize("analysis_status", ["SUCCEEDED", "REQUESTED", "FAILED"])
def test_import_report_distinguishes_enrichment_from_analysis_status(
    tmp_path, monkeypatch, ica_v3_payload, analysis_status,
):
    ica_v3_payload["run"]["status"] = analysis_status
    session = FakeSession([FakeResponse(ica_v3_payload["project"]), FakeResponse(ica_v3_payload["run"])])
    monkeypatch.setenv("ILLUMINA_ICA_API_KEY", "test-key")
    monkeypatch.setattr(requests, "Session", lambda: session)
    output_dir = tmp_path / "enriched"
    illumina_bridge.import_bundle(
        bundle_dir=DEMO_BUNDLE, output_dir=output_dir, metadata_provider_name="ica",
        ica_project_id=ICA_PROJECT_ID, ica_run_id=ICA_ANALYSIS_ID,
    )
    payload = json.loads((output_dir / "result.json").read_text())
    metadata = payload["data"]["metadata_enrichment"]
    report = (output_dir / "report.md").read_text()
    assert payload["summary"]["metadata_status"] == "enriched"
    assert metadata["status"] == "enriched"
    assert metadata["project"]["active"] is True
    assert metadata["project"]["status"] == ""
    assert metadata["run"]["name"] == "DRAGEN Germline Demo"
    assert metadata["run"]["status"] == analysis_status
    assert metadata["samples"] == []
    assert metadata["merge"]["samples_enriched"] == 0
    assert "**Analysis**: DRAGEN Germline Demo" in report
    assert f"**Analysis status**: {analysis_status}" in report
    assert "sample-level enrichment is unavailable" in report
    for call in session.calls:
        assert call["headers"] == {"X-API-Key": "test-key"}


def test_ica_run_id_help_describes_analysis_id(capsys):
    with pytest.raises(SystemExit) as exc_info:
        illumina_bridge.parse_args(["--help"])
    assert exc_info.value.code == 0
    help_text = " ".join(capsys.readouterr().out.split())
    assert "ICA analysis ID" in help_text
    assert "not a sequencing run ID" in help_text


def test_import_bundle_demo_creates_standard_outputs(tmp_path):
    output_dir = tmp_path / "demo_output"
    result = illumina_bridge.import_bundle(
        bundle_dir=DEMO_BUNDLE,
        output_dir=output_dir,
        metadata_provider_name="none",
        allow_mock_metadata=True,
    )
    assert result["summary"]["platform"] == "illumina"
    assert (output_dir / "report.md").exists()
    assert (output_dir / "result.json").exists()
    assert (output_dir / "tables" / "sample_manifest.csv").exists()
    assert (output_dir / "reproducibility" / "commands.sh").exists()


def test_import_bundle_with_basespace_style_inputs_creates_standard_outputs(tmp_path):
    bundle = tmp_path / "basespace_bundle"
    (bundle / "Results" / "TumorA" / "TumorA_dna").mkdir(parents=True)
    (bundle / "samplesheet.csv").write_text(BASESPACE_SAMPLE_SHEET, encoding="utf-8")
    (bundle / "MetricsOutput.tsv").write_text(METRICS_OUTPUT_TSV, encoding="utf-8")
    (bundle / "Results" / "TumorA" / "TumorA_dna" / "TumorA_dna.hard-filtered.vcf").write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
        encoding="utf-8",
    )

    output_dir = tmp_path / "basespace_output"
    result = illumina_bridge.import_bundle(
        bundle_dir=bundle,
        output_dir=output_dir,
        metadata_provider_name="none",
        allow_mock_metadata=False,
    )
    assert result["summary"]["platform"] == "illumina"
    assert result["summary"]["sample_count"] == 2
    assert (output_dir / "report.md").exists()
    manifest_rows = list(csv.DictReader((output_dir / "tables" / "sample_manifest.csv").open(encoding="utf-8")))
    assert manifest_rows[0]["sample_type"] == "DNA"
    assert manifest_rows[0]["library_name"] == "TumorA_dna_AAAACCCC_GGGGTTTT"


def test_clawbio_run_illumina_input_bundle_completes(tmp_path):
    output_dir = tmp_path / "cli_output"
    proc = subprocess.run(
        [
            sys.executable,
            str(PROJECT_ROOT / "clawbio.py"),
            "run",
            "illumina",
            "--input",
            str(DEMO_BUNDLE),
            "--output",
            str(output_dir),
        ],
        capture_output=True,
        text=True,
        cwd=str(PROJECT_ROOT),
        check=False,
    )
    assert proc.returncode == 0, proc.stderr
    assert (output_dir / "report.md").exists()
    payload = json.loads((output_dir / "result.json").read_text())
    assert payload["data"]["platform"] == "illumina"


def test_clawbio_run_illumina_forwards_ica_flags_offline(tmp_path):
    output_dir = tmp_path / "ica_cli_output"
    proc = subprocess.run(
        [sys.executable, str(PROJECT_ROOT / "clawbio.py"), "run", "illumina", "--demo",
         "--metadata-provider", "ica", "--ica-project-id", ICA_PROJECT_ID,
         "--ica-run-id", ICA_ANALYSIS_ID, "--output", str(output_dir)],
        capture_output=True, text=True, cwd=str(PROJECT_ROOT), check=False,
    )
    assert proc.returncode == 0, proc.stderr
    payload = json.loads((output_dir / "result.json").read_text())
    metadata = payload["data"]["metadata_enrichment"]
    assert payload["summary"]["metadata_status"] == "mocked-demo"
    assert metadata["project"]["active"] is True
    assert metadata["run"]["name"] == "DRAGEN Germline Demo"
    assert metadata["samples"] == []
    commands = (output_dir / "reproducibility" / "commands.sh").read_text()
    assert "\nILLUMINA_ICA_API_KEY= python clawbio.py run illumina --demo" in commands
    assert f"--ica-project-id {ICA_PROJECT_ID}" in commands
    assert f"--ica-run-id {ICA_ANALYSIS_ID}" in commands


def test_orchestrator_routes_illumina_keywords():
    skill, _ = orchestrator.detect_skill_with_hint_from_query(
        "Import this Illumina DRAGEN sample sheet bundle and add ICA metadata"
    )
    assert skill == "illumina-bridge"


def test_orchestrator_routes_illumina_bundle_directory():
    assert is_recognizable_illumina_bundle(DEMO_BUNDLE)
    assert orchestrator.detect_skill_from_file(DEMO_BUNDLE) == "illumina-bridge"


def test_security_filter_rejects_unsupported_flags_for_illumina(tmp_path):
    output_dir = tmp_path / "secure_output"
    result = clawbio_runner.run_skill(
        skill_name="illumina",
        demo=True,
        output_dir=str(output_dir),
        extra_args=["--bogus", "nope", "--metadata-provider", "none"],
    )
    assert result["success"] is True
    assert (output_dir / "result.json").exists()


def test_sample_manifest_csv_contains_expected_columns(tmp_path):
    output_dir = tmp_path / "manifest_output"
    result = illumina_bridge.import_bundle(
        bundle_dir=DEMO_BUNDLE,
        output_dir=output_dir,
        metadata_provider_name="ica",
        ica_project_id="ica-project-demo",
        ica_run_id="ica-run-demo",
        allow_mock_metadata=True,
    )
    manifest = output_dir / "tables" / "sample_manifest.csv"
    with manifest.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert rows[0]["sample_id"] == "DEMO_SAMPLE_01"
    assert result["summary"]["metadata_status"] == "mocked-demo"
    for row in rows:
        for field in ("ica_sample_id", "ica_analysis_status", "ica_cohort", "ica_notes"):
            assert row[field] == ""
