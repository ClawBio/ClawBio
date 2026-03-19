"""Tests for the variant-annotation skill."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Any, cast

import pytest
import requests

_PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent.parent
_SKILL_DIR = Path(__file__).resolve().parent.parent

if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))
if str(_SKILL_DIR) not in sys.path:
    sys.path.insert(0, str(_SKILL_DIR))

import variant_annotation as va


DEMO_VCF = _SKILL_DIR / "example_data" / "synthetic_clinvar_panel.vcf"


@pytest.fixture
def parsed_demo_records() -> list[dict]:
    return va.parse_vcf(DEMO_VCF)


@pytest.fixture
def record_template() -> dict:
    return {
        "input_id": "rsTest",
        "chrom": "1",
        "pos": 100,
        "ref": "G",
        "alt": "A",
        "genotype": "G/A",
    }


@pytest.fixture
def vep_entry() -> dict:
    return {
        "id": "rsTest",
        "input": "1 100 100 G/A +",
        "variant_class": "SNV",
        "most_severe_consequence": "missense_variant",
        "warnings": ["pick_order_warning"],
        "transcript_consequences": [
            {
                "transcript_id": "ENST000001",
                "gene_symbol": "GENE1",
                "gene_id": "ENSG000001",
                "consequence_terms": ["synonymous_variant"],
                "impact": "LOW",
                "canonical": 1,
                "hgvsc": "ENST000001:c.100G>A",
                "hgvsp": "ENSP000001:p.Gly34Ser",
            },
            {
                "transcript_id": "ENST000002",
                "gene_symbol": "GENE1",
                "gene_id": "ENSG000001",
                "consequence_terms": ["stop_gained", "missense_variant"],
                "impact": "HIGH",
                "mane_select": "NM_000001.1",
                "hgvsc": "ENST000002:c.100G>A",
                "hgvsp": "ENSP000002:p.Trp34Ter",
            },
        ],
        "colocated_variants": [
            {
                "id": "rsTest",
                "clin_sig": ["pathogenic", "likely_pathogenic"],
                "clinvar_accession": "VCV000000001",
                "frequencies": {
                    "A": {
                        "gnomad": 0.0005,
                        "gnomad_afr": 0.0008,
                        "1000GENOMES:phase_3:ALL": 0.004,
                    }
                },
                "gnomad_af": 0.0005,
                "afr_af": 0.02,
            }
        ],
    }


class DummyResponse:
    def __init__(
        self, payload: list[dict], status_code: int = 200, headers: dict | None = None
    ):
        self._payload = payload
        self.status_code = status_code
        self.headers = headers or {}

    def json(self) -> list[dict]:
        return self._payload

    def raise_for_status(self) -> None:
        if self.status_code >= 400:
            raise requests.HTTPError(f"status {self.status_code}")


class DummySession:
    def __init__(self, responses: list[DummyResponse]):
        self.responses = list(responses)
        self.headers: dict[str, str] = {}
        self.calls: list[dict] = []

    def post(
        self,
        url: str,
        params: dict | None = None,
        json: dict | None = None,
        timeout: int | None = None,
    ):
        self.calls.append(
            {"url": url, "params": params, "json": json, "timeout": timeout}
        )
        return self.responses.pop(0)


def test_parse_vcf_reads_demo_variants(parsed_demo_records):
    assert len(parsed_demo_records) == 20

    ids = {record["input_id"] for record in parsed_demo_records}
    assert {"rs1801133", "rs1799853", "rs1057910", "rs9923231"}.issubset(ids)

    vkorc1 = next(
        record for record in parsed_demo_records if record["input_id"] == "rs9923231"
    )
    apoe = next(
        record for record in parsed_demo_records if record["input_id"] == "rs7412"
    )
    assert vkorc1["genotype"] == "A/A"
    assert apoe["genotype"] == "C/C"


def test_parse_vcf_raises_without_pysam(monkeypatch):
    monkeypatch.setattr(va, "pysam", None)
    with pytest.raises(ImportError, match="pysam"):
        va.parse_vcf(DEMO_VCF)


def test_batch_variants_and_region_string(parsed_demo_records):
    batches = va.batch_variants(parsed_demo_records, batch_size=7)
    assert [len(batch) for batch in batches] == [7, 7, 6]

    vep_string = va.build_vep_variant_string(parsed_demo_records[0])
    assert vep_string == "1 11794419 11794419 T/A +"


def test_annotation_helpers_prioritize_severity_and_extract_fields(
    record_template, vep_entry
):
    best = va.choose_best_annotation(vep_entry)
    assert best["consequence"] == "stop_gained"
    assert best["impact"] == "HIGH"
    assert best["feature_id"] == "ENST000002"

    clinvar = va.extract_clinvar_fields(vep_entry)
    assert clinvar["clinvar_significance"] == "pathogenic;likely_pathogenic"
    assert clinvar["clinvar_accessions"] == "VCV000000001;rsTest"

    freqs = va.extract_frequency_fields(vep_entry, "A")
    assert freqs["gnomad_af"] == pytest.approx(0.0005)
    assert freqs["global_af"] == pytest.approx(0.0005)
    assert freqs["max_af"] == pytest.approx(0.02)
    assert freqs["highest_frequency_population"] == "afr_af"
    assert freqs["population_outlier_flag"] is False

    annotation = va.normalize_annotation(record_template, vep_entry)
    assert annotation["gene"] == "GENE1"
    assert annotation["impact_tier"] == "HIGH"
    assert annotation["clinvar_bucket"] == "pathogenic_or_likely_pathogenic"
    assert annotation["is_rare"] is True
    assert annotation["clinically_relevant"] is True
    assert annotation["priority_bucket"] == "Tier 1"
    assert annotation["priority_score"] > 0
    assert annotation["warning"] == "pick_order_warning"


def test_pathogenic_and_list_helpers():
    assert va.ensure_list(None) == []
    assert va.ensure_list(("a", "b")) == ["a", "b"]
    assert va.join_unique(["x", "x", "y", ""]) == "x;y"
    assert va.to_float("0.5") == pytest.approx(0.5)
    assert va.to_float(".") is None
    assert va.has_pathogenic_clinvar("benign/pathogenic") is True
    assert va.has_pathogenic_clinvar("drug_response") is False
    assert va.classify_clinvar_bucket("drug_response") == "drug_response"
    assert va.classify_clinvar_bucket("Conflicting classifications") == "conflicting"
    assert va.impact_tier("HIGH") == "HIGH"
    assert va.impact_tier("") == "UNKNOWN"
    assert (
        va.is_clinically_relevant(
            {"gnomad_af": 0.01, "clinvar_significance": "pathogenic"}
        )
        is False
    )


def test_choose_best_annotation_falls_back_to_top_level():
    entry = {"most_severe_consequence": "intergenic_variant", "impact": "MODIFIER"}
    best = va.choose_best_annotation(entry)
    assert best["consequence"] == "intergenic_variant"
    assert best["impact"] == "MODIFIER"
    assert best["gene"] == ""


def test_annotate_variants_success_missing_and_error(record_template, vep_entry):
    records = [
        record_template,
        {**record_template, "input_id": "rsMissing", "pos": 101},
        {**record_template, "input_id": "rsError", "pos": 102},
    ]

    class FakeClient:
        def __init__(self):
            self.calls = 0

        def annotate_batch(self, variant_strings, assembly="GRCh38"):
            self.calls += 1
            if self.calls == 1:
                return [vep_entry]
            raise RuntimeError("backend unavailable")

    annotations, failures, metadata = va.annotate_variants(
        records,
        cast(Any, FakeClient()),
        batch_size=2,
    )
    assert len(annotations) == 3
    assert len(failures) == 2
    assert metadata["batches_sent"] == 2
    assert metadata["failed_batches"] == 1

    statuses = {row["input_id"]: row["vep_status"] for row in annotations}
    assert statuses["rsTest"] == "ok"
    assert statuses["rsMissing"] == "missing"
    assert statuses["rsError"] == "error"


def test_vep_client_caches_and_retries(tmp_path, monkeypatch):
    responses = [
        DummyResponse([], status_code=429, headers={"Retry-After": "0"}),
        DummyResponse([{"id": "rs1"}], status_code=200),
    ]
    session = DummySession(responses)
    client = va.VEPClient(cache_dir=tmp_path, use_cache=True, rate_interval=0)
    client.session = session  # type: ignore[assignment]

    sleep_calls = []
    monkeypatch.setattr(va.time, "sleep", lambda seconds: sleep_calls.append(seconds))

    first = client.annotate_batch(["1 100 100 G/A +"])
    second = client.annotate_batch(["1 100 100 G/A +"])

    assert first == [{"id": "rs1"}]
    assert second == first
    assert len(session.calls) == 2
    assert sleep_calls == [0.0]
    assert any(tmp_path.iterdir())


def test_write_outputs_and_report_generation(tmp_path, record_template, vep_entry):
    annotation = va.normalize_annotation(record_template, vep_entry)
    tsv_path = va.write_tsv(tmp_path, [annotation])
    report = va.generate_markdown_report(
        input_path=DEMO_VCF,
        output_dir=tmp_path,
        annotations=[annotation],
        failures=[],
        metadata={
            "assembly": "GRCh38",
            "batch_size": 200,
            "batches_sent": 1,
            "failed_batches": 0,
        },
        tsv_path=tsv_path,
    )

    assert tsv_path.exists()
    with tsv_path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    assert reader.fieldnames == va.ANNOTATED_COLUMNS
    assert len(rows) == 1
    assert rows[0]["input_id"] == "rsTest"
    assert rows[0]["gene"] == "GENE1"
    assert rows[0]["priority_bucket"] == "Tier 1"
    assert rows[0]["clinvar_bucket"] == "pathogenic_or_likely_pathogenic"
    assert rows[0]["clinically_relevant"] == "true"
    assert "Variant Annotation Report" in report
    assert "Top Prioritized Variants" in report
    assert "Clinically Relevant Findings" in report
    assert "Impact Tiers" in report
    assert "research and educational tool" in report


def test_write_reproducibility_bundle(tmp_path):
    args = argparse.Namespace(
        demo=False,
        assembly="GRCh38",
        batch_size=123,
        no_cache=False,
        cache_dir=str(tmp_path / "cache"),
    )
    commands_path = va.write_reproducibility_bundle(tmp_path, args, DEMO_VCF)
    text = commands_path.read_text()
    assert "--input" in text
    assert str(DEMO_VCF) in text
    assert "--batch-size 123" in text
    assert "--cache-dir" in text


def test_run_pipeline_writes_expected_files(tmp_path, monkeypatch):
    flagged_ids = {"rs1801133", "rs6025", "rs334"}

    def make_annotation(record: dict) -> dict:
        is_flagged = record["input_id"] in flagged_ids
        clinvar_bucket = (
            "pathogenic_or_likely_pathogenic" if is_flagged else "drug_response"
        )
        reasons = "rare_pathogenic_clinvar" if is_flagged else "drug_response"
        return {
            "input_id": record["input_id"],
            "chrom": record["chrom"],
            "pos": record["pos"],
            "ref": record["ref"],
            "alt": record["alt"],
            "genotype": record["genotype"],
            "gene": f"GENE_{record['input_id']}",
            "gene_id": f"ENSG_{record['input_id']}",
            "feature_type": "transcript",
            "feature_id": f"ENST_{record['input_id']}",
            "consequence": "missense_variant",
            "all_consequences": "missense_variant;synonymous_variant",
            "impact": "MODERATE",
            "impact_tier": "MODERATE",
            "consequence_severity_rank": 12,
            "hgvsc": f"c.{record['pos']}{record['ref']}>{record['alt']}",
            "hgvsp": "p.MockVariant",
            "clinvar_significance": "pathogenic" if is_flagged else "drug_response",
            "clinvar_bucket": clinvar_bucket,
            "clinvar_accessions": f"VCV_{record['input_id']};{record['input_id']}",
            "is_rare": is_flagged,
            "gnomad_af": 0.0002 if is_flagged else 0.02,
            "global_af": 0.0002 if is_flagged else 0.02,
            "max_af": 0.0008 if is_flagged else 0.05,
            "min_af": 0.0001 if is_flagged else 0.01,
            "highest_frequency_population": "gnomad_nfe" if is_flagged else "gnomad",
            "population_frequency_spread": 0.001 if is_flagged else 0.02,
            "population_freq_summary": "gnomad_nfe=0.0008;gnomad=0.0002"
            if is_flagged
            else "gnomad=0.02;gnomad_afr=0.01",
            "population_outlier_flag": False,
            "variant_class": "SNV",
            "existing_variation": record["input_id"],
            "priority_score": 92 if is_flagged else 41,
            "priority_bucket": "Tier 1" if is_flagged else "Tier 3",
            "review_reasons": reasons,
            "clinically_relevant": is_flagged,
            "vep_status": "ok",
            "warning": "",
        }

    def fake_annotate_variants(records, client, assembly="GRCh38", batch_size=200):
        annotations = [make_annotation(record) for record in records]
        failures = []
        metadata = {
            "assembly": assembly,
            "batch_size": batch_size,
            "batches_sent": 1,
            "failed_batches": 0,
        }
        return annotations, failures, metadata

    monkeypatch.setattr(va, "annotate_variants", fake_annotate_variants)

    args = argparse.Namespace(
        input=str(DEMO_VCF),
        output=str(tmp_path / "out"),
        demo=False,
        assembly="GRCh38",
        batch_size=50,
        cache_dir=str(tmp_path / "cache"),
        no_cache=True,
    )
    result = va.run_pipeline(args)

    output_dir = Path(result["output_dir"])
    assert (output_dir / "report.md").exists()
    assert (output_dir / "tables" / "annotated_variants.tsv").exists()
    assert (output_dir / "reproducibility" / "commands.sh").exists()

    payload = json.loads((output_dir / "result.json").read_text())
    assert payload["skill"] == "variant-annotation"
    assert payload["summary"]["variants_in_vcf"] == 20
    assert payload["summary"]["variants_annotated"] == 20
    assert payload["summary"]["clinically_relevant_findings"] == 3
    assert payload["summary"]["tier_1_variants"] == 3
    assert payload["data"]["files"]["annotated_tsv"] == "tables/annotated_variants.tsv"
    assert len(payload["data"]["top_variants"]) > 0
    assert payload["data"]["summary_tables"]["priority_buckets"]["Tier 1"] == 3
    assert {
        row["input_id"] for row in payload["data"]["flagged_findings"]
    } == flagged_ids

    with (output_dir / "tables" / "annotated_variants.tsv").open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    assert reader.fieldnames == va.ANNOTATED_COLUMNS
    assert len(rows) == 20
    assert rows[0]["priority_bucket"] in {"Tier 1", "Tier 3"}
    assert {
        row["input_id"] for row in rows if row["clinically_relevant"] == "true"
    } == flagged_ids


def test_main_exits_without_input(monkeypatch, capsys):
    monkeypatch.setattr(
        va,
        "parse_args",
        lambda: argparse.Namespace(input=None, demo=False, batch_size=10),
    )
    with pytest.raises(SystemExit) as exc:
        va.main()
    assert exc.value.code == 1
    assert "Provide --input" in capsys.readouterr().err


def test_main_exits_on_pipeline_error(monkeypatch, capsys):
    monkeypatch.setattr(
        va,
        "parse_args",
        lambda: argparse.Namespace(input=str(DEMO_VCF), demo=False, batch_size=10),
    )
    monkeypatch.setattr(
        va, "run_pipeline", lambda args: (_ for _ in ()).throw(RuntimeError("boom"))
    )

    with pytest.raises(SystemExit) as exc:
        va.main()
    assert exc.value.code == 1
    assert "ERROR: boom" in capsys.readouterr().err
