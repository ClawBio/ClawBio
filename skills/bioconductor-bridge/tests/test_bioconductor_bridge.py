from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

SKILL_DIR = Path(__file__).resolve().parents[1]
PROJECT_ROOT = SKILL_DIR.parents[1]
sys.path.insert(0, str(SKILL_DIR))

import bioconductor_bridge as bridge  # noqa: E402


SAMPLE_CATALOG = {
    "version": "live",
    "source": "test-fixture",
    "package_count": 5,
    "packages": [
        {
            "name": "SummarizedExperiment",
            "title": "SummarizedExperiment container",
            "description": "Container for assay matrices and sample metadata.",
            "version": "1.0.0",
            "biocViews": ["Infrastructure", "Transcriptomics"],
            "repository": "https://bioconductor.org/packages/3.23/bioc/src/contrib",
            "repo_name": "BioCsoft",
            "official_url": "https://bioconductor.org/packages/SummarizedExperiment/",
            "domains": ["bulk_rnaseq"],
            "containers": ["SummarizedExperiment"],
            "input_formats": [".counts_matrix"],
            "workflow_role": "container",
            "modalities": ["bulk-rnaseq"],
            "curated_priority": 20,
            "aliases": ["summarizedexperiment"],
        },
        {
            "name": "DESeq2",
            "title": "Differential gene expression for RNA-seq",
            "description": "Bulk RNA-seq differential expression modeling from count data.",
            "version": "1.0.0",
            "biocViews": ["RNASeq", "DifferentialExpression"],
            "repository": "https://bioconductor.org/packages/3.23/bioc/src/contrib",
            "repo_name": "BioCsoft",
            "official_url": "https://bioconductor.org/packages/DESeq2/",
            "domains": ["bulk_rnaseq"],
            "containers": ["SummarizedExperiment"],
            "input_formats": [".counts_matrix"],
            "workflow_role": "primary",
            "modalities": ["bulk-rnaseq"],
            "curated_priority": 22,
            "aliases": ["deseq2"],
        },
        {
            "name": "VariantAnnotation",
            "title": "Annotation of genetic variants",
            "description": "Read and annotate VCF files with Bioconductor.",
            "version": "1.0.0",
            "biocViews": ["Annotation", "VariantAnnotation"],
            "repository": "https://bioconductor.org/packages/3.23/bioc/src/contrib",
            "repo_name": "BioCsoft",
            "official_url": "https://bioconductor.org/packages/VariantAnnotation/",
            "domains": ["variant_annotation"],
            "containers": ["VCF"],
            "input_formats": [".vcf", ".vcf.gz"],
            "workflow_role": "primary",
            "modalities": ["variants"],
            "curated_priority": 22,
            "aliases": ["variantannotation"],
        },
        {
            "name": "AnnotationHub",
            "title": "Client to access annotation resources",
            "description": "Search current annotation resources from Bioconductor.",
            "version": "1.0.0",
            "biocViews": ["Annotation", "ExperimentHub"],
            "repository": "https://bioconductor.org/packages/3.23/bioc/src/contrib",
            "repo_name": "BioCsoft",
            "official_url": "https://bioconductor.org/packages/AnnotationHub/",
            "domains": ["resource_hubs", "variant_annotation"],
            "containers": ["AnnotationHub"],
            "input_formats": [],
            "workflow_role": "data_access",
            "modalities": ["annotation"],
            "curated_priority": 18,
            "aliases": ["annotationhub"],
        },
        {
            "name": "GenomicRanges",
            "title": "Representation and manipulation of genomic intervals",
            "description": "Work with genomic ranges and overlaps using GRanges.",
            "version": "1.0.0",
            "biocViews": ["Genomics", "Annotation"],
            "repository": "https://bioconductor.org/packages/3.23/bioc/src/contrib",
            "repo_name": "BioCsoft",
            "official_url": "https://bioconductor.org/packages/GenomicRanges/",
            "domains": ["genomic_ranges", "variant_annotation"],
            "containers": ["GRanges"],
            "input_formats": [".bed", ".gtf"],
            "workflow_role": "primary",
            "modalities": ["genomics"],
            "curated_priority": 20,
            "aliases": ["genomicranges"],
        },
    ],
}


def write_fixture(tmp_path: Path) -> Path:
    fixture = tmp_path / "live_catalog_fixture.json"
    fixture.write_text(json.dumps(SAMPLE_CATALOG), encoding="utf-8")
    return fixture


def test_invalid_primary_mode_combo_raises() -> None:
    with pytest.raises(SystemExit):
        bridge.main(["--demo", "--recommend", "bulk RNA-seq"])


def test_install_requires_explicit_packages() -> None:
    with pytest.raises(SystemExit):
        bridge.validate_args(bridge.build_parser().parse_args(["--install", " , "]))


def test_setup_reports_live_validation_field(monkeypatch: pytest.MonkeyPatch) -> None:
    outputs = {
        'cat(R.version.string)': subprocess.CompletedProcess([], 0, stdout="R 4.5.0\n", stderr=""),
        'cat(requireNamespace("BiocManager", quietly = TRUE))': subprocess.CompletedProcess([], 0, stdout="TRUE\n", stderr=""),
    }

    def fake_run(expr: str, timeout: int = 120):
        if "packageVersion" in expr:
            return subprocess.CompletedProcess([], 0, stdout="biocmanager_version=1.30.27\nbioconductor_version=3.23\n", stderr="")
        if "BiocManager::valid()" in expr:
            return subprocess.CompletedProcess([], 0, stdout="remote_validation=ok\n", stderr="")
        if "requireNamespace(pkg" in expr:
            return subprocess.CompletedProcess([], 0, stdout="DESeq2=TRUE\n", stderr="")
        return outputs.get(expr, subprocess.CompletedProcess([], 0, stdout="", stderr=""))

    monkeypatch.setattr(bridge, "_rscript_path", lambda: "/usr/local/bin/Rscript")
    monkeypatch.setattr(bridge, "_run_rscript", fake_run)
    status = bridge.inspect_setup(["DESeq2"])
    assert status["remote_validation"] == "validated live with BiocManager::valid()"
    assert status["bioconductor_version"] == "3.23"


def test_setup_detects_bioc_versions_from_r_output(monkeypatch: pytest.MonkeyPatch) -> None:
    def fake_run(expr: str, timeout: int = 120):
        if expr == 'cat(R.version.string)':
            return subprocess.CompletedProcess([], 0, stdout="R 4.5.0\n", stderr="")
        if expr == 'cat(requireNamespace("BiocManager", quietly = TRUE))':
            return subprocess.CompletedProcess([], 0, stdout="TRUE\n", stderr="")
        if 'packageVersion("BiocManager")' in expr:
            return subprocess.CompletedProcess([], 0, stdout="biocmanager_version=1.30.27\nbioconductor_version=3.23\n", stderr="")
        if "BiocManager::valid()" in expr:
            return subprocess.CompletedProcess([], 0, stdout="remote_validation=ok\n", stderr="")
        if "requireNamespace(pkg" in expr:
            return subprocess.CompletedProcess([], 0, stdout="MotifPeeker=TRUE\n", stderr="")
        return subprocess.CompletedProcess([], 1, stdout="", stderr="unexpected expr")

    monkeypatch.setattr(bridge, "_rscript_path", lambda: "/usr/local/bin/Rscript")
    monkeypatch.setattr(bridge, "_run_rscript", fake_run)
    status = bridge.inspect_setup(["MotifPeeker"])
    assert status["biocmanager_version"] == "1.30.27"
    assert status["bioconductor_version"] == "3.23"


def test_install_command_construction() -> None:
    expr = bridge.build_install_expr(["DESeq2", "ComplexHeatmap"])
    assert 'BiocManager::install(c("DESeq2", "ComplexHeatmap"), ask = FALSE, update = FALSE)' in expr


def test_demo_creates_report_and_result_json(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(bridge, "load_catalog", lambda: SAMPLE_CATALOG)
    monkeypatch.setattr(
        bridge,
        "inspect_setup",
        lambda packages: {
            "rscript_path": "/usr/local/bin/Rscript",
            "r_available": True,
            "r_version": "R 4.5.0",
            "r_status": "release",
            "warnings": [],
            "biocmanager_installed": True,
            "biocmanager_version": "1.30.27",
            "bioconductor_version": "3.23",
            "installed_packages": {pkg: False for pkg in packages},
            "remote_validation": "validated live with BiocManager::valid()",
        },
    )
    monkeypatch.setattr(bridge, "capture_session_info", lambda: "sessionInfo()\n")

    exit_code = bridge.main(["--demo", "--output", str(tmp_path)])
    assert exit_code == 0
    result = json.loads((tmp_path / "result.json").read_text(encoding="utf-8"))
    assert result["summary"]["mode"] == "demo"
    assert result["summary"]["detected_domain"] == "bulk_rnaseq"


def test_search_handles_package_without_inferred_domain(monkeypatch: pytest.MonkeyPatch) -> None:
    args = bridge.build_parser().parse_args(["--search", "MotifPeeker"])
    monkeypatch.setattr(
        bridge,
        "search_catalog_with_docs",
        lambda query, catalog, max_results=5: [
            {
                "name": "MotifPeeker",
                "title": "Motif discovery package",
                "description": "Live package with no inferred workflow domain.",
                "domains": [],
                "biocViews": ["MotifDiscovery"],
                "official_url": "https://bioconductor.org/packages/MotifPeeker/",
                "workflow_role": "supporting",
                "containers": [],
                "input_formats": [],
                "modalities": [],
                "score": 42,
                "explanation": "exact package name match",
            }
        ],
    )
    monkeypatch.setattr(bridge, "suggest_workflow", lambda *args, **kwargs: None)

    payload = bridge.execute_mode(args, {"packages": []}, {"input_format": None, "modality": None, "container": None, "input_path": None})
    assert payload["recommendations"][0]["name"] == "MotifPeeker"
    assert payload["detected_domain"] is None


def test_package_docs_mode_includes_documentation(monkeypatch: pytest.MonkeyPatch) -> None:
    args = bridge.build_parser().parse_args(["--package-docs", "DESeq2"])
    monkeypatch.setattr(bridge, "inspect_setup", lambda packages: {"r_available": False, "r_version": "", "r_status": "missing", "biocmanager_installed": False, "biocmanager_version": "", "bioconductor_version": "", "installed_packages": {}, "remote_validation": "not attempted", "warnings": []})
    monkeypatch.setattr(
        bridge,
        "fetch_package_documentation",
        lambda package: {
            "official_url": package["official_url"],
            "subtitle": "Differential gene expression analysis based on the negative binomial distribution",
            "documentation_summary": ["Analyzing RNA-seq data with DESeq2 HTML R Script", "Reference Manual PDF"],
            "documentation_links": [{"label": "HTML", "url": "https://bioconductor.org/packages/DESeq2/inst/doc/DESeq2.html"}],
        },
    )

    payload = bridge.execute_mode(args, SAMPLE_CATALOG, {"input_format": None, "modality": None, "container": None, "input_path": None})
    assert payload["documentation"]["subtitle"].startswith("Differential gene expression")
    assert payload["recommendations"][0]["documentation_summary"][0].startswith("Analyzing RNA-seq data")


def test_docs_search_mode_uses_doc_search(monkeypatch: pytest.MonkeyPatch) -> None:
    args = bridge.build_parser().parse_args(["--docs-search", "ATAC analysis"])
    monkeypatch.setattr(
        bridge,
        "docs_search_catalog",
        lambda query, catalog, max_results=5: [
            {
                "name": "ATACseqQC",
                "title": "Quality assessment for ATAC-seq data",
                "description": "ATAC docs match",
                "domains": [],
                "biocViews": ["ATACSeq"],
                "official_url": "https://bioconductor.org/packages/ATACseqQC/",
                "workflow_role": "supporting",
                "containers": [],
                "input_formats": [".bam"],
                "modalities": [],
                "score": 44,
                "documentation_score": 20,
                "documentation_summary": ["ATAC-seq quality control", "Open chromatin accessibility analysis"],
                "documentation_links": [],
                "explanation": "documentation matched query terms: atac, analysis",
            }
        ],
    )
    monkeypatch.setattr(bridge, "inspect_setup", lambda packages: {"r_available": False, "r_version": "", "r_status": "missing", "biocmanager_installed": False, "biocmanager_version": "", "bioconductor_version": "", "installed_packages": {}, "remote_validation": "not attempted", "warnings": []})
    monkeypatch.setattr(bridge, "suggest_workflow", lambda *args, **kwargs: None)

    payload = bridge.execute_mode(args, SAMPLE_CATALOG, {"input_format": None, "modality": None, "container": None, "input_path": None})
    assert payload["mode"] == "docs_search"
    assert payload["recommendations"][0]["name"] == "ATACseqQC"


def test_clawbio_run_bioc_demo_smoke(tmp_path: Path) -> None:
    fixture = write_fixture(tmp_path)
    env = os.environ.copy()
    env["CLAWBIO_BIOC_FIXTURE_JSON"] = str(fixture)
    proc = subprocess.run(
        [
            sys.executable,
            str(PROJECT_ROOT / "clawbio.py"),
            "run",
            "bioc",
            "--demo",
            "--output",
            str(tmp_path),
        ],
        cwd=str(PROJECT_ROOT),
        capture_output=True,
        text=True,
        timeout=180,
        env=env,
    )
    assert proc.returncode == 0, proc.stderr
    assert (tmp_path / "report.md").exists()
    assert (tmp_path / "result.json").exists()


def test_clawbio_run_bioc_recommend_smoke(tmp_path: Path) -> None:
    fixture = write_fixture(tmp_path)
    env = os.environ.copy()
    env["CLAWBIO_BIOC_FIXTURE_JSON"] = str(fixture)
    proc = subprocess.run(
        [
            sys.executable,
            str(PROJECT_ROOT / "clawbio.py"),
            "run",
            "bioc",
            "--recommend",
            "annotate this VCF",
            "--format",
            ".vcf",
            "--output",
            str(tmp_path),
        ],
        cwd=str(PROJECT_ROOT),
        capture_output=True,
        text=True,
        timeout=180,
        env=env,
    )
    assert proc.returncode == 0, proc.stderr
    result = json.loads((tmp_path / "result.json").read_text(encoding="utf-8"))
    names = [pkg["name"] for pkg in result["data"]["recommendations"]]
    assert "VariantAnnotation" in names


def test_bio_orchestrator_routes_bioc_queries() -> None:
    source = (PROJECT_ROOT / "skills" / "bio-orchestrator" / "orchestrator.py").read_text(encoding="utf-8")
    assert '"bioconductor": "bioconductor-bridge"' in source
    assert '"which bioconductor package": "bioconductor-bridge"' in source
    assert '"bioconductor-bridge": "bioc"' in source
