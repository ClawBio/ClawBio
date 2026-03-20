from __future__ import annotations

from pathlib import Path
import sys

SKILL_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SKILL_DIR))

import bioc_recommender as rec  # noqa: E402


CATALOG = {
    "version": "live",
    "source": "test-fixture",
    "package_count": 10,
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
            "name": "edgeR",
            "title": "Empirical differential expression for count data",
            "description": "Bulk RNA-seq count-based differential expression.",
            "version": "1.0.0",
            "biocViews": ["RNASeq", "DifferentialExpression"],
            "repository": "https://bioconductor.org/packages/3.23/bioc/src/contrib",
            "repo_name": "BioCsoft",
            "official_url": "https://bioconductor.org/packages/edgeR/",
            "domains": ["bulk_rnaseq"],
            "containers": ["SummarizedExperiment"],
            "input_formats": [".counts_matrix"],
            "workflow_role": "alternative",
            "modalities": ["bulk-rnaseq"],
            "curated_priority": 18,
            "aliases": ["edger"],
        },
        {
            "name": "SingleCellExperiment",
            "title": "Single-cell container extending SummarizedExperiment",
            "description": "Single-cell assay container with reduced dimensions.",
            "version": "1.0.0",
            "biocViews": ["SingleCell", "Infrastructure"],
            "repository": "https://bioconductor.org/packages/3.23/bioc/src/contrib",
            "repo_name": "BioCsoft",
            "official_url": "https://bioconductor.org/packages/SingleCellExperiment/",
            "domains": ["single_cell"],
            "containers": ["SingleCellExperiment"],
            "input_formats": [".mtx", ".h5ad"],
            "workflow_role": "container",
            "modalities": ["single-cell"],
            "curated_priority": 20,
            "aliases": ["singlecellexperiment"],
        },
        {
            "name": "scater",
            "title": "Single-cell QC and visualization toolkit",
            "description": "Quality control and exploratory analysis for single-cell RNA-seq.",
            "version": "1.0.0",
            "biocViews": ["SingleCell", "QualityControl", "Visualization"],
            "repository": "https://bioconductor.org/packages/3.23/bioc/src/contrib",
            "repo_name": "BioCsoft",
            "official_url": "https://bioconductor.org/packages/scater/",
            "domains": ["single_cell", "visualization"],
            "containers": ["SingleCellExperiment"],
            "input_formats": [".mtx", ".h5ad"],
            "workflow_role": "primary",
            "modalities": ["single-cell"],
            "curated_priority": 18,
            "aliases": ["scater"],
        },
        {
            "name": "scran",
            "title": "Single-cell normalization and analysis",
            "description": "Single-cell normalization, HVG detection, and clustering support.",
            "version": "1.0.0",
            "biocViews": ["SingleCell", "Normalization"],
            "repository": "https://bioconductor.org/packages/3.23/bioc/src/contrib",
            "repo_name": "BioCsoft",
            "official_url": "https://bioconductor.org/packages/scran/",
            "domains": ["single_cell"],
            "containers": ["SingleCellExperiment"],
            "input_formats": [".mtx", ".h5ad"],
            "workflow_role": "primary",
            "modalities": ["single-cell"],
            "curated_priority": 18,
            "aliases": ["scran"],
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
            "name": "ATACseqQC",
            "title": "Quality assessment for ATAC-seq data",
            "description": "ATAC-seq quality control and analysis utilities.",
            "version": "1.0.0",
            "biocViews": ["ATACSeq", "QualityControl", "Sequencing"],
            "repository": "https://bioconductor.org/packages/3.23/bioc/src/contrib",
            "repo_name": "BioCsoft",
            "official_url": "https://bioconductor.org/packages/ATACseqQC/",
            "domains": [],
            "containers": [],
            "input_formats": [".bam"],
            "workflow_role": "supporting",
            "modalities": [],
            "curated_priority": 5,
            "aliases": ["atacseqqc"],
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
            "name": "ComplexHeatmap",
            "title": "Make complex heatmaps",
            "description": "Visualization package for richly annotated heatmaps.",
            "version": "1.0.0",
            "biocViews": ["Visualization"],
            "repository": "https://bioconductor.org/packages/3.23/bioc/src/contrib",
            "repo_name": "BioCsoft",
            "official_url": "https://bioconductor.org/packages/ComplexHeatmap/",
            "domains": ["visualization", "bulk_rnaseq"],
            "containers": ["SummarizedExperiment"],
            "input_formats": [".csv", ".tsv"],
            "workflow_role": "visualization",
            "modalities": ["visualization"],
            "curated_priority": 15,
            "aliases": ["complexheatmap"],
        },
    ],
}


def test_search_finds_annotationhub() -> None:
    results = rec.search_catalog("annotation resources hub", CATALOG, max_results=5)
    names = [pkg["name"] for pkg in results]
    assert "AnnotationHub" in names


def test_recommend_bulk_rnaseq_prefers_deseq2_family() -> None:
    results = rec.recommend_packages(
        "bulk RNA-seq differential expression from counts",
        CATALOG,
        input_format=".counts_matrix",
        modality="bulk-rnaseq",
        container="SummarizedExperiment",
        max_results=5,
    )
    names = [pkg["name"] for pkg in results]
    assert "DESeq2" in names
    assert "SummarizedExperiment" in names
    assert "edgeR" in names


def test_recommend_single_cell_prefers_sce_stack() -> None:
    results = rec.recommend_packages(
        "single-cell QC and clustering",
        CATALOG,
        input_format=".mtx",
        modality="single-cell",
        container="SingleCellExperiment",
        max_results=5,
    )
    names = [pkg["name"] for pkg in results]
    assert "SingleCellExperiment" in names
    assert "scater" in names
    assert "scran" in names


def test_recommend_vcf_prefers_variantannotation() -> None:
    results = rec.recommend_packages(
        "annotate this VCF",
        CATALOG,
        input_format=".vcf",
        modality="variants",
        container="VCF",
        max_results=5,
    )
    assert results[0]["name"] == "VariantAnnotation"


def test_recommend_ranges_prefers_genomicranges() -> None:
    results = rec.recommend_packages(
        "how do I work with genomic intervals",
        CATALOG,
        input_format=".bed",
        modality="genomics",
        container="GRanges",
        max_results=5,
    )
    assert results[0]["name"] == "GenomicRanges"


def test_recommend_prefers_literal_atac_match_over_curated_bulk_packages() -> None:
    results = rec.recommend_packages("atac analysis", CATALOG, max_results=5)
    assert results[0]["name"] == "ATACseqQC"
    assert "DESeq2" not in [pkg["name"] for pkg in results[:3]]


def test_recommend_with_docs_uses_documentation_to_refine_ranking(monkeypatch) -> None:
    def fake_docs(package, timeout=30):
        if package["name"] == "ATACseqQC":
            return {
                "subtitle": "ATAC-seq documentation",
                "documentation_summary": ["ATAC-seq quality control", "Open chromatin accessibility analysis"],
                "documentation_links": [],
                "doc_text": "atac open chromatin accessibility analysis peak qc",
            }
        return {
            "subtitle": f"{package['name']} docs",
            "documentation_summary": [],
            "documentation_links": [],
            "doc_text": "bulk rna differential expression counts",
        }

    monkeypatch.setattr(rec, "fetch_package_documentation", fake_docs)
    results = rec.recommend_packages_with_docs("ATAC accessibility analysis", CATALOG, max_results=5, max_fetch=5)
    assert results[0]["name"] == "ATACseqQC"
    assert results[0]["documentation_score"] > 0


def test_docs_search_prefers_doc_matches(monkeypatch) -> None:
    def fake_docs(package, timeout=30):
        text = "motif discovery transcription factor binding analysis" if package["name"] == "ATACseqQC" else "generic package docs"
        return {
            "subtitle": f"{package['name']} docs",
            "documentation_summary": ["Motif discovery workflow"] if package["name"] == "ATACseqQC" else [],
            "documentation_links": [],
            "doc_text": text,
        }

    monkeypatch.setattr(rec, "fetch_package_documentation", fake_docs)
    results = rec.docs_search_catalog("ATAC motif discovery", CATALOG, max_results=3, max_fetch=5)
    assert results[0]["name"] == "ATACseqQC"
    assert results[0]["documentation_score"] > 0


def test_workflow_template_is_deterministic() -> None:
    workflow = rec.suggest_workflow("single-cell RNA-seq preprocessing", CATALOG, input_format=".mtx")
    assert workflow is not None
    assert workflow["id"] == "single_cell"
    assert [step["package"] for step in workflow["steps"][:3]] == [
        "SingleCellExperiment",
        "scater",
        "scran",
    ]


def test_render_starter_script_for_variant_workflow() -> None:
    workflow = rec.suggest_workflow("annotate variants from a VCF", CATALOG, input_format=".vcf")
    script = rec.render_starter_script(workflow)
    assert "VariantAnnotation" in script
    assert "AnnotationHub" in script


def test_infer_vcf_context_without_existing_file() -> None:
    ctx = rec.infer_input_context("sample.vcf.gz")
    assert ctx["domain"] == "variant_annotation"
    assert ctx["container"] == "VCF"


def test_infer_mtx_context_without_existing_file() -> None:
    ctx = rec.infer_input_context("matrix.mtx.gz")
    assert ctx["domain"] == "single_cell"
    assert ctx["container"] == "SingleCellExperiment"


def test_infer_granges_context_without_existing_file() -> None:
    ctx = rec.infer_input_context("regions.bed")
    assert ctx["domain"] == "genomic_ranges"
    assert ctx["container"] == "GRanges"


def test_infer_count_matrix_context_from_csv(tmp_path: Path) -> None:
    counts = tmp_path / "counts.csv"
    counts.write_text("gene,s1,s2,s3\nTP53,10,11,12\n", encoding="utf-8")
    ctx = rec.infer_input_context(counts)
    assert ctx["domain"] == "bulk_rnaseq"
    assert ctx["container"] == "SummarizedExperiment"
