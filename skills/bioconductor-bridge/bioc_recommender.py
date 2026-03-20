#!/usr/bin/env python3
"""
bioc_recommender.py — Live Bioconductor package recommendation
==============================================================
Fetches package metadata from Bioconductor at runtime using
BiocManager-resolved repositories and the official VIEWS indexes.
"""

from __future__ import annotations

import csv
import io
import json
import os
import re
import shutil
import subprocess
import sys
from html.parser import HTMLParser
from pathlib import Path
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.parse import urljoin
from urllib.request import Request, urlopen

SKILL_DIR = Path(__file__).resolve().parent

DOMAIN_DEFINITIONS: dict[str, dict[str, Any]] = {
    "bulk_rnaseq": {
        "display_name": "Bulk RNA-seq Differential Expression",
        "aliases": [
            "bulk rna", "bulk rnaseq", "bulk rna-seq", "rna-seq", "rnaseq",
            "differential expression", "count matrix", "deseq2", "edger",
            "limma", "voom", "gene expression",
        ],
        "input_formats": [".csv", ".tsv", ".counts_matrix"],
        "containers": ["SummarizedExperiment"],
        "modalities": ["bulk-rnaseq"],
        "representative_packages": ["SummarizedExperiment", "DESeq2", "edgeR", "limma"],
    },
    "single_cell": {
        "display_name": "Single-cell RNA-seq",
        "aliases": [
            "single-cell", "single cell", "scrna", "10x", "matrix market",
            "singlecellexperiment", "scater", "scran", "batch correction",
            "integration",
        ],
        "input_formats": [".mtx", ".mtx.gz", ".h5ad"],
        "containers": ["SingleCellExperiment"],
        "modalities": ["single-cell"],
        "representative_packages": ["SingleCellExperiment", "scater", "scran", "scuttle", "batchelor"],
    },
    "genomic_ranges": {
        "display_name": "Genomic Ranges and Interval Operations",
        "aliases": [
            "granges", "genomicranges", "genomic ranges", "intervals",
            "genomic intervals", "bed", "gtf", "gff", "find overlaps",
            "peaks", "coordinates",
        ],
        "input_formats": [".bed", ".gtf", ".gff", ".gff3", ".bigwig", ".bw"],
        "containers": ["GRanges"],
        "modalities": ["genomics"],
        "representative_packages": ["GenomicRanges", "IRanges", "rtracklayer", "GenomicFeatures"],
    },
    "variant_annotation": {
        "display_name": "Variant Annotation and VCF Handling",
        "aliases": [
            "variantannotation", "annotate variants", "vcf", "variants",
            "variant effect", "functional annotation", "snp", "indel",
            "consequence",
        ],
        "input_formats": [".vcf", ".vcf.gz", ".bcf", ".bam", ".cram"],
        "containers": ["VCF"],
        "modalities": ["variants"],
        "representative_packages": ["VariantAnnotation", "Rsamtools", "BSgenome"],
    },
    "enrichment": {
        "display_name": "Pathway and Enrichment Analysis",
        "aliases": [
            "enrichment", "pathway", "gsea", "gene ontology", "go analysis",
            "kegg", "reactome", "pathway analysis",
        ],
        "input_formats": [".csv", ".tsv"],
        "containers": [".gene_list"],
        "modalities": ["enrichment"],
        "representative_packages": ["clusterProfiler", "fgsea", "AnnotationDbi"],
    },
    "methylation": {
        "display_name": "Methylation and Epigenomics",
        "aliases": [
            "methylation", "epigenomics", "epigenetics", "450k", "epic",
            "wgbs", "bisulfite", "beta values", "idat",
        ],
        "input_formats": [".idat", ".cov", ".bam"],
        "containers": ["RGChannelSet", "BSseq"],
        "modalities": ["methylation"],
        "representative_packages": ["minfi", "bsseq"],
    },
    "resource_hubs": {
        "display_name": "Annotation and Resource Hubs",
        "aliases": [
            "annotationhub", "experimenthub", "biocfilecache", "biomart",
            "resource hub", "annotation resources", "hub", "reference data",
        ],
        "input_formats": [],
        "containers": ["AnnotationHub", "ExperimentHub"],
        "modalities": ["annotation"],
        "representative_packages": ["AnnotationHub", "ExperimentHub", "BiocFileCache", "biomaRt"],
    },
    "visualization": {
        "display_name": "Visualization and Reporting",
        "aliases": [
            "heatmap", "visualization", "visualisation", "genome tracks",
            "gviz", "ggbio", "complexheatmap", "plot", "reporting",
        ],
        "input_formats": [".csv", ".tsv", ".bed", ".gtf", ".gff", ".vcf"],
        "containers": ["SummarizedExperiment", "SingleCellExperiment", "GRanges"],
        "modalities": ["visualization"],
        "representative_packages": ["ComplexHeatmap", "Gviz", "ggbio"],
    },
}

FORMAT_HINTS: dict[str, dict[str, str]] = {
    ".vcf": {"domain": "variant_annotation", "modality": "variants", "container": "VCF"},
    ".vcf.gz": {"domain": "variant_annotation", "modality": "variants", "container": "VCF"},
    ".bcf": {"domain": "variant_annotation", "modality": "variants", "container": "VCF"},
    ".bam": {"domain": "variant_annotation", "modality": "variants", "container": "VCF"},
    ".cram": {"domain": "variant_annotation", "modality": "variants", "container": "VCF"},
    ".mtx": {"domain": "single_cell", "modality": "single-cell", "container": "SingleCellExperiment"},
    ".mtx.gz": {"domain": "single_cell", "modality": "single-cell", "container": "SingleCellExperiment"},
    ".h5ad": {"domain": "single_cell", "modality": "single-cell", "container": "SingleCellExperiment"},
    ".bed": {"domain": "genomic_ranges", "modality": "genomics", "container": "GRanges"},
    ".gtf": {"domain": "genomic_ranges", "modality": "genomics", "container": "GRanges"},
    ".gff": {"domain": "genomic_ranges", "modality": "genomics", "container": "GRanges"},
    ".gff3": {"domain": "genomic_ranges", "modality": "genomics", "container": "GRanges"},
    ".bw": {"domain": "genomic_ranges", "modality": "genomics", "container": "GRanges"},
    ".bigwig": {"domain": "genomic_ranges", "modality": "genomics", "container": "GRanges"},
    ".idat": {"domain": "methylation", "modality": "methylation", "container": "RGChannelSet"},
    ".cov": {"domain": "methylation", "modality": "methylation", "container": "BSseq"},
    ".counts_matrix": {"domain": "bulk_rnaseq", "modality": "bulk-rnaseq", "container": "SummarizedExperiment"},
}

WORKFLOW_TEMPLATES: dict[str, dict[str, Any]] = {
    "bulk_rnaseq": {
        "name": "Bulk RNA-seq Differential Expression Workflow",
        "description": "Build a SummarizedExperiment from counts and metadata, then fit a count-aware differential expression model.",
        "container": "SummarizedExperiment",
        "steps": [
            {"package": "SummarizedExperiment", "role": "container", "description": "Construct a container from counts, row metadata, and sample metadata."},
            {"package": "DESeq2", "role": "primary", "description": "Run a default count-based differential expression workflow."},
            {"package": "edgeR", "role": "alternative", "description": "Use quasi-likelihood modeling when edgeR idioms are preferred."},
            {"package": "limma", "role": "alternative", "description": "Use voom + limma for linear-model-oriented workflows."},
            {"package": "ComplexHeatmap", "role": "visualization", "description": "Create reproducible expression heatmaps and sample annotations."},
            {"package": "clusterProfiler", "role": "follow_up", "description": "Perform pathway enrichment from ranked differential expression results."},
        ],
        "starter_template": "bulk_rnaseq",
    },
    "single_cell": {
        "name": "Single-cell RNA-seq Workflow",
        "description": "Organize assays in SingleCellExperiment, perform QC and normalization, and optionally integrate batches.",
        "container": "SingleCellExperiment",
        "steps": [
            {"package": "SingleCellExperiment", "role": "container", "description": "Create the canonical single-cell container for counts, metadata, and reduced dimensions."},
            {"package": "scater", "role": "primary", "description": "Run quality control, exploratory plots, and dimensionality reduction."},
            {"package": "scran", "role": "primary", "description": "Estimate normalization factors and support downstream clustering analyses."},
            {"package": "scuttle", "role": "supporting", "description": "Apply helper preprocessing operations and aggregation utilities."},
            {"package": "batchelor", "role": "optional", "description": "Apply batch-correction or integration methods when needed."},
        ],
        "starter_template": "single_cell",
    },
    "genomic_ranges": {
        "name": "Genomic Interval Workflow",
        "description": "Import annotation tracks, convert to GRanges, then perform overlap and feature annotation operations.",
        "container": "GRanges",
        "steps": [
            {"package": "rtracklayer", "role": "primary", "description": "Import BED, GTF, or BigWig files into Bioconductor objects."},
            {"package": "GenomicRanges", "role": "primary", "description": "Represent loci as GRanges and run overlap / nearest-feature operations."},
            {"package": "IRanges", "role": "supporting", "description": "Use interval infrastructure for coverage and range arithmetic."},
            {"package": "GenomicFeatures", "role": "supporting", "description": "Annotate overlaps against transcripts and gene models."},
            {"package": "Gviz", "role": "visualization", "description": "Create coordinate-aware genomic track visualizations."},
        ],
        "starter_template": "genomic_ranges",
    },
    "variant_annotation": {
        "name": "Variant Annotation Workflow",
        "description": "Load VCF data with VariantAnnotation, add reference context, and access annotation resources as needed.",
        "container": "VCF",
        "steps": [
            {"package": "VariantAnnotation", "role": "primary", "description": "Read VCF, filter variants, and compute annotation-aware summaries."},
            {"package": "Rsamtools", "role": "supporting", "description": "Use indexed genomic file access for BAM, Tabix, and VCF inputs."},
            {"package": "BSgenome", "role": "supporting", "description": "Provide reference genome context for sequence-aware annotation tasks."},
            {"package": "AnnotationHub", "role": "data_access", "description": "Retrieve curated annotation resources when needed."},
            {"package": "GenomicRanges", "role": "supporting", "description": "Manipulate annotated variant loci as genomic ranges."},
        ],
        "starter_template": "variant_annotation",
    },
    "enrichment": {
        "name": "Enrichment Analysis Workflow",
        "description": "Convert ranked or selected genes into pathway-level summaries with annotation support and visualization.",
        "container": ".gene_list",
        "steps": [
            {"package": "clusterProfiler", "role": "primary", "description": "Run over-representation or gene-set enrichment analyses."},
            {"package": "fgsea", "role": "primary", "description": "Use fast preranked GSEA when you already have ranked statistics."},
            {"package": "AnnotationDbi", "role": "supporting", "description": "Map identifiers and annotation databases for enrichment inputs."},
            {"package": "ComplexHeatmap", "role": "visualization", "description": "Visualize pathway-associated gene patterns when needed."},
        ],
        "starter_template": "enrichment",
    },
    "methylation": {
        "name": "Methylation and Epigenomics Workflow",
        "description": "Choose array-focused or sequencing-focused methylation tooling based on the input modality.",
        "container": "RGChannelSet",
        "steps": [
            {"package": "minfi", "role": "primary", "description": "Analyze Illumina methylation arrays and preprocess IDAT-based studies."},
            {"package": "bsseq", "role": "primary", "description": "Analyze bisulfite sequencing data and smoothed methylation profiles."},
            {"package": "GenomicRanges", "role": "supporting", "description": "Manipulate loci and methylation-associated genomic coordinates."},
        ],
        "starter_template": "methylation",
    },
    "resource_hubs": {
        "name": "Annotation and Resource Hub Workflow",
        "description": "Discover curated annotation or experiment resources and cache them locally for reproducible use.",
        "container": "AnnotationHub",
        "steps": [
            {"package": "AnnotationHub", "role": "primary", "description": "Search official annotation resources by organism, genome build, or keyword."},
            {"package": "ExperimentHub", "role": "primary", "description": "Locate experiment datasets and examples suitable for method development."},
            {"package": "BiocFileCache", "role": "supporting", "description": "Persist remote resources locally in a reproducible cache."},
            {"package": "biomaRt", "role": "supporting", "description": "Retrieve identifier mappings and gene annotations when BioMart access is appropriate."},
        ],
        "starter_template": "resource_hubs",
    },
    "visualization": {
        "name": "Bioconductor Visualization Workflow",
        "description": "Use matrix-aware or coordinate-aware visualization packages depending on the object model.",
        "container": "SummarizedExperiment",
        "steps": [
            {"package": "ComplexHeatmap", "role": "primary", "description": "Build high-quality heatmaps from assay matrices and sample annotations."},
            {"package": "Gviz", "role": "primary", "description": "Create track-based genome coordinate plots."},
            {"package": "ggbio", "role": "supporting", "description": "Use grammar-of-graphics-style genomics visualization for selected tasks."},
        ],
        "starter_template": "visualization",
    },
}

STARTER_TEMPLATES: dict[str, str] = {
    "bulk_rnaseq": """# Bulk RNA-seq starter workflow
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("SummarizedExperiment", "DESeq2", "ComplexHeatmap"), ask = FALSE, update = FALSE)

library(SummarizedExperiment)
library(DESeq2)
library(ComplexHeatmap)
""",
    "single_cell": """# Single-cell starter workflow
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("SingleCellExperiment", "scater", "scran"), ask = FALSE, update = FALSE)

library(SingleCellExperiment)
library(scater)
library(scran)
""",
    "genomic_ranges": """# Genomic ranges starter workflow
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GenomicRanges", "rtracklayer", "GenomicFeatures"), ask = FALSE, update = FALSE)

library(GenomicRanges)
library(rtracklayer)
library(GenomicFeatures)
""",
    "variant_annotation": """# Variant annotation starter workflow
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("VariantAnnotation", "AnnotationHub"), ask = FALSE, update = FALSE)

library(VariantAnnotation)
library(AnnotationHub)
""",
    "enrichment": """# Enrichment starter workflow
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "AnnotationDbi"), ask = FALSE, update = FALSE)

library(clusterProfiler)
""",
    "methylation": """# Methylation starter workflow
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("minfi", "bsseq"), ask = FALSE, update = FALSE)

library(minfi)
""",
    "resource_hubs": """# Annotation and resource hub starter workflow
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("AnnotationHub", "ExperimentHub", "BiocFileCache"), ask = FALSE, update = FALSE)

library(AnnotationHub)
library(ExperimentHub)
""",
    "visualization": """# Visualization starter workflow
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("ComplexHeatmap", "Gviz"), ask = FALSE, update = FALSE)

library(ComplexHeatmap)
library(Gviz)
""",
}

WORKFLOW_ROLE_MAP: dict[str, str] = {
    step["package"]: step["role"]
    for template in WORKFLOW_TEMPLATES.values()
    for step in template["steps"]
}

PACKAGE_DOMAIN_SIGNATURES: dict[str, list[str]] = {
    "bulk_rnaseq": ["rnaseq", "differential expression", "count data", "summarizedexperiment", "deseq2", "edger", "limma"],
    "single_cell": ["singlecell", "single-cell", "single cell", "singlecellexperiment", "scater", "scran", "scuttle", "batchelor"],
    "genomic_ranges": ["genomicranges", "granges", "iranges", "rtracklayer", "genomicfeatures", "interval"],
    "variant_annotation": ["variantannotation", "vcf", "rsamtools", "bsgenome"],
    "enrichment": ["clusterprofiler", "fgsea", "genesetenrichment", "pathway"],
    "methylation": ["methylation", "bisulfite", "wgbs", "minfi", "bsseq"],
    "resource_hubs": ["annotationhub", "experimenthub", "biocfilecache", "biomart", "hub"],
    "visualization": ["complexheatmap", "gviz", "ggbio", "heatmap", "visualization"],
}


def run_process(cmd: list[str], **kwargs) -> subprocess.CompletedProcess[str]:
    return subprocess.run(cmd, **kwargs)


def _normalize(text: str) -> str:
    return re.sub(r"[^a-z0-9.+_-]+", " ", (text or "").lower()).strip()


def _tokens(text: str) -> set[str]:
    return {tok for tok in _normalize(text).split() if tok}


GENERIC_QUERY_TOKENS = {
    "a",
    "an",
    "analysis",
    "analyze",
    "annotation",
    "annotations",
    "best",
    "data",
    "dataset",
    "datasets",
    "do",
    "for",
    "how",
    "i",
    "in",
    "is",
    "method",
    "methods",
    "of",
    "package",
    "packages",
    "recommend",
    "recommended",
    "show",
    "suggest",
    "tool",
    "tools",
    "use",
    "using",
    "what",
    "which",
    "with",
    "workflow",
    "workflows",
}

DOC_SECTION_STOP_LINES = {
    "documentation",
    "details",
    "package archives",
    "installation",
    "need some help? ask on the bioconductor support site!",
}

_DOC_CACHE: dict[tuple[str, str], dict[str, Any] | None] = {}


def _rscript_path() -> str | None:
    return shutil.which("Rscript")


def _run_rscript(expr: str, timeout: int = 600) -> subprocess.CompletedProcess[str]:
    rscript = _rscript_path()
    if not rscript:
        raise RuntimeError("Rscript was not found on PATH.")
    return run_process(
        [rscript, "--vanilla", "-e", expr],
        capture_output=True,
        text=True,
        timeout=timeout,
        cwd=str(SKILL_DIR),
    )


def _split_biocviews(value: str) -> list[str]:
    if not value:
        return []
    return [item.strip() for item in value.split(",") if item.strip()]


def _fixture_catalog() -> dict[str, Any] | None:
    fixture_path = os.environ.get("CLAWBIO_BIOC_FIXTURE_JSON")
    if not fixture_path:
        return None
    return json.loads(Path(fixture_path).read_text(encoding="utf-8"))


def _infer_domains_for_package(name: str, title: str, description: str, biocviews: list[str]) -> list[str]:
    searchable = " ".join([name, title, " ".join(biocviews)])
    query_norm = _normalize(searchable)
    domains: list[str] = []
    for domain_id, meta in DOMAIN_DEFINITIONS.items():
        if name in meta["representative_packages"]:
            domains.append(domain_id)
            continue
        if any(_normalize(sig) in query_norm for sig in PACKAGE_DOMAIN_SIGNATURES[domain_id]):
            domains.append(domain_id)
    return domains


def _infer_containers(name: str, title: str, description: str) -> list[str]:
    searchable = _normalize(" ".join([name, title, description]))
    containers = []
    for container in ["SummarizedExperiment", "SingleCellExperiment", "GRanges", "VCF", "AnnotationHub", "ExperimentHub"]:
        if _normalize(container) in searchable or name == container:
            containers.append(container)
    return containers


def _input_formats_for_domains(domains: list[str]) -> list[str]:
    formats: list[str] = []
    for domain in domains:
        for file_format in DOMAIN_DEFINITIONS[domain]["input_formats"]:
            if file_format not in formats:
                formats.append(file_format)
    return formats


def _modalities_for_domains(domains: list[str]) -> list[str]:
    modalities: list[str] = []
    for domain in domains:
        for modality in DOMAIN_DEFINITIONS[domain]["modalities"]:
            if modality not in modalities:
                modalities.append(modality)
    return modalities


def _curated_priority(name: str, domains: list[str]) -> int:
    priority = 5
    for domain in domains:
        if name in DOMAIN_DEFINITIONS[domain]["representative_packages"]:
            priority += 12
    if name in WORKFLOW_ROLE_MAP:
        priority += 8
    return priority


def _official_url(name: str) -> str:
    return f"https://bioconductor.org/packages/{name}/"


def load_catalog(timeout: int = 600) -> dict[str, Any]:
    """Fetch live package metadata from Bioconductor VIEWS indexes."""
    fixture = _fixture_catalog()
    if fixture is not None:
        return fixture

    expr = r'''
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  stop("BiocManager is required for live Bioconductor search. Install it with install.packages('BiocManager').")
}

repos <- BiocManager::repositories()
repos <- repos[grepl("bioconductor.org", repos, fixed = TRUE)]
repo_names <- names(repos)
records <- list()

for (i in seq_along(repos)) {
  repo <- repos[[i]]
  view_url <- paste0(sub("/+$", "", repo), "/VIEWS")
  view_data <- tryCatch(
    read.dcf(url(view_url), fields = c("Package", "Version", "Title", "Description", "biocViews")),
    error = function(e) NULL
  )
  if (is.null(view_data) || nrow(view_data) == 0) {
    next
  }
  frame <- as.data.frame(view_data, stringsAsFactors = FALSE)
  frame$Repository <- repo
  frame$RepoName <- if (length(repo_names) >= i) repo_names[[i]] else ""
  records[[length(records) + 1]] <- frame
}

if (length(records) == 0) {
  stop("No live Bioconductor package metadata could be retrieved from VIEWS indexes.")
}

catalog <- do.call(rbind, records)
write.table(catalog, file = "", sep = "\t", row.names = FALSE, quote = TRUE)
'''
    proc = _run_rscript(expr, timeout=timeout)
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr.strip() or proc.stdout.strip() or "Failed to retrieve live Bioconductor metadata.")

    reader = csv.DictReader(io.StringIO(proc.stdout), delimiter="\t")
    packages_by_name: dict[str, dict[str, Any]] = {}
    for row in reader:
        name = (row.get("Package") or "").strip()
        if not name:
            continue
        title = (row.get("Title") or "").strip()
        description = (row.get("Description") or "").strip()
        biocviews = _split_biocviews((row.get("biocViews") or "").strip())
        domains = _infer_domains_for_package(name, title, description, biocviews)
        package = {
            "name": name,
            "title": title,
            "description": description,
            "version": (row.get("Version") or "").strip(),
            "biocViews": biocviews,
            "repository": (row.get("Repository") or "").strip(),
            "repo_name": (row.get("RepoName") or "").strip(),
            "official_url": _official_url(name),
            "domains": domains,
            "containers": _infer_containers(name, title, description),
            "input_formats": _input_formats_for_domains(domains),
            "workflow_role": WORKFLOW_ROLE_MAP.get(name, "supporting"),
            "modalities": _modalities_for_domains(domains),
            "curated_priority": _curated_priority(name, domains),
            "aliases": [name.lower(), _normalize(name).replace(" ", "")],
        }
        if name not in packages_by_name or len(package["domains"]) > len(packages_by_name[name]["domains"]):
            packages_by_name[name] = package

    packages = sorted(packages_by_name.values(), key=lambda pkg: pkg["name"].lower())
    return {
        "version": "live",
        "source": "bioconductor_views",
        "package_count": len(packages),
        "packages": packages,
    }


def infer_input_context(input_path: str | Path | None) -> dict[str, str | None]:
    if not input_path:
        return {"input_format": None, "domain": None, "modality": None, "container": None}
    path = Path(input_path)
    if path.is_dir():
        if (path / "matrix.mtx").exists() or (path / "matrix.mtx.gz").exists():
            return {"input_format": ".mtx.gz" if (path / "matrix.mtx.gz").exists() else ".mtx", "domain": "single_cell", "modality": "single-cell", "container": "SingleCellExperiment"}
        return {"input_format": None, "domain": None, "modality": None, "container": None}
    suffixes = "".join(path.suffixes).lower()
    suffix = path.suffix.lower()
    fmt = suffixes if suffixes in FORMAT_HINTS else suffix
    if fmt in FORMAT_HINTS:
        return {"input_format": fmt, **FORMAT_HINTS[fmt]}
    if suffix in {".csv", ".tsv"} and path.exists():
        sep = "\t" if suffix == ".tsv" else ","
        with path.open("r", encoding="utf-8") as handle:
            reader = csv.reader(handle, delimiter=sep)
            header = next(reader, [])
            second = next(reader, [])
        headers = [_normalize(cell) for cell in header]
        if headers and headers[0] in {"gene", "gene_id", "symbol", "ensembl_id"} and len(headers) >= 4:
            numeric_count = 0
            for value in second[1:]:
                try:
                    float(value)
                    numeric_count += 1
                except ValueError:
                    continue
            if numeric_count >= 3:
                return {"input_format": ".counts_matrix", "domain": "bulk_rnaseq", "modality": "bulk-rnaseq", "container": "SummarizedExperiment"}
    return {"input_format": fmt if fmt else None, "domain": None, "modality": None, "container": None}


def detect_domain(query: str = "", input_format: str | None = None, modality: str | None = None, container: str | None = None) -> str | None:
    scores = {domain: 0.0 for domain in DOMAIN_DEFINITIONS}
    query_norm = _normalize(query)
    for domain, meta in DOMAIN_DEFINITIONS.items():
        for alias in meta["aliases"]:
            if _normalize(alias) and _normalize(alias) in query_norm:
                scores[domain] += 20
        if input_format and input_format in meta["input_formats"]:
            scores[domain] += 18
        if modality and modality in meta["modalities"]:
            scores[domain] += 14
        if container and container in meta["containers"]:
            scores[domain] += 16
    if input_format in FORMAT_HINTS:
        scores[FORMAT_HINTS[input_format]["domain"]] += 20
    best_domain, best_score = max(scores.items(), key=lambda item: item[1])
    return best_domain if best_score > 0 else None


def get_package_details(package_name: str, catalog: dict[str, Any]) -> dict[str, Any] | None:
    needle = _normalize(package_name)
    for package in catalog.get("packages", []):
        if _normalize(package["name"]) == needle:
            return package
        if needle in {_normalize(alias) for alias in package.get("aliases", [])}:
            return package
    return None


def _package_searchable_text(package: dict[str, Any]) -> str:
    return " ".join(
        [
            package.get("name", ""),
            package.get("title", ""),
            package.get("description", ""),
            " ".join(package.get("biocViews", [])),
            " ".join(package.get("domains", [])),
            " ".join(package.get("containers", [])),
        ]
    )


def _specific_query_tokens(text: str) -> set[str]:
    return {
        token
        for token in _tokens(text)
        if token not in GENERIC_QUERY_TOKENS and len(token) >= 3
    }


def _matched_tokens(query_tokens: set[str], searchable: str) -> set[str]:
    return {token for token in query_tokens if token in searchable}


class _PackagePageParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.lines: list[str] = []
        self._buffer: list[str] = []
        self._current_href: str | None = None
        self._current_link_text: list[str] = []
        self.links: list[dict[str, str]] = []
        self.h1_text: list[str] = []
        self.h2_text: list[str] = []
        self._heading_level: str | None = None
        self._skip_depth = 0

    def _flush(self) -> None:
        text = " ".join(part for part in self._buffer if part).strip()
        if text:
            self.lines.append(re.sub(r"\s+", " ", text))
        self._buffer = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        attrs_dict = dict(attrs)
        if tag in {"script", "style"}:
            self._skip_depth += 1
            return
        if tag in {"p", "div", "section", "article", "li", "ul", "ol", "br", "tr", "td", "table"}:
            self._flush()
        if tag in {"h1", "h2", "h3", "h4"}:
            self._flush()
            self._heading_level = tag
        if tag == "a":
            self._current_href = attrs_dict.get("href")
            self._current_link_text = []

    def handle_endtag(self, tag: str) -> None:
        if tag in {"script", "style"}:
            self._skip_depth = max(0, self._skip_depth - 1)
            return
        if self._skip_depth:
            return
        if tag in {"p", "div", "section", "article", "li", "ul", "ol", "br", "tr", "td", "table"}:
            self._flush()
        if tag in {"h1", "h2", "h3", "h4"}:
            self._flush()
            self._heading_level = None
        if tag == "a":
            text = " ".join(self._current_link_text).strip()
            if self._current_href and text:
                self.links.append({"text": re.sub(r"\s+", " ", text), "href": self._current_href})
            self._current_href = None
            self._current_link_text = []

    def handle_data(self, data: str) -> None:
        if self._skip_depth:
            return
        text = re.sub(r"\s+", " ", data).strip()
        if not text:
            return
        self._buffer.append(text)
        if self._current_href is not None:
            self._current_link_text.append(text)
        if self._heading_level == "h1":
            self.h1_text.append(text)
        elif self._heading_level == "h2":
            self.h2_text.append(text)

    def close(self) -> None:
        self._flush()
        super().close()


def _fetch_url_text(url: str, timeout: int = 30) -> str:
    request = Request(url, headers={"User-Agent": "ClawBio bioconductor-bridge/0.1"})
    with urlopen(request, timeout=timeout) as response:
        charset = response.headers.get_content_charset() or "utf-8"
        return response.read().decode(charset, errors="replace")


def _extract_doc_section(lines: list[str]) -> list[str]:
    start_index: int | None = None
    end_index: int | None = None
    for index, line in enumerate(lines):
        normalized = _normalize(line)
        if normalized == "documentation":
            start_index = index + 1
            continue
        if start_index is not None and normalized in {"details", "package archives"}:
            end_index = index
            break
    if start_index is None:
        return []
    section = [line for line in lines[start_index:end_index] if line.strip()]
    return section


def _doc_links_from_page(links: list[dict[str, str]], base_url: str) -> list[dict[str, str]]:
    doc_links: list[dict[str, str]] = []
    seen: set[str] = set()
    for link in links:
        href = link.get("href", "").strip()
        text = link.get("text", "").strip()
        if not href or not text:
            continue
        absolute_url = urljoin(base_url, href)
        normalized_text = _normalize(text)
        if (
            "/inst/doc/" in absolute_url
            or absolute_url.endswith(".pdf")
            or "support.bioconductor.org" in absolute_url
            or normalized_text in {"html", "pdf", "r script", "news", "text", "reference manual"}
        ):
            key = f"{text}|{absolute_url}"
            if key not in seen:
                seen.add(key)
                doc_links.append({"label": text, "url": absolute_url})
    return doc_links


def _summarize_documentation_lines(lines: list[str]) -> list[str]:
    summary: list[str] = []
    seen: set[str] = set()
    for line in lines:
        cleaned = re.sub(r"\s+", " ", line).strip()
        normalized = _normalize(cleaned)
        if not cleaned or normalized in DOC_SECTION_STOP_LINES:
            continue
        if normalized.startswith("to view documentation for the version of this package installed"):
            continue
        if normalized.startswith("browsevignettes"):
            continue
        if cleaned not in seen:
            seen.add(cleaned)
            summary.append(cleaned)
        if len(summary) >= 6:
            break
    return summary


def fetch_package_documentation(package: dict[str, Any], timeout: int = 30) -> dict[str, Any] | None:
    url = package.get("official_url")
    if not url:
        return None
    cache_key = (package.get("name", ""), package.get("version", ""))
    if cache_key in _DOC_CACHE:
        return _DOC_CACHE[cache_key]
    try:
        html = _fetch_url_text(url, timeout=timeout)
    except (HTTPError, URLError, TimeoutError, ValueError):
        _DOC_CACHE[cache_key] = None
        return None

    parser = _PackagePageParser()
    parser.feed(html)
    parser.close()
    lines = [line.strip() for line in parser.lines if line.strip()]
    documentation_lines = _extract_doc_section(lines)
    summary_lines = _summarize_documentation_lines(documentation_lines)
    combined_text = " ".join(
        [
            " ".join(parser.h1_text[:1]),
            " ".join(parser.h2_text[:1]),
            " ".join(lines[:80]),
        ]
    )
    profile = {
        "package": package.get("name", ""),
        "official_url": url,
        "page_title": " ".join(parser.h1_text[:1]).strip() or package.get("name", ""),
        "subtitle": " ".join(parser.h2_text[:1]).strip(),
        "documentation_lines": documentation_lines,
        "documentation_summary": summary_lines,
        "documentation_links": _doc_links_from_page(parser.links, url),
        "doc_text": _normalize(combined_text),
    }
    _DOC_CACHE[cache_key] = profile
    return profile


def _documentation_match(query: str, doc_profile: dict[str, Any] | None) -> tuple[float, list[str]]:
    if not doc_profile:
        return 0.0, []
    query_norm = _normalize(query)
    query_tokens = _specific_query_tokens(query)
    doc_text = doc_profile.get("doc_text", "")
    score = 0.0
    reasons: list[str] = []
    if query_norm and query_norm in doc_text:
        score += 20
        reasons.append("documentation contains the full query phrase")
    matched_tokens = _matched_tokens(query_tokens, doc_text)
    if matched_tokens:
        score += min(len(matched_tokens) * 10, 40)
        reasons.append(
            f"documentation matched query term{'s' if len(matched_tokens) != 1 else ''}: {', '.join(sorted(matched_tokens))}"
        )
    return score, reasons


def _merge_candidate_lists(*package_lists: list[dict[str, Any]]) -> list[dict[str, Any]]:
    merged: dict[str, dict[str, Any]] = {}
    for package_list in package_lists:
        for package in package_list:
            name = package.get("name")
            if not name:
                continue
            existing = merged.get(name)
            if existing is None or float(package.get("score", 0)) > float(existing.get("score", 0)):
                merged[name] = dict(package)
                continue
            combined = dict(existing)
            explanation_bits = [bit for bit in [existing.get("explanation", ""), package.get("explanation", "")] if bit]
            if explanation_bits:
                combined["explanation"] = "; ".join(dict.fromkeys(explanation_bits))
            for key in ["domains", "containers", "input_formats", "modalities", "biocViews", "aliases"]:
                combined_values = []
                for value in existing.get(key, []) + package.get(key, []):
                    if value not in combined_values:
                        combined_values.append(value)
                combined[key] = combined_values
            combined["score"] = max(float(existing.get("score", 0)), float(package.get("score", 0)))
            merged[name] = combined
    return sorted(merged.values(), key=lambda item: (-float(item.get("score", 0)), item.get("name", "")))


def rerank_packages_with_documentation(
    query: str,
    packages: list[dict[str, Any]],
    max_results: int = 5,
    max_fetch: int = 10,
    require_doc_match: bool = False,
) -> list[dict[str, Any]]:
    reranked: list[tuple[float, float, dict[str, Any]]] = []
    doc_match_found = False
    for index, package in enumerate(packages):
        enriched = dict(package)
        doc_score = 0.0
        doc_reasons: list[str] = []
        doc_profile: dict[str, Any] | None = None
        if index < max_fetch:
            doc_profile = fetch_package_documentation(package)
            doc_score, doc_reasons = _documentation_match(query, doc_profile)
            if doc_profile:
                enriched["documentation_summary"] = doc_profile.get("documentation_summary", [])
                enriched["documentation_links"] = doc_profile.get("documentation_links", [])
                enriched["documentation_subtitle"] = doc_profile.get("subtitle", "")
            if doc_score > 0:
                doc_match_found = True
        if require_doc_match and index < max_fetch and doc_score <= 0:
            continue
        base_score = float(enriched.get("score", 0))
        total_score = base_score + doc_score
        if doc_reasons:
            explanation_bits = [enriched.get("explanation", ""), *doc_reasons]
            enriched["explanation"] = "; ".join(bit for bit in dict.fromkeys(explanation_bits) if bit)
        enriched["documentation_score"] = round(doc_score, 2)
        enriched["score"] = round(total_score, 2)
        reranked.append((total_score, doc_score, enriched))

    if require_doc_match and doc_match_found:
        reranked = [item for item in reranked if item[1] > 0]

    reranked.sort(key=lambda item: (-item[0], -item[1], item[2]["name"]))
    return [item for _, _, item in reranked[:max_results]]


def search_catalog(query: str, catalog: dict[str, Any], max_results: int = 20) -> list[dict[str, Any]]:
    query_norm = _normalize(query)
    query_tokens = _specific_query_tokens(query)
    scored: list[tuple[float, dict[str, Any]]] = []
    for package in catalog.get("packages", []):
        searchable = _normalize(_package_searchable_text(package))
        score = 0.0
        reasons: list[str] = []
        if _normalize(package["name"]) == query_norm:
            score += 35
            reasons.append("exact package name match")
        if query_norm and query_norm in searchable:
            score += 18
            reasons.append("query phrase found in live Bioconductor metadata")
        matched_tokens = _matched_tokens(query_tokens, searchable)
        if matched_tokens:
            score += min(len(matched_tokens) * 7, 21)
            reasons.append(
                f"matched query term{'s' if len(matched_tokens) != 1 else ''}: {', '.join(sorted(matched_tokens))}"
            )
        if score > 0:
            enriched = dict(package)
            enriched["score"] = round(score, 2)
            enriched["explanation"] = "; ".join(reasons)
            scored.append((score, enriched))
    scored.sort(key=lambda item: (-item[0], item[1]["name"]))
    return [item for _, item in scored[:max_results]]


def recommend_packages(task: str, catalog: dict[str, Any], input_format: str | None = None, modality: str | None = None, container: str | None = None, max_results: int = 5) -> list[dict[str, Any]]:
    detected_domain = detect_domain(task, input_format=input_format, modality=modality, container=container)
    query_norm = _normalize(task)
    query_tokens = _specific_query_tokens(task)
    results: list[tuple[float, dict[str, Any]]] = []
    for package in catalog.get("packages", []):
        score = 0.0
        reasons: list[str] = []
        searchable = _normalize(_package_searchable_text(package))
        evidence_score = 0.0
        context_score = 0.0

        normalized_name = _normalize(package["name"])
        normalized_aliases = {_normalize(alias) for alias in package.get("aliases", [])}
        if query_norm and (normalized_name in query_norm or query_norm == normalized_name):
            evidence_score += 35
            reasons.append("package name mentioned directly")
        elif query_norm and query_norm in normalized_aliases:
            evidence_score += 33
            reasons.append("package alias mentioned directly")
        if query_norm and query_norm in searchable:
            evidence_score += 18
            reasons.append("query phrase found in live Bioconductor metadata")

        matched_tokens = _matched_tokens(query_tokens, searchable)
        if matched_tokens:
            evidence_score += min(len(matched_tokens) * 8, 32)
            reasons.append(
                f"matched query term{'s' if len(matched_tokens) != 1 else ''}: {', '.join(sorted(matched_tokens))}"
            )

        if detected_domain and detected_domain in package.get("domains", []):
            context_score += 10
            reasons.append(f"matched {DOMAIN_DEFINITIONS[detected_domain]['display_name']}")
        if container and container in package.get("containers", []):
            context_score += 8
            reasons.append(f"fits the {container} object model")
        if modality and modality in package.get("modalities", []):
            context_score += 6
            reasons.append(f"matches {modality} modality")
        if input_format and input_format in package.get("input_formats", []):
            context_score += 5
            reasons.append(f"supports {input_format} inputs")
        bioc_hits = sum(1 for view in package.get("biocViews", []) if _normalize(view) in query_norm)
        if bioc_hits:
            evidence_score += bioc_hits * 6
            reasons.append("matched biocViews terms")

        has_context = bool(context_score)
        has_evidence = bool(evidence_score)
        if not has_evidence and not has_context:
            continue
        if detected_domain and not has_evidence and detected_domain not in package.get("domains", []):
            continue

        score = evidence_score + context_score
        if has_evidence:
            score += min(float(package.get("curated_priority", 0)) / 10.0, 3.0)
            score += {
                "primary": 2.0,
                "container": 1.5,
                "alternative": 1.0,
                "supporting": 0.5,
                "data_access": 0.5,
                "visualization": 0.5,
                "follow_up": 0.25,
                "optional": 0.25,
            }.get(package.get("workflow_role", ""), 0.0)

        if score > 0:
            enriched = dict(package)
            enriched["score"] = round(score, 2)
            enriched["detected_domain"] = detected_domain
            enriched["explanation"] = "; ".join(dict.fromkeys(reasons)) or "live Bioconductor metadata match"
            results.append((score, enriched))
    results.sort(key=lambda item: (-item[0], item[1]["name"]))
    return [item for _, item in results[:max_results]]


def search_catalog_with_docs(query: str, catalog: dict[str, Any], max_results: int = 20, max_fetch: int = 10) -> list[dict[str, Any]]:
    candidates = search_catalog(query, catalog, max_results=max(max_results * 3, max_fetch))
    return rerank_packages_with_documentation(
        query=query,
        packages=candidates,
        max_results=max_results,
        max_fetch=max_fetch,
    )


def docs_search_catalog(query: str, catalog: dict[str, Any], max_results: int = 10, max_fetch: int = 12) -> list[dict[str, Any]]:
    candidates = _merge_candidate_lists(
        search_catalog(query, catalog, max_results=max(max_results * 4, max_fetch)),
        recommend_packages(query, catalog, max_results=max(max_results * 3, max_fetch)),
    )
    enriched = rerank_packages_with_documentation(
        query=query,
        packages=candidates,
        max_results=max_results,
        max_fetch=max_fetch,
        require_doc_match=True,
    )
    if enriched:
        return enriched
    return rerank_packages_with_documentation(
        query=query,
        packages=candidates,
        max_results=max_results,
        max_fetch=max_fetch,
    )


def recommend_packages_with_docs(
    task: str,
    catalog: dict[str, Any],
    input_format: str | None = None,
    modality: str | None = None,
    container: str | None = None,
    max_results: int = 5,
    max_fetch: int = 10,
) -> list[dict[str, Any]]:
    candidates = _merge_candidate_lists(
        recommend_packages(
            task=task,
            catalog=catalog,
            input_format=input_format,
            modality=modality,
            container=container,
            max_results=max(max_results * 4, max_fetch),
        ),
        search_catalog(task, catalog, max_results=max(max_results * 4, max_fetch)),
    )
    return rerank_packages_with_documentation(
        query=task,
        packages=candidates,
        max_results=max_results,
        max_fetch=max_fetch,
    )


def suggest_workflow(task: str, catalog: dict[str, Any], input_format: str | None = None, modality: str | None = None, container: str | None = None) -> dict[str, Any] | None:
    domain = detect_domain(task, input_format=input_format, modality=modality, container=container)
    if not domain:
        return None
    template = WORKFLOW_TEMPLATES[domain]
    workflow = {
        "id": domain,
        "domain": domain,
        "display_name": DOMAIN_DEFINITIONS[domain]["display_name"],
        "name": template["name"],
        "description": template["description"],
        "container": template["container"],
        "starter_template": template["starter_template"],
        "steps": [],
    }
    for step in template["steps"]:
        details = get_package_details(step["package"], catalog)
        step_info = dict(step)
        if details:
            step_info["official_url"] = details["official_url"]
            step_info["biocViews"] = details.get("biocViews", [])
            step_info["version"] = details.get("version", "")
        workflow["steps"].append(step_info)
    return workflow


def render_starter_script(workflow: dict[str, Any] | None) -> str:
    if not workflow:
        return """# Generic Bioconductor starter script
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
browseVignettes()
"""
    return STARTER_TEMPLATES.get(workflow.get("starter_template"), STARTER_TEMPLATES["resource_hubs"])


def list_domains() -> list[dict[str, Any]]:
    return [
        {
            "id": domain_id,
            "display_name": meta["display_name"],
            "representative_packages": meta["representative_packages"],
            "containers": meta["containers"],
            "input_formats": meta["input_formats"],
        }
        for domain_id, meta in DOMAIN_DEFINITIONS.items()
    ]
