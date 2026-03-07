#!/usr/bin/env python3
"""
tool_recommender.py — Intelligent Galaxy tool selection for ClawBio
====================================================================
Given a task description (and optional input format), recommends the best
Galaxy tool(s) and suggests multi-step pipelines.  Uses pre-built lookup
tables — no ML, no API calls, no external dependencies.

Key features:
  - Version deduplication: 8,182 tools → ~2,300 unique
  - EDAM ontology resolution: raw IDs → human-readable labels
  - 15-category task taxonomy with natural-language aliases
  - Multi-signal scoring across 7 dimensions
  - 6 pre-defined workflow templates for common pipelines
  - Input format awareness for file-type-driven recommendations
"""

from __future__ import annotations

import json
import math
import re
import sys
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parent

# ---------------------------------------------------------------------------
# EDAM ontology lookup (108 entries — all IDs found in the Galaxy catalog)
# ---------------------------------------------------------------------------

EDAM_TOPICS: dict[str, str] = {
    "topic_0080": "Sequence analysis",
    "topic_0082": "Structure analysis",
    "topic_0085": "Functional genomics",
    "topic_0091": "Bioinformatics",
    "topic_0092": "Data visualisation",
    "topic_0102": "Mapping",
    "topic_0121": "Pharmacology",
    "topic_0157": "Sequence composition, complexity and repeats",
    "topic_0194": "Structure determination",
    "topic_0196": "Sequence assembly",
    "topic_0199": "Genetic variation",
    "topic_0203": "Gene expression",
    "topic_0610": "Cell biology",
    "topic_0622": "Genomics",
    "topic_0632": "Taxonomy",
    "topic_0797": "Comparative genomics",
    "topic_0798": "Mobile genetic elements",
    "topic_2269": "Statistics and probability",
    "topic_3047": "Molecular interactions, pathways and networks",
    "topic_3050": "Biodiversity",
    "topic_3168": "Sequencing",
    "topic_3170": "RNA-seq",
    "topic_3172": "DNA methylation",
    "topic_3174": "Epigenetics",
    "topic_3295": "Phylogenetics",
    "topic_3301": "Microbiology",
    "topic_3307": "Computational biology",
    "topic_3308": "Transcriptomics",
    "topic_3318": "Phylogenomics",
    "topic_3320": "DNA binding sites",
    "topic_3332": "Computational chemistry",
    "topic_3383": "Biological imaging",
    "topic_3511": "Nucleic acid sequence analysis",
    "topic_3512": "Gene transcripts",
    "topic_3520": "Proteomics",
    "topic_3570": "Ecology",
    "topic_3673": "Whole genome sequencing",
    "topic_3697": "Metagenomics",
    "topic_3855": "Environmental science",
    "topic_4028": "ChIP-seq",
}

EDAM_OPERATIONS: dict[str, str] = {
    "operation_0224": "Query and retrieval",
    "operation_0231": "Sequence editing",
    "operation_0233": "Sequence conversion",
    "operation_0237": "Repeat sequence analysis",
    "operation_0239": "Sequence motif recognition",
    "operation_0262": "Genome assembly",
    "operation_0284": "Phylogenetic tree generation",
    "operation_0286": "Profile-profile alignment",
    "operation_0292": "Sequence alignment",
    "operation_0324": "Phylogenetic tree analysis",
    "operation_0335": "Formatting",
    "operation_0336": "Format validation",
    "operation_0337": "Visualisation",
    "operation_0362": "Genome annotation",
    "operation_0368": "Sequence masking",
    "operation_0369": "Sequence cutting",
    "operation_0436": "Coding region prediction",
    "operation_0474": "Protein structure prediction",
    "operation_0491": "Pairwise sequence alignment",
    "operation_0496": "Global alignment",
    "operation_0525": "Genome assembly",
    "operation_0564": "Sequence visualisation",
    "operation_0567": "Phylogenetic tree visualisation",
    "operation_0573": "Map drawing",
    "operation_1777": "Protein function prediction",
    "operation_1812": "Parsing",
    "operation_2238": "Statistical calculation",
    "operation_2422": "Data retrieval",
    "operation_2426": "Modelling and simulation",
    "operation_2428": "Validation",
    "operation_2436": "Gene-set enrichment analysis",
    "operation_2478": "Nucleic acid sequence analysis",
    "operation_2495": "Expression analysis",
    "operation_2940": "Classification",
    "operation_2945": "Analysis",
    "operation_3087": "Protein structure comparison",
    "operation_3096": "Editing",
    "operation_3180": "Sequence assembly validation",
    "operation_3187": "Sequence alignment refinement",
    "operation_3192": "Sequence trimming",
    "operation_3197": "Sequence composition calculation",
    "operation_3200": "Community profiling",
    "operation_3207": "Genome visualisation",
    "operation_3218": "Sequencing quality control",
    "operation_3223": "Differential gene expression profiling",
    "operation_3227": "Variant calling",
    "operation_3357": "Format detection",
    "operation_3359": "Splitting",
    "operation_3429": "Generation",
    "operation_3430": "Nucleic acid sequence composition analysis",
    "operation_3432": "Clustering",
    "operation_3434": "Conversion",
    "operation_3435": "Standardisation and normalisation",
    "operation_3436": "Aggregation",
    "operation_3443": "Image analysis",
    "operation_3460": "Taxonomic classification",
    "operation_3463": "Expression profile comparison",
    "operation_3465": "Correlation",
    "operation_3482": "Antimicrobial resistance prediction",
    "operation_3557": "Imputation",
    "operation_3563": "RNA-seq quantification",
    "operation_3627": "Mass spectra calibration",
    "operation_3629": "Deisotoping",
    "operation_3672": "Gene functional annotation",
    "operation_3695": "Filtering",
    "operation_3800": "RNA-seq read count analysis",
    "operation_3860": "Spectrum calculation",
    "operation_3946": "RNA-seq time series analysis",
}


# ---------------------------------------------------------------------------
# Task taxonomy — 15 bioinformatics task categories
# ---------------------------------------------------------------------------

TASK_TAXONOMY: dict[str, dict] = {
    "quality_control": {
        "display_name": "Read Quality Control",
        "aliases": [
            "qc", "quality control", "quality check", "read quality",
            "sequence quality", "check my reads", "fastq quality",
            "trim reads", "adapter removal", "quality filter",
            "q30", "phred", "base quality",
        ],
        "sections": ["FASTQ Quality Control", "FASTA/FASTQ"],
        "preferred_tools": ["fastqc", "fastp", "trimmomatic", "cutadapt", "multiqc"],
        "input_formats": [".fastq", ".fq", ".fastq.gz", ".fq.gz"],
        "description": "Assess and improve the quality of sequencing reads",
    },
    "read_mapping": {
        "display_name": "Read Mapping / Alignment",
        "aliases": [
            "mapping", "alignment", "align reads", "map reads",
            "align to reference", "read alignment", "map to genome",
            "short read alignment", "long read alignment",
            "bwa", "bowtie", "hisat",
        ],
        "sections": ["Mapping"],
        "preferred_tools": ["bwa", "bwa_mem2", "bowtie2", "hisat2", "minimap2", "star"],
        "input_formats": [".fastq", ".fq", ".fastq.gz", ".fq.gz", ".fasta", ".fa"],
        "description": "Align sequencing reads to a reference genome",
    },
    "variant_calling": {
        "display_name": "Variant Calling",
        "aliases": [
            "variant calling", "call variants", "snp calling",
            "find mutations", "detect variants", "genotyping",
            "somatic variants", "germline variants", "indels",
            "gatk", "haplotypecaller", "freebayes", "deepvariant",
            "snp", "indel", "variant discovery",
        ],
        "sections": ["Variant Calling", "VCF/BCF"],
        "preferred_tools": ["freebayes", "gatk4", "bcftools", "deepvariant", "lofreq"],
        "input_formats": [".bam", ".cram"],
        "description": "Identify genetic variants (SNPs, indels) from aligned reads",
    },
    "variant_annotation": {
        "display_name": "Variant Annotation & Filtering",
        "aliases": [
            "annotate variants", "variant annotation", "variant effect",
            "predict variant impact", "clinvar", "vep", "snpeff",
            "functional annotation", "pathogenic", "benign",
            "variant filtering", "filter variants", "annovar",
            "missense", "stopgain", "splicing", "frameshift",
            "dbsnp", "gnomad", "exac", "cadd",
        ],
        "sections": ["VCF/BCF", "Annotation"],
        "preferred_tools": ["snpsift", "snpeff", "vep", "gemini"],
        "input_formats": [".vcf", ".vcf.gz", ".bcf"],
        "description": "Annotate variants with functional impact, population frequency, and clinical significance",
    },
    "wes_wgs_analysis": {
        "display_name": "Whole Exome / Whole Genome Analysis",
        "aliases": [
            "exome", "wes", "whole exome", "wgs", "whole genome",
            "exome sequencing", "capture", "enrichment",
            "coverage", "depth", "on-target", "off-target",
            "germline analysis", "clinical exome",
            "novogene", "illumina", "novaseq",
        ],
        "sections": ["Variant Calling", "VCF/BCF", "SAM/BAM", "Picard"],
        "preferred_tools": ["gatk4", "freebayes", "deepvariant", "picard", "samtools", "bcftools", "snpsift", "snpeff"],
        "input_formats": [".bam", ".vcf", ".vcf.gz", ".fastq", ".fastq.gz"],
        "description": "End-to-end analysis of whole exome or whole genome sequencing data",
    },
    "rnaseq": {
        "display_name": "RNA-seq Analysis",
        "aliases": [
            "rna-seq", "rnaseq", "rna seq", "differential expression",
            "gene expression", "transcriptomics", "deseq", "edger",
            "count reads", "transcript quantification",
        ],
        "sections": ["RNA-seq", "RNA Analysis"],
        "preferred_tools": [
            "featurecounts", "htseq_count", "deseq2", "edger",
            "stringtie", "salmon", "kallisto",
        ],
        "input_formats": [".bam", ".fastq", ".fq", ".fastq.gz"],
        "description": "Quantify and compare gene expression from RNA-seq data",
    },
    "metagenomics": {
        "display_name": "Metagenomics / Microbiome",
        "aliases": [
            "metagenomics", "microbiome", "16s", "amplicon",
            "taxonomic classification", "species identification",
            "microbial community", "metagenomic profiling",
            "classify species", "identify organisms",
        ],
        "sections": ["Metagenomic Analysis", "QIIME2", "Mothur"],
        "preferred_tools": [
            "kraken2", "metaphlan", "humann", "bracken",
            "diamond", "megahit",
        ],
        "input_formats": [".fastq", ".fq", ".fastq.gz", ".fq.gz", ".fasta"],
        "description": "Profile microbial communities from metagenomic sequences",
    },
    "genome_assembly": {
        "display_name": "Genome Assembly",
        "aliases": [
            "assembly", "assemble", "de novo assembly",
            "genome assembly", "scaffold", "contig",
            "assemble reads", "build genome",
        ],
        "sections": ["Assembly"],
        "preferred_tools": ["spades", "flye", "unicycler", "megahit", "hifiasm"],
        "input_formats": [".fastq", ".fq", ".fastq.gz", ".fq.gz"],
        "description": "Assemble sequencing reads into contigs and scaffolds",
    },
    "genome_annotation": {
        "display_name": "Genome Annotation",
        "aliases": [
            "annotation", "annotate genome", "gene prediction",
            "gene finding", "orf prediction", "prokaryotic annotation",
        ],
        "sections": ["Annotation"],
        "preferred_tools": ["prokka", "augustus", "maker", "bakta"],
        "input_formats": [".fasta", ".fa", ".gbk", ".gff"],
        "description": "Predict and annotate genes in assembled genomes",
    },
    "phylogenetics": {
        "display_name": "Phylogenetics / Evolution",
        "aliases": [
            "phylogenetics", "phylogeny", "evolutionary analysis",
            "tree building", "multiple alignment", "phylogenetic tree",
            "molecular evolution", "msa",
        ],
        "sections": ["Multiple Alignments", "Evolution", "HyPhy"],
        "preferred_tools": ["iqtree", "raxml", "mafft", "muscle", "clustalw"],
        "input_formats": [".fasta", ".fa", ".phy", ".nex"],
        "description": "Construct phylogenetic trees and analyse molecular evolution",
    },
    "chip_seq": {
        "display_name": "ChIP-seq / Epigenetics",
        "aliases": [
            "chip-seq", "chipseq", "peak calling", "epigenetics",
            "histone modification", "transcription factor binding",
            "methylation", "atac-seq", "atacseq",
        ],
        "sections": ["ChIP-seq", "Epigenetics", "deepTools"],
        "preferred_tools": ["macs2", "deeptools", "diffbind"],
        "input_formats": [".bam", ".bed", ".bigwig", ".bw"],
        "description": "Analyse protein-DNA interactions and epigenetic modifications",
    },
    "single_cell": {
        "display_name": "Single-cell Analysis",
        "aliases": [
            "single-cell", "single cell", "scrna", "scrna-seq",
            "10x genomics", "cell clustering", "trajectory",
            "cell type", "marker genes",
        ],
        "sections": ["Single-cell", "HCA-Scanpy", "Seurat", "Monocle3"],
        "preferred_tools": ["scanpy", "cellranger", "seurat"],
        "input_formats": [".h5ad", ".h5", ".mtx", ".loom"],
        "description": "Analyse single-cell RNA-seq data for cell type identification",
    },
    "proteomics": {
        "display_name": "Proteomics",
        "aliases": [
            "proteomics", "mass spectrometry", "protein identification",
            "peptide search", "ms/ms", "protein quantification",
        ],
        "sections": ["Proteomics"],
        "preferred_tools": ["maxquant", "searchgui", "peptideshaker"],
        "input_formats": [".mzml", ".mzxml", ".raw", ".mgf"],
        "description": "Identify and quantify proteins from mass spectrometry data",
    },
    "nanopore": {
        "display_name": "Nanopore / Long-read Sequencing",
        "aliases": [
            "nanopore", "long reads", "oxford nanopore", "ont",
            "pacbio", "long-read", "basecalling", "hifi",
        ],
        "sections": ["Nanopore"],
        "preferred_tools": ["nanoplot", "medaka", "nanopolish", "minimap2", "flye"],
        "input_formats": [".fastq", ".fast5", ".pod5", ".bam"],
        "description": "Process and analyse long-read sequencing data",
    },
    "bam_processing": {
        "display_name": "BAM/SAM Processing",
        "aliases": [
            "bam", "sam", "sort bam", "index bam", "merge bam",
            "filter bam", "markduplicates", "bam stats",
            "sambamba", "duplicate", "duplicates",
        ],
        "sections": ["SAM/BAM", "Picard"],
        "preferred_tools": ["samtools", "picard", "bamtools", "bedtools"],
        "input_formats": [".bam", ".sam", ".cram"],
        "description": "Manipulate, filter, and analyse aligned read files",
    },
}


# ---------------------------------------------------------------------------
# Input format → compatible task categories
# ---------------------------------------------------------------------------

FORMAT_TO_TASKS: dict[str, list[str]] = {
    ".fastq":    ["quality_control", "read_mapping", "metagenomics", "genome_assembly", "nanopore", "rnaseq", "wes_wgs_analysis"],
    ".fq":       ["quality_control", "read_mapping", "metagenomics", "genome_assembly", "nanopore", "rnaseq", "wes_wgs_analysis"],
    ".fastq.gz": ["quality_control", "read_mapping", "metagenomics", "genome_assembly", "nanopore", "rnaseq", "wes_wgs_analysis"],
    ".fq.gz":    ["quality_control", "read_mapping", "metagenomics", "genome_assembly", "nanopore", "rnaseq", "wes_wgs_analysis"],
    ".fasta":    ["read_mapping", "genome_annotation", "phylogenetics", "metagenomics"],
    ".fa":       ["read_mapping", "genome_annotation", "phylogenetics", "metagenomics"],
    ".bam":      ["variant_calling", "bam_processing", "rnaseq", "chip_seq", "wes_wgs_analysis"],
    ".sam":      ["bam_processing"],
    ".cram":     ["variant_calling", "bam_processing", "wes_wgs_analysis"],
    ".vcf":      ["variant_annotation", "wes_wgs_analysis"],
    ".vcf.gz":   ["variant_annotation", "wes_wgs_analysis"],
    ".bcf":      ["variant_annotation"],
    ".bed":      ["chip_seq", "bam_processing"],
    ".bigwig":   ["chip_seq"],
    ".bw":       ["chip_seq"],
    ".h5ad":     ["single_cell"],
    ".h5":       ["single_cell"],
    ".mtx":      ["single_cell"],
    ".loom":     ["single_cell"],
    ".gff":      ["genome_annotation"],
    ".gtf":      ["rnaseq", "genome_annotation"],
    ".mzml":     ["proteomics"],
    ".mzxml":    ["proteomics"],
    ".mgf":      ["proteomics"],
    ".fast5":    ["nanopore"],
    ".pod5":     ["nanopore"],
    ".phy":      ["phylogenetics"],
    ".nex":      ["phylogenetics"],
    ".tsv":      ["rnaseq"],
    ".csv":      ["rnaseq"],
}


# ---------------------------------------------------------------------------
# Workflow templates
# ---------------------------------------------------------------------------

WORKFLOW_TEMPLATES: dict[str, dict] = {
    "wes_germline": {
        "name": "WES/WGS Germline Variant Analysis",
        "description": "Full pipeline from raw exome/genome reads to annotated pathogenic variants",
        "aliases": [
            "wes", "exome", "whole exome", "wgs", "whole genome",
            "germline", "clinical exome", "novogene", "exome sequencing",
            "variant analysis", "pathogenic variants",
        ],
        "input_format": ".fastq",
        "steps": [
            {"tool": "fastqc",       "purpose": "Quality assessment of raw reads"},
            {"tool": "trimmomatic",   "purpose": "Adapter removal and quality trimming"},
            {"tool": "bwa_mem2",      "purpose": "Align reads to hg38 reference genome"},
            {"tool": "samtools",      "purpose": "Sort, index, and mark duplicates"},
            {"tool": "gatk4",         "purpose": "GATK HaplotypeCaller for germline variant calling"},
            {"tool": "snpeff",        "purpose": "Functional annotation of variants (missense, stopgain, splicing)"},
            {"tool": "snpsift",       "purpose": "Filter by ClinVar, gnomAD frequency, CADD score"},
        ],
    },
    "wes_annotation": {
        "name": "WES Variant Annotation & Prioritisation",
        "description": "Annotate and prioritise variants from existing VCF files (e.g. Novogene GATK output)",
        "aliases": [
            "annotate vcf", "variant annotation", "filter variants",
            "pathogenic", "clinvar", "gnomad", "cadd",
            "prioritise variants", "exome annotation", "annovar",
            "missense", "stopgain", "frameshift",
        ],
        "input_format": ".vcf",
        "steps": [
            {"tool": "bcftools",    "purpose": "Normalise and merge multi-sample VCFs"},
            {"tool": "snpeff",      "purpose": "Predict functional effect (missense, nonsense, splicing)"},
            {"tool": "snpsift",     "purpose": "Filter: gnomAD AF < 0.01, ClinVar pathogenic/likely pathogenic"},
            {"tool": "gemini",      "purpose": "Query variants by inheritance model (dominant, recessive, de novo)"},
            {"tool": "vep",         "purpose": "Ensembl VEP: SIFT, PolyPhen, CADD scores"},
        ],
    },
    "rnaseq_de": {
        "name": "RNA-seq Differential Expression",
        "description": "Standard RNA-seq pipeline from raw reads to differentially expressed genes",
        "aliases": [
            "rna-seq", "rnaseq", "differential expression", "gene expression",
            "transcriptomics", "deseq2",
        ],
        "input_format": ".fastq",
        "steps": [
            {"tool": "fastqc",        "purpose": "Quality assessment of raw reads"},
            {"tool": "trimmomatic",    "purpose": "Adapter removal and quality trimming"},
            {"tool": "hisat2",         "purpose": "Splice-aware alignment to reference genome"},
            {"tool": "featurecounts",  "purpose": "Count reads per gene"},
            {"tool": "deseq2",         "purpose": "Statistical test for differential expression"},
        ],
    },
    "metagenomics_profiling": {
        "name": "Metagenomics Community Profiling",
        "description": "Classify microbial species and estimate relative abundance",
        "aliases": [
            "metagenomics", "microbiome", "taxonomic classification",
            "species identification", "microbial community", "16s",
        ],
        "input_format": ".fastq",
        "steps": [
            {"tool": "fastqc",   "purpose": "Quality assessment of raw reads"},
            {"tool": "fastp",    "purpose": "Quality filtering and adapter removal"},
            {"tool": "kraken2",  "purpose": "Taxonomic classification of reads"},
            {"tool": "bracken",  "purpose": "Bayesian re-estimation of species abundance"},
        ],
    },
    "variant_calling_pipeline": {
        "name": "Germline Variant Calling from BAM",
        "description": "Discover SNPs and indels from aligned reads (BAM input)",
        "aliases": [
            "variant calling", "snp calling", "find variants",
            "germline", "genotyping", "call variants from bam",
        ],
        "input_format": ".bam",
        "steps": [
            {"tool": "samtools",  "purpose": "Verify BAM is sorted and indexed"},
            {"tool": "freebayes", "purpose": "Bayesian variant calling"},
            {"tool": "bcftools",  "purpose": "Filter low-quality variants"},
            {"tool": "snpeff",    "purpose": "Functional annotation"},
            {"tool": "snpsift",   "purpose": "Filter by frequency and clinical significance"},
        ],
    },
    "chipseq_pipeline": {
        "name": "ChIP-seq Peak Calling",
        "description": "Identify protein-DNA binding sites from ChIP-seq data",
        "aliases": [
            "chip-seq", "chipseq", "peak calling", "histone modification",
            "transcription factor",
        ],
        "input_format": ".fastq",
        "steps": [
            {"tool": "fastqc",     "purpose": "Quality assessment of raw reads"},
            {"tool": "trimmomatic", "purpose": "Adapter removal and quality trimming"},
            {"tool": "bowtie2",    "purpose": "Align reads to reference genome"},
            {"tool": "macs2",      "purpose": "Peak calling for enriched regions"},
            {"tool": "deeptools",  "purpose": "Signal visualisation and normalisation"},
        ],
    },
    "nanopore_assembly": {
        "name": "Nanopore Long-read Assembly",
        "description": "Assemble and polish a genome from Oxford Nanopore reads",
        "aliases": [
            "nanopore", "long read assembly", "ont assembly", "minion",
        ],
        "input_format": ".fastq",
        "steps": [
            {"tool": "nanoplot",  "purpose": "Long-read quality statistics"},
            {"tool": "flye",      "purpose": "De novo assembly from long reads"},
            {"tool": "medaka",    "purpose": "Polish assembly with neural network"},
            {"tool": "quast",     "purpose": "Assembly quality assessment"},
        ],
    },
    "genome_assembly_pipeline": {
        "name": "De Novo Genome Assembly (Short Reads)",
        "description": "Assemble a genome from Illumina short reads",
        "aliases": [
            "assembly", "de novo", "assemble genome", "scaffold",
        ],
        "input_format": ".fastq",
        "steps": [
            {"tool": "fastqc",  "purpose": "Quality assessment of raw reads"},
            {"tool": "fastp",   "purpose": "Quality filtering"},
            {"tool": "spades",  "purpose": "De novo assembly"},
            {"tool": "quast",   "purpose": "Assembly quality assessment"},
            {"tool": "prokka",  "purpose": "Gene annotation of assembled contigs"},
        ],
    },
}


# ---------------------------------------------------------------------------
# Version deduplication
# ---------------------------------------------------------------------------


def _tool_base_key(tool_id: str) -> str:
    """Extract version-independent key from a Galaxy tool ID.

    'toolshed.g2.bx.psu.edu/repos/devteam/fastqc/fastqc/0.74+galaxy1'
    → 'toolshed.g2.bx.psu.edu/repos/devteam/fastqc/fastqc'
    """
    if not tool_id.startswith("toolshed."):
        return tool_id
    parts = tool_id.strip("/").split("/")
    if len(parts) >= 6:
        return "/".join(parts[:-1])
    return tool_id


def _version_sort_key(version_str: str) -> tuple:
    """Parse '0.74+galaxy1' into sortable tuple."""
    base, _, galaxy_suffix = version_str.partition("+galaxy")
    galaxy_rev = int(galaxy_suffix) if galaxy_suffix.isdigit() else 0
    segments = []
    for part in re.split(r"[.\-]", base):
        try:
            segments.append(int(part))
        except ValueError:
            segments.append(0)
    while len(segments) < 4:
        segments.append(0)
    segments.append(galaxy_rev)
    return tuple(segments)


def _tool_slug(tool_id: str) -> str:
    """Extract short tool name from a Galaxy tool ID.

    'toolshed.../repos/devteam/fastqc/fastqc/0.74' → 'fastqc'
    """
    parts = tool_id.strip("/").split("/")
    if len(parts) >= 2:
        return parts[-2].lower()
    return parts[-1].lower()


def build_deduped_index(tools: list[dict]) -> dict[str, dict]:
    """Group tools by base key, keep only the latest version of each.

    Returns dict mapping base_key → latest tool entry (with version_count).
    """
    groups: dict[str, list[dict]] = {}
    for tool in tools:
        key = _tool_base_key(tool.get("id", ""))
        groups.setdefault(key, []).append(tool)

    index = {}
    for key, group in groups.items():
        best = max(group, key=lambda t: _version_sort_key(t.get("version", "0")))
        best = dict(best)  # don't mutate original
        best["version_count"] = len(group)
        best["all_versions"] = [
            t.get("version", "")
            for t in sorted(group, key=lambda t: _version_sort_key(t.get("version", "0")), reverse=True)
        ]
        index[key] = best
    return index


# ---------------------------------------------------------------------------
# Task matching
# ---------------------------------------------------------------------------


def _match_tasks(description: str) -> list[tuple[str, float]]:
    """Match a natural-language description against the task taxonomy.

    Returns list of (task_id, confidence) sorted by confidence desc.
    """
    desc_lower = description.lower()
    desc_words = set(desc_lower.split())

    matches: list[tuple[str, float]] = []

    for task_id, task in TASK_TAXONOMY.items():
        score = 0.0

        # Check aliases (phrase matching)
        for alias in task["aliases"]:
            if alias in desc_lower:
                score += len(alias.split()) * 2.0

        # Check display name
        if task["display_name"].lower() in desc_lower:
            score += 5.0

        # Word overlap with aliases
        alias_words = set()
        for alias in task["aliases"]:
            alias_words.update(alias.split())
        overlap = desc_words & alias_words
        score += len(overlap) * 0.5

        if score > 0:
            matches.append((task_id, score))

    matches.sort(key=lambda x: -x[1])
    return matches


# ---------------------------------------------------------------------------
# EDAM resolution
# ---------------------------------------------------------------------------


def _resolve_edam(tool: dict) -> list[str]:
    """Resolve EDAM IDs to human-readable labels."""
    labels = []
    for tid in tool.get("edam_topics", []):
        label = EDAM_TOPICS.get(tid)
        if label:
            labels.append(label)
    for oid in tool.get("edam_operations", []):
        label = EDAM_OPERATIONS.get(oid)
        if label:
            labels.append(label)
    return labels


# ---------------------------------------------------------------------------
# Input format detection
# ---------------------------------------------------------------------------


def detect_format(filename: str) -> str | None:
    """Detect format from filename, handling double extensions."""
    name = filename.lower()
    for ext in [".fastq.gz", ".fq.gz", ".vcf.gz"]:
        if name.endswith(ext):
            return ext
    from pathlib import PurePath
    return PurePath(name).suffix or None


# ---------------------------------------------------------------------------
# Multi-signal scoring
# ---------------------------------------------------------------------------


def _score_tool(
    tool: dict,
    description: str,
    matched_tasks: list[tuple[str, float]],
    input_format: str | None,
) -> tuple[float, str]:
    """Score a single tool against task + description + format.

    Signals (weights):
      1. Section matches task taxonomy     (30 pts max)
      2. Preferred tool bonus              (20 pts)
      3. Exact tool name match             (15 pts)
      4. Keyword match in name/desc        (15 pts max)
      5. EDAM topic/operation match        (10 pts max)
      6. Input format compatibility        (10 pts)
      7. Version maturity bonus            (5 pts max)
    """
    score = 0.0
    reasons = []

    tool_name = (tool.get("name") or "").lower()
    tool_desc = (tool.get("description") or "").lower()
    tool_section = (tool.get("section") or "").lower()
    tool_id = (tool.get("id") or "").lower()
    tool_sl = _tool_slug(tool_id)
    desc_lower = description.lower()

    max_task_conf = max((t[1] for t in matched_tasks), default=1.0)

    # 1. Section matches task taxonomy
    for task_id, task_conf in matched_tasks[:3]:
        task = TASK_TAXONOMY[task_id]
        for section in task["sections"]:
            if section.lower() == tool_section:
                bonus = 30.0 * (task_conf / max_task_conf)
                score += bonus
                reasons.append(f"Category: {task['display_name']}")
                break

    # 2. Preferred tool bonus
    for task_id, task_conf in matched_tasks[:3]:
        task = TASK_TAXONOMY[task_id]
        if tool_sl in task["preferred_tools"]:
            score += 20.0
            reasons.append(f"Recommended for {task['display_name']}")
            break

    # 3. Exact tool name match
    if tool_name and tool_name in desc_lower:
        score += 15.0
        reasons.append("Exact name match")

    # 4. Keyword match in name/desc
    kw_score = 0.0
    for word in desc_lower.split():
        if len(word) < 3:
            continue
        if word in tool_name:
            kw_score += 5.0
        elif word in tool_desc:
            kw_score += 2.0
    score += min(kw_score, 15.0)

    # 5. EDAM match
    edam_topics_str = " ".join(
        EDAM_TOPICS.get(t, "").lower() for t in tool.get("edam_topics", [])
    )
    edam_ops_str = " ".join(
        EDAM_OPERATIONS.get(o, "").lower() for o in tool.get("edam_operations", [])
    )
    edam_score = 0.0
    for word in desc_lower.split():
        if len(word) < 3:
            continue
        if word in edam_topics_str:
            edam_score += 3.0
        if word in edam_ops_str:
            edam_score += 3.0
    if edam_score > 0:
        score += min(edam_score, 10.0)
        reasons.append("EDAM ontology match")

    # 6. Input format compatibility
    if input_format:
        format_tasks = FORMAT_TO_TASKS.get(input_format.lower(), [])
        for task_id in format_tasks:
            task = TASK_TAXONOMY.get(task_id, {})
            for section in task.get("sections", []):
                if section.lower() == tool_section:
                    score += 10.0
                    reasons.append(f"Accepts {input_format}")
                    break
            else:
                continue
            break

    # 7. Version maturity (log scale, capped at 5)
    vc = tool.get("version_count", 1)
    if vc > 1:
        maturity = min(5.0, math.log2(vc) * 1.5)
        score += maturity
        if vc >= 5:
            reasons.append(f"Mature ({vc} versions)")

    explanation = "; ".join(reasons[:4]) if reasons else "Keyword match"
    return score, explanation


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

# Cache for deduped index
_deduped_cache: dict[str, dict] | None = None


def recommend_tool(
    task_description: str,
    catalog: dict,
    input_format: str | None = None,
    max_results: int = 5,
) -> list[dict]:
    """Recommend Galaxy tools for a task description.

    Args:
        task_description: Natural language, e.g. "I need to annotate variants from WES"
        catalog: The loaded galaxy_catalog.json dict
        input_format: Optional file extension like ".fastq" or ".vcf"
        max_results: Number of recommendations to return

    Returns:
        List of dicts with: tool_id, name, version, section, description,
        score, explanation, edam_labels, version_count
    """
    global _deduped_cache
    if _deduped_cache is None:
        _deduped_cache = build_deduped_index(catalog.get("tools", []))
    deduped = _deduped_cache

    # Match task
    matched_tasks = _match_tasks(task_description)

    # If input_format provided but no task matched, infer from format
    if not matched_tasks and input_format:
        format_tasks = FORMAT_TO_TASKS.get(input_format.lower(), [])
        matched_tasks = [(tid, 1.0) for tid in format_tasks]

    # Score each deduped tool
    scored: list[tuple[float, dict, str]] = []
    for _base_key, tool in deduped.items():
        s, explanation = _score_tool(tool, task_description, matched_tasks, input_format)
        if s > 0:
            scored.append((s, tool, explanation))

    scored.sort(key=lambda x: -x[0])

    results = []
    for s, tool, explanation in scored[:max_results]:
        results.append({
            "tool_id": tool["id"],
            "name": tool.get("name", "?"),
            "version": tool.get("version", "?"),
            "section": tool.get("section", ""),
            "description": tool.get("description", ""),
            "score": round(s, 2),
            "explanation": explanation,
            "edam_labels": _resolve_edam(tool),
            "version_count": tool.get("version_count", 1),
        })
    return results


def suggest_workflow(
    task_description: str,
    input_format: str | None = None,
) -> list[dict]:
    """Match task description to pre-defined workflow templates.

    Returns list of matching workflow dicts with name, steps, match_score.
    """
    desc_lower = task_description.lower()
    results = []

    for wf_id, wf in WORKFLOW_TEMPLATES.items():
        score = 0.0
        for alias in wf["aliases"]:
            if alias in desc_lower:
                score += len(alias.split()) * 2.0

        # Format compatibility bonus
        if input_format and wf.get("input_format") == input_format:
            score += 3.0

        if score > 0:
            results.append({
                "workflow_id": wf_id,
                "name": wf["name"],
                "description": wf["description"],
                "steps": wf["steps"],
                "match_score": round(score, 2),
            })

    results.sort(key=lambda x: -x["match_score"])
    return results


def get_task_categories() -> list[dict]:
    """Return all task categories for display."""
    return [
        {
            "id": tid,
            "name": task["display_name"],
            "description": task["description"],
            "sections": task["sections"],
            "input_formats": task["input_formats"],
        }
        for tid, task in TASK_TAXONOMY.items()
    ]
