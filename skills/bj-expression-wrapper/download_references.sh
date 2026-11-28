#!/usr/bin/env bash
# Download the public inputs bj_reference_builder.py needs into an EMPTY directory.
#   usage: download_references.sh <empty-dir>
# Sources: NCBI (GRCh38 no-alt analysis-set genome) and Ensembl release 104
# (GTF, cDNA, ncRNA). Checksum files are saved beside the data; the builder
# verifies the genome against NCBI's md5.
set -euo pipefail
[ -n "${1:-}" ] && [ -d "$1" ] && [ -z "$(ls -A "$1")" ] || { echo "usage: $0 <existing EMPTY directory>" >&2; exit 2; }
D="$1"; cd "$D"
N=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GRCh38_major_release_seqs_for_alignment_pipelines
E=https://ftp.ensembl.org/pub/release-104
get() { curl -fsSL --retry 3 -C - -o "$2" "$1" && echo "ok $2"; }
get $N/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz genome.fna.gz
get $N/md5checksums.txt ncbi_md5checksums.txt
get $E/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.gtf.gz ensembl104.gtf.gz
get $E/gtf/homo_sapiens/CHECKSUMS ensembl_gtf_CHECKSUMS
get $E/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz ensembl104.cdna.fa.gz
get $E/fasta/homo_sapiens/cdna/CHECKSUMS ensembl_cdna_CHECKSUMS
get $E/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz ensembl104.ncrna.fa.gz
get $E/fasta/homo_sapiens/ncrna/CHECKSUMS ensembl_ncrna_CHECKSUMS
echo ALL_DONE
