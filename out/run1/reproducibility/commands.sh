#!/bin/sh
# Regenerate every number in report.md from a clean checkout.
set -eu

python3 skills/abstention-ledger/fetch_evidence.py \
    --input ../data/challenge1/challenge1-b37-segregation.tsv \
    --output out/vep_grch37_cache.json

python3 skills/abstention-ledger/abstention_ledger.py \
    --input ../data/challenge1/challenge1-b37-segregation.tsv \
    --vcf ../data/challenge1/challenge1-b37-segregation.vcf.gz \
    --evidence out/vep_grch37_cache.json \
    --output ../out/run1
