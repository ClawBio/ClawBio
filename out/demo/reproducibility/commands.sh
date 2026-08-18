#!/bin/sh
# Regenerate every number in report.md from a clean checkout.
set -eu

python3 skills/abstention-ledger/fetch_evidence.py \
    --input /Users/m/Claude/clawbio/ClawBio/skills/abstention-ledger/examples/demo_segregation.tsv \
    --output out/vep_grch37_cache.json

python3 skills/abstention-ledger/abstention_ledger.py \
    --input /Users/m/Claude/clawbio/ClawBio/skills/abstention-ledger/examples/demo_segregation.tsv \
    --vcf <vcf> \
    --evidence out/vep_grch37_cache.json \
    --output ../out/demo
