# 🦖 Lit Synthesizer

**ClawBio skill** — Search PubMed and bioRxiv, synthesise results, build citation graphs.

Part of the [ClawBio](https://github.com/ClawBio/ClawBio) bioinformatics AI agent skill library.

## Quick Start

```bash
# Demo (no network, no setup needed)
python skills/lit-synthesizer/lit_synthesizer.py --demo --output /tmp/demo

# Real search
python skills/lit-synthesizer/lit_synthesizer.py \
    --query "CRISPR off-target effects" \
    --output report/
```

## What it produces

```
output/
├── report.md                    # Full synthesis report
├── results.json                 # All papers as JSON
├── citation_graph.json          # Node-edge citation graph
├── tables/papers.csv            # Tabular paper list
└── reproducibility/
    ├── commands.sh              # Reproduce this exact run
    ├── environment.yml          # Python environment
    └── checksums.sha256         # SHA-256 of all outputs
```

## Run Tests

```bash
python -m pytest skills/lit-synthesizer/tests/ -v
```

## Author

Contributed by [Sooraj](https://github.com/sooraj-codes) — resolves [Issue #7](https://github.com/ClawBio/ClawBio/issues/7).
