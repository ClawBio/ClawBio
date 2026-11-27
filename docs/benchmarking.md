# Validation and benchmarking

ClawBio v0.5.0 introduces the platform's first systematic validation infrastructure. Skills are scored against curated ground truth with objective metrics.

**What's in the benchmark suite:**

| Component | Description |
|-----------|-------------|
| **AD Ground Truth** | 34 Alzheimer's disease genes across 3 evidence tiers (Mendelian, GWAS-replicated, novel), 20 negative controls, 10 lead variants |
| **Mock API Server** | Deterministic endpoints for Ensembl, GWAS Catalog, ClinPGx. Offline CI without rate limits |
| **Benchmark Scorer** | Gene recovery rate, FDR, precision, recall, F1, tier-weighted composite score |
| **Swappable Fine-Mapping** | ABF vs SuSiE head-to-head on same data. Method registry pattern for adding FINEMAP, PolyFun |
| **Nightly Sweep** | Runs every skill in demo mode, scores outputs against ground truth |

```bash
# Run the fine-mapping benchmark
python tests/benchmark/finemapping_benchmark.py --output /tmp/fm_bench

# Score a gene list against AD ground truth
python tests/benchmark/benchmark_scorer.py --genes "APP,BIN1,CLU,TREM2,GAPDH"

# Start mock API server for offline testing
python tests/benchmark/mock_api_server.py &
```

**74 benchmark tests** at v0.5.0 baseline, all green. The public leaderboard now tracks **161 / 182 tests passing (88.5%)** across 10 audited skills, up from 80 / 140 (57.1%) at the original audit. See [benchmarks.html](https://clawbio.ai/benchmarks.html) for the live leaderboard and [CHANGELOG.md](../CHANGELOG.md) for full details.
