#!/usr/bin/env python3
"""
benchmark_scorer.py — Score pipeline outputs against AD ground truth.

Evaluates gene recovery rate, false discovery rate, precision, recall, F1,
and weighted score against the curated AD benchmark set.

Usage:
    from tests.benchmark.benchmark_scorer import BenchmarkScorer

    scorer = BenchmarkScorer()
    results = scorer.score(pipeline_genes=["BIN1", "CLU", "GAPDH", "APP"])
    print(results["f1"], results["gene_recovery_rate"])

CLI:
    python tests/benchmark/benchmark_scorer.py --genes "BIN1,CLU,CR1,GAPDH,APP"
    python tests/benchmark/benchmark_scorer.py --json pipeline_output.json --gene-field "gene_symbol"
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

GROUND_TRUTH_PATH = Path(__file__).parent / "ad_ground_truth.json"


class BenchmarkScorer:
    """Score a set of pipeline-discovered genes against AD ground truth."""

    def __init__(self, ground_truth_path: Path = GROUND_TRUTH_PATH):
        with open(ground_truth_path, encoding="utf-8") as f:
            self.ground_truth = json.load(f)

        pos = self.ground_truth["positive_genes"]
        self.tier1 = {g["gene"] for g in pos["tier1_causal"]}
        self.tier2 = {g["gene"] for g in pos["tier2_gwas_replicated"]}
        self.tier3 = {g["gene"] for g in pos["tier3_novel_bellenguez"]}
        self.all_positive = self.tier1 | self.tier2 | self.tier3
        self.negative = {
            g["gene"] for g in self.ground_truth["negative_genes"]["genes"]
        }

        scoring = self.ground_truth["scoring"]
        self.tier_weights = scoring["tier_weighting"]
        self.minimums = scoring["minimum_acceptable"]

    def score(self, pipeline_genes: list[str]) -> dict:
        """Score a list of genes discovered by a pipeline.

        Returns dict with gene_recovery_rate, precision, recall, f1,
        weighted_score, tier breakdown, and pass/fail status.
        """
        found = set(pipeline_genes)

        # True positives by tier
        tp_tier1 = found & self.tier1
        tp_tier2 = found & self.tier2
        tp_tier3 = found & self.tier3
        tp_all = found & self.all_positive

        # False positives (in negative set)
        fp = found & self.negative

        # Genes not in either set (unknown, not scored)
        unknown = found - self.all_positive - self.negative

        # Metrics
        tp_count = len(tp_all)
        fp_count = len(fp)
        fn_count = len(self.all_positive - found)

        precision = tp_count / (tp_count + fp_count) if (tp_count + fp_count) > 0 else 0.0
        recall = tp_count / len(self.all_positive) if len(self.all_positive) > 0 else 0.0
        f1 = (2 * precision * recall / (precision + recall)) if (precision + recall) > 0 else 0.0

        gene_recovery_rate = recall  # Same as recall for this benchmark
        fdr = fp_count / (tp_count + fp_count) if (tp_count + fp_count) > 0 else 0.0

        # Weighted score (tier1 genes count 3x, tier2 2x, tier3 1x)
        w1 = self.tier_weights["tier1_causal"]
        w2 = self.tier_weights["tier2_gwas_replicated"]
        w3 = self.tier_weights["tier3_novel_bellenguez"]
        weighted_tp = len(tp_tier1) * w1 + len(tp_tier2) * w2 + len(tp_tier3) * w3
        weighted_max = len(self.tier1) * w1 + len(self.tier2) * w2 + len(self.tier3) * w3
        weighted_score = weighted_tp / weighted_max if weighted_max > 0 else 0.0

        # Pass/fail against minimums
        passes = (
            gene_recovery_rate >= self.minimums["gene_recovery_rate"]
            and precision >= self.minimums["precision"]
            and f1 >= self.minimums["f1"]
        )

        return {
            "pipeline_genes_count": len(found),
            "true_positives": tp_count,
            "false_positives": fp_count,
            "false_negatives": fn_count,
            "unknown_genes": len(unknown),
            "gene_recovery_rate": round(gene_recovery_rate, 4),
            "false_discovery_rate": round(fdr, 4),
            "precision": round(precision, 4),
            "recall": round(recall, 4),
            "f1": round(f1, 4),
            "weighted_score": round(weighted_score, 4),
            "tier_breakdown": {
                "tier1_found": sorted(tp_tier1),
                "tier1_missed": sorted(self.tier1 - found),
                "tier2_found": sorted(tp_tier2),
                "tier2_missed": sorted(self.tier2 - found),
                "tier3_found": sorted(tp_tier3),
                "tier3_missed": sorted(self.tier3 - found),
            },
            "false_positive_genes": sorted(fp),
            "unknown_genes_list": sorted(unknown),
            "passes_minimum": passes,
            "minimums": self.minimums,
        }

    def score_variants(self, pipeline_variants: list[dict], rsid_field: str = "rsid") -> dict:
        """Score by lead variant recovery instead of gene recovery.

        pipeline_variants: list of dicts, each with at least an rsid field.
        """
        lead = self.ground_truth.get("lead_variants", {}).get("variants", [])
        lead_rsids = {v["rsid"] for v in lead}
        lead_genes = {v["rsid"]: v["gene"] for v in lead}

        found_rsids = {v[rsid_field] for v in pipeline_variants if rsid_field in v}
        recovered = found_rsids & lead_rsids
        missed = lead_rsids - found_rsids

        recovery_rate = len(recovered) / len(lead_rsids) if lead_rsids else 0.0

        return {
            "lead_variants_total": len(lead_rsids),
            "lead_variants_recovered": len(recovered),
            "lead_variant_recovery_rate": round(recovery_rate, 4),
            "recovered": sorted(recovered),
            "recovered_genes": sorted({lead_genes[r] for r in recovered}),
            "missed": sorted(missed),
            "missed_genes": sorted({lead_genes[r] for r in missed}),
        }

    def summary_markdown(self, result: dict) -> str:
        """Generate a markdown summary of benchmark results."""
        status = "PASS" if result["passes_minimum"] else "FAIL"
        lines = [
            f"# AD Benchmark Results [{status}]",
            "",
            f"**Pipeline genes evaluated**: {result['pipeline_genes_count']}",
            f"**True positives**: {result['true_positives']}",
            f"**False positives**: {result['false_positives']}",
            f"**False negatives**: {result['false_negatives']}",
            "",
            "## Metrics",
            "",
            "| Metric | Value | Minimum |",
            "|--------|-------|---------|",
            f"| Gene recovery rate | {result['gene_recovery_rate']:.4f} | {result['minimums']['gene_recovery_rate']} |",
            f"| Precision | {result['precision']:.4f} | {result['minimums']['precision']} |",
            f"| F1 | {result['f1']:.4f} | {result['minimums']['f1']} |",
            f"| FDR | {result['false_discovery_rate']:.4f} | - |",
            f"| Weighted score | {result['weighted_score']:.4f} | - |",
            "",
            "## Tier Breakdown",
            "",
            f"- **Tier 1 (causal)**: {len(result['tier_breakdown']['tier1_found'])}/{len(result['tier_breakdown']['tier1_found']) + len(result['tier_breakdown']['tier1_missed'])}",
            f"- **Tier 2 (GWAS replicated)**: {len(result['tier_breakdown']['tier2_found'])}/{len(result['tier_breakdown']['tier2_found']) + len(result['tier_breakdown']['tier2_missed'])}",
            f"- **Tier 3 (novel Bellenguez)**: {len(result['tier_breakdown']['tier3_found'])}/{len(result['tier_breakdown']['tier3_found']) + len(result['tier_breakdown']['tier3_missed'])}",
        ]

        if result["false_positive_genes"]:
            lines.append("")
            lines.append(f"## False Positives: {', '.join(result['false_positive_genes'])}")

        lines.append("")
        return "\n".join(lines)


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Score pipeline genes against AD ground truth")
    parser.add_argument("--genes", type=str, help="Comma-separated gene list")
    parser.add_argument("--json", type=str, help="JSON file with pipeline output")
    parser.add_argument("--gene-field", type=str, default="gene", help="Field name for gene in JSON")
    parser.add_argument("--output", type=str, help="Output directory for results")
    args = parser.parse_args()

    scorer = BenchmarkScorer()

    if args.genes:
        genes = [g.strip() for g in args.genes.split(",") if g.strip()]
    elif args.json:
        with open(args.json) as f:
            data = json.load(f)
        if isinstance(data, list):
            genes = [item[args.gene_field] for item in data if args.gene_field in item]
        elif isinstance(data, dict) and "genes" in data:
            genes = data["genes"]
        else:
            print("Cannot extract genes from JSON. Use --gene-field.", file=sys.stderr)
            sys.exit(1)
    else:
        print("Provide --genes or --json", file=sys.stderr)
        sys.exit(1)

    result = scorer.score(genes)
    print(scorer.summary_markdown(result))

    if args.output:
        out = Path(args.output)
        out.mkdir(parents=True, exist_ok=True)
        (out / "benchmark_results.json").write_text(json.dumps(result, indent=2))
        (out / "benchmark_results.md").write_text(scorer.summary_markdown(result))
        print(f"\nResults written to {out}/")


if __name__ == "__main__":
    main()
