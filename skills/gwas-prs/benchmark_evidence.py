#!/usr/bin/env python3
"""Offline synthetic evidence-policy contract benchmark; not clinical validation."""
from __future__ import annotations

import argparse
import copy
import importlib.util
import json
from pathlib import Path
import tempfile

SKILL = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("prs_evidence", SKILL / "prs_evidence.py")
policy = importlib.util.module_from_spec(spec)
spec.loader.exec_module(policy)

DISCLAIMER = ("ClawBio is a research and educational tool. It is not a medical device "
              "and does not provide clinical diagnoses. Consult a healthcare "
              "professional before making any medical decisions.")


def run_benchmark():
    """Deliberately fabricated evidence exercises policy mechanics only."""
    with tempfile.TemporaryDirectory(prefix="prs-evidence-") as tmp:
        root = Path(tmp)
        score = root / "score.txt"
        inp = root / "input.txt"
        score.write_text("rsID\teffect_allele\tother_allele\teffect_weight\nrs1\tA\tG\t0.2\n")
        inp.write_text("# synthetic only\nrs1\t1\t1\tAG\n")
        context = {"population": "synthetic-population", "context": "synthetic-adults"}
        base = {
            "score_id": "SYNTHETIC", "scoring_path": score, "input_path": inp,
            "build": "GRCh37", "raw_score": 0.2, "genotypes": {"rs1": "AG"},
            "variants": [{"rsid": "rs1", "effect_allele": "A", "other_allele": "G", "effect_weight": 0.2}],
            "evidence": {
                "schema_version": 1,
                "input": dict(context, sha256=policy.digest(inp), build="GRCh37",
                              strand="forward", source="synthetic manifest"),
                "scores": {"SYNTHETIC": {
                    "sha256": policy.digest(score), "build": "GRCh37", "version": "1",
                    "source": "synthetic score", "variant_count": 1,
                    "validation": dict(context, independent=True, n=1000,
                        source="fabricated validation", metric="r2", estimate=0.03, ci_lower=0.01, ci_upper=0.05),
                    "reference": dict(context, n=1000, source="fabricated reference",
                        mean=0.2, sd=0.1, distribution="normal", missingness="complete_only"),
                }},
            },
        }
        cases = [
            ("synthetic_personality_supported", (), None, "supported"),
            ("synthetic_metabolic_supported", ("evidence", "scores", "SYNTHETIC", "trait"), "metabolic", "supported"),
            ("missing_evidence", ("evidence",), None, "not_established"),
            ("missing_score", ("evidence", "scores"), {}, "not_established"),
            ("different_population", ("evidence", "input", "population"), "another", "not_established"),
            ("different_context", ("evidence", "input", "context"), "children", "not_established"),
            ("no_independent_validation", ("evidence", "scores", "SYNTHETIC", "validation", "independent"), False, "not_established"),
            ("missing_reference", ("evidence", "scores", "SYNTHETIC", "reference"), None, "not_established"),
            ("different_weights", ("evidence", "scores", "SYNTHETIC", "sha256"), "0"*64, "incompatible"),
            ("different_build", ("evidence", "input", "build"), "GRCh38", "incompatible"),
            ("missing_variant", ("genotypes",), {}, "not_established"),
            ("invalid_genotype", ("genotypes", "rs1"), "00", "incompatible"),
            ("wrong_sum", ("raw_score",), 999, "incompatible"),
            ("unknown_strand", ("evidence", "input", "strand"), None, "not_established"),
        ]
        rows = []
        for name, path, value, expected in cases:
            args = copy.deepcopy(base)
            if path:
                obj = args
                for key in path[:-1]:
                    obj = obj[key]
                obj[path[-1]] = value
            result = policy.assess(**args)
            rows.append({"case": name, "expected": expected, "observed": result["status"],
                         "percentile": result["percentile"],
                         "reasons": [r["code"] for r in result["reasons"]],
                         "passed": result["status"] == expected
                         and (result["percentile"] == 50 if expected == "supported"
                              else result["percentile"] is None)})
    return {
        "evidence_type": "synthetic_software_contract", "policy_version": policy.POLICY_VERSION,
        "all_passed": all(row["passed"] for row in rows),
        "supported_cases": sum(row["expected"] == "supported" for row in rows),
        "withheld_cases": sum(row["expected"] != "supported" for row in rows),
        "false_supported": sum(row["expected"] != "supported" and row["observed"] == "supported" for row in rows),
        "false_withheld": sum(row["expected"] == "supported" and row["observed"] != "supported" for row in rows),
        "cases": rows,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--demo", action="store_true", required=True)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists() and any(args.output.iterdir()):
        parser.error("Choose a new or empty output directory.")
    result = run_benchmark()
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    lines = ["# PRS evidence policy: synthetic software benchmark", "",
             "Fabricated evidence and weights. Not a reproduction of the Nature paper,",
             "an independently adjudicated benchmark, or evidence of clinical validity.", "",
             "| Case | Expected | Observed | Percentile |", "|---|---|---|---|"]
    for row in result["cases"]:
        lines.append(f"| {row['case']} | {row['expected']} | {row['observed']} | {row['percentile']} |")
    lines += ["", DISCLAIMER, ""]
    (args.output / "report.md").write_text("\n".join(lines))
    print(json.dumps({k: v for k, v in result.items() if k != "cases"}, indent=2))
    return 0 if result["all_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
