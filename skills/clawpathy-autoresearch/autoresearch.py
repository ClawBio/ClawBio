#!/usr/bin/env python3
"""clawpathy-autoresearch: iterative skill improvement loop.

Inspired by Karpathy's autoresearch. Instead of modifying train.py,
this loop modifies SKILL.md files to optimise agent performance on
a defined task.

Usage:
    python autoresearch.py --task tasks/gwas_reproduction/task.yaml \
        --output /tmp/autoresearch --iterations 20
    python autoresearch.py --demo --output /tmp/autoresearch_demo
"""
from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

# ---------------------------------------------------------------------------
# Ensure project root is importable
# ---------------------------------------------------------------------------
_PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from skills.clawpathy_autoresearch.task import TaskDefinition, load_task
from skills.clawpathy_autoresearch.scorer import reproduction_error
from skills.clawpathy_autoresearch.plotter import plot_progress, ExperimentRecord
from skills.clawpathy_autoresearch.skill_manager import SkillManager


@dataclass
class ExperimentResult:
    """Result of a single experiment iteration.

    Score is mean reproduction error (lower is better, 0 = perfect).
    """

    experiment: int
    score: float  # mean reproduction error, lower is better
    kept: bool
    label: str
    skill_changes: list[str]
    skills_created: list[str]
    per_paper_errors: dict[str, float]
    error_breakdown: dict[str, float]
    timestamp: str = ""

    def __post_init__(self):
        if not self.timestamp:
            self.timestamp = datetime.now(timezone.utc).isoformat()

    def to_dict(self) -> dict[str, Any]:
        return {
            "experiment": self.experiment,
            "score": self.score,
            "kept": self.kept,
            "label": self.label,
            "skill_changes": self.skill_changes,
            "skills_created": self.skills_created,
            "per_paper_errors": self.per_paper_errors,
            "error_breakdown": self.error_breakdown,
            "timestamp": self.timestamp,
        }


def save_experiment_log(results: list[ExperimentResult], path: Path) -> None:
    """Save experiment history to JSON."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    data = [r.to_dict() for r in results]
    path.write_text(json.dumps(data, indent=2))


def load_experiment_log(path: Path) -> list[dict]:
    """Load experiment history from JSON."""
    return json.loads(Path(path).read_text())


class AutoResearchLoop:
    """The core iterative skill improvement loop.

    Score = mean reproduction error. Lower is better. 0 = perfect.
    """

    def __init__(
        self,
        task_path: Path,
        output_dir: Path,
        max_iterations: int = 80,
    ) -> None:
        self.task = load_task(task_path)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.max_iterations = max_iterations
        self.skill_mgr = SkillManager(Path(self.task.skills_dir))
        self.history: list[ExperimentResult] = []
        self.best_score: float = float("inf")  # lower is better

    def run_agent_on_paper(self, paper_id: str) -> dict[str, Any]:
        """Run the agent on a single paper using current skills.

        Override this method to plug in different agent backends.
        """
        return {
            "variants_found": [],
            "total_loci_reported": 0,
        }

    def propose_skill_changes(
        self, feedback: str, current_score: float
    ) -> tuple[str, list[str], list[str]]:
        """Ask the agent to propose skill modifications based on feedback.

        Returns (label, modified_skill_names, created_skill_names).
        Override this for real agent integration.
        """
        return "no changes proposed", [], []

    def run_iteration(self, iteration: int) -> ExperimentResult:
        """Execute one iteration of the loop."""
        snapshot = self.skill_mgr.snapshot()

        paper_errors = {}
        total_breakdown = {}

        for paper in self.task.papers:
            agent_output = self.run_agent_on_paper(paper.id)
            error = reproduction_error(paper.ground_truth or {}, agent_output)
            paper_errors[paper.id] = error.total

            # Accumulate breakdowns
            for k, v in error.to_dict().items():
                if k != "total":
                    total_breakdown[k] = total_breakdown.get(k, 0.0) + v

        n = max(len(self.task.papers), 1)
        avg_error = sum(paper_errors.values()) / n
        avg_breakdown = {k: v / n for k, v in total_breakdown.items()}

        # Lower is better: keep if error decreased
        kept = avg_error < self.best_score
        if kept:
            self.best_score = avg_error
            label = "baseline" if iteration == 1 else "improvement"
        else:
            self.skill_mgr.restore(snapshot)
            label = "reverted"

        diff = self.skill_mgr.diff_from_snapshot(snapshot) if kept else {}
        modified = [k for k, v in diff.items() if v == "modified"]
        created = [k for k, v in diff.items() if v == "created"]

        result = ExperimentResult(
            experiment=iteration,
            score=avg_error,
            kept=kept,
            label=label,
            skill_changes=modified,
            skills_created=created,
            per_paper_errors=paper_errors,
            error_breakdown=avg_breakdown,
        )
        self.history.append(result)
        return result

    def run(self) -> list[ExperimentResult]:
        """Run the full loop for max_iterations."""
        print(f"Starting autoresearch: {self.task.name}")
        print(f"Papers: {self.task.paper_count}, Iterations: {self.max_iterations}")
        print(f"Skills dir: {self.task.skills_dir}")
        print("-" * 60)

        for i in range(1, self.max_iterations + 1):
            result = self.run_iteration(i)
            status = "KEPT" if result.kept else "DISCARDED"
            print(
                f"[{i}/{self.max_iterations}] Error: {result.score:.4f} "
                f"({status}) — {result.label}"
            )

            save_experiment_log(
                self.history, self.output_dir / "experiment_log.json"
            )
            records = [
                ExperimentRecord(
                    experiment=r.experiment,
                    score=r.score,
                    kept=r.kept,
                    label=r.label,
                )
                for r in self.history
            ]
            plot_progress(records, self.output_dir / "progress.png")

        print("-" * 60)
        print(f"Done. Best error: {self.best_score:.4f}")
        print(f"Results: {self.output_dir}")
        return self.history


def run_demo(output_dir: Path) -> None:
    """Run a demo with synthetic data to show a descending error plot.

    Simulates 83 experiments like Karpathy's autoresearch, with error
    starting high (~1.0) and being driven towards zero.
    """
    import random

    random.seed(42)

    history: list[ExperimentRecord] = []
    best = 0.98  # starting error (high)

    labels_kept = [
        "baseline",
        "added neg_log10_p normalisation to gwas-lookup",
        "created variant-resolution skill",
        "added OR confidence-interval cross-check",
        "gwas-lookup: chain to fine-mapping for lead SNPs",
        "expanded PheWAS cross-check workflow",
        "effect-direction: use beta sign not OR",
        "locus-count: parse supplementary tables",
        "added allele-frequency population matching",
        "variant-resolution: fuzzy rsID matching (merged SNPs)",
        "gwas-lookup: add EBI GWAS Catalog fallback",
        "fine-mapping: SuSiE credible set → lead variant",
        "created multi-ancestry reconciliation skill",
        "locus-count: exclude HLA region duplicates",
        "p-value: handle extreme values (1e-300+) via log",
        "OR: cap error at 2.0 for protective alleles",
        "added proxy-SNP LD lookup (r2 > 0.8)",
        "effect-freq: gnomAD v4 supersedes ExAC",
        "variant-resolution: handle tri-allelic sites",
        "direction: infer from beta when OR missing",
        "created GWAS-catalog-scraper skill",
        "locus-count: merge LD-linked loci (r2 > 0.1)",
        "p-value: use exact -log10 from summary stats",
        "fine-mapping: weight by posterior probability",
        "OR: impute from beta + SE when not reported",
        "added ancestry-aware allele frequency skill",
        "variant-resolution: dbSNP merge history lookup",
        "p-value: condition on lead SNP for secondary",
        "created effect-size-harmoniser skill",
        "locus-count: count independent signals not loci",
    ]
    labels_disc = [
        "removed validation step (broke OR checks)",
        "aggressive gotcha pruning",
        "reordered workflow (lost chaining context)",
        "over-specified trigger conditions",
        "merged unrelated skills (gwas + prs)",
        "stripped safety checks",
        "added redundant API calls (rate limited)",
        "weakened scoring rubric",
        "skipped allele-frequency normalisation",
        "removed LD-based deduplication",
        "used p-value instead of -log10(p)",
        "hardcoded EUR frequencies (ancestry bias)",
        "dropped protective variant handling",
    ]

    kept_idx = 0
    disc_idx = 0

    for i in range(1, 84):
        # Higher improvement rate early, diminishing returns later
        p_improve = max(0.12, 0.55 - i * 0.005)
        improves = random.random() < p_improve

        if improves and kept_idx < len(labels_kept):
            # Phase-based reductions: big early, moderate mid, smaller late
            if i < 15:
                reduction = random.uniform(0.04, 0.10)
            elif i < 35:
                reduction = random.uniform(0.02, 0.06)
            elif i < 60:
                reduction = random.uniform(0.01, 0.03)
            else:
                reduction = random.uniform(0.005, 0.015)
            score = best - reduction
            score = max(score, 0.008)  # floor near zero
            label = labels_kept[kept_idx]
            kept_idx += 1
            best = score
            history.append(ExperimentRecord(i, score, True, label))
        else:
            # Failed attempt: error goes up from current best
            score = best + random.uniform(0.002, 0.02)
            label = labels_disc[disc_idx % len(labels_disc)]
            disc_idx += 1
            history.append(ExperimentRecord(i, score, False, label))

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    plot_progress(history, output_dir / "progress.png")

    log = [
        {
            "experiment": r.experiment,
            "score": round(r.score, 4),
            "kept": r.kept,
            "label": r.label,
        }
        for r in history
    ]
    (output_dir / "experiment_log.json").write_text(json.dumps(log, indent=2))

    n_kept = sum(1 for r in history if r.kept)
    print(f"Demo complete. {len(history)} experiments, {n_kept} kept.")
    print(f"Final error: {best:.4f}")
    print(f"Plot: {output_dir / 'progress.png'}")
    print(f"Log: {output_dir / 'experiment_log.json'}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="clawpathy-autoresearch: iterative skill improvement loop"
    )
    parser.add_argument("--task", type=Path, help="Path to task.yaml")
    parser.add_argument("--output", type=Path, default=Path("/tmp/autoresearch"))
    parser.add_argument("--iterations", type=int, default=80)
    parser.add_argument("--demo", action="store_true", help="Run with synthetic demo data")
    args = parser.parse_args()

    if args.demo:
        run_demo(args.output)
        return

    if not args.task:
        parser.error("--task is required (or use --demo)")

    loop = AutoResearchLoop(
        task_path=args.task,
        output_dir=args.output,
        max_iterations=args.iterations,
    )
    loop.run()


if __name__ == "__main__":
    main()
