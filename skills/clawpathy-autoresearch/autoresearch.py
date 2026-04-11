#!/usr/bin/env python3
"""clawpathy-autoresearch: iterative skill improvement loop.

General-purpose meta-skill that creates and iteratively improves SKILL.md
files for any reproducible task. Domain-agnostic: all domain knowledge
enters through the workspace.

Usage:
    python autoresearch.py --task /path/to/workspace --iterations 80
    python autoresearch.py --demo --output /tmp/autoresearch_demo
"""
from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

_PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from skills.clawpathy_autoresearch.workspace import load_workspace, Workspace
from skills.clawpathy_autoresearch.scorer import load_scorer
from skills.clawpathy_autoresearch.plotter import plot_progress, ExperimentRecord
from skills.clawpathy_autoresearch.skill_manager import SkillManager


@dataclass
class ExperimentResult:
    """Result of a single experiment iteration."""

    experiment: int
    score: float
    kept: bool
    label: str
    skill_diff: dict[str, str]
    error_breakdown: dict[str, Any]
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
            "skill_diff": self.skill_diff,
            "error_breakdown": self.error_breakdown,
            "timestamp": self.timestamp,
        }


def save_experiment_log(results: list[ExperimentResult], path: Path) -> None:
    """Save experiment history to JSON."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps([r.to_dict() for r in results], indent=2))


def load_experiment_log(path: Path) -> list[dict]:
    """Load experiment history from JSON."""
    return json.loads(Path(path).read_text())


class AutoResearchLoop:
    """The core iterative skill improvement loop.

    Operates on a workspace directory containing task.json, ground_truth.json,
    scorer.py, and skill/SKILL.md. Score = output of workspace scorer.
    Lower is better. 0 = perfect.
    """

    def __init__(
        self,
        workspace_dir: Path,
        output_dir: Path,
        max_iterations: int | None = None,
        early_stop_n: int | None = None,
    ) -> None:
        self.workspace = load_workspace(workspace_dir)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.max_iterations = max_iterations or self.workspace.max_iterations
        self.early_stop_n = early_stop_n or self.workspace.early_stop_n
        self.score_fn = load_scorer(self.workspace.scorer_path)
        self.skill_mgr = SkillManager(self.workspace.workspace_dir)
        self.history: list[ExperimentResult] = []
        self.best_score: float = float("inf")
        self._consecutive_non_improvements = 0

    def run_agent_on_skill(
        self, skill_path: Path, sources_dir: Path
    ) -> dict[str, Any]:
        """Run an agent that reads the SKILL.md and produces output.

        Override this method to plug in different agent backends.
        Returns a dict matching the ground truth schema for scoring.
        """
        return {}

    def propose_skill_changes(
        self, score: float, error_breakdown: dict[str, Any], history: list[ExperimentResult]
    ) -> str:
        """Propose a single targeted edit to the SKILL.md.

        Override this for real agent integration.
        Returns a label describing the proposed change.
        """
        return "no changes proposed"

    def run_iteration(self, iteration: int) -> ExperimentResult:
        """Execute one iteration of the loop."""
        snapshot = self.skill_mgr.snapshot()

        skill_path = self.workspace.skill_dir / "SKILL.md"
        output = self.run_agent_on_skill(skill_path, self.workspace.sources_dir)
        score = self.score_fn(output, self.workspace.ground_truth)

        kept = score < self.best_score
        if kept:
            self.best_score = score
            self._consecutive_non_improvements = 0
            label = "baseline" if iteration == 1 else "improvement"
        else:
            self.skill_mgr.restore(snapshot)
            self._consecutive_non_improvements += 1
            label = "discarded"

        diff = self.skill_mgr.diff_from_snapshot(snapshot) if kept else {}

        result = ExperimentResult(
            experiment=iteration,
            score=score,
            kept=kept,
            label=label,
            skill_diff=diff,
            error_breakdown={},
        )
        self.history.append(result)
        return result

    def run(self) -> list[ExperimentResult]:
        """Run the full loop with early stopping."""
        print(f"Starting autoresearch: {self.workspace.name}")
        print(f"Max iterations: {self.max_iterations}, Early stop: {self.early_stop_n}")
        print("-" * 60)

        for i in range(1, self.max_iterations + 1):
            result = self.run_iteration(i)
            status = "KEPT" if result.kept else "DISCARDED"
            print(
                f"[{i}/{self.max_iterations}] Score: {result.score:.4f} "
                f"({status}) — {result.label}"
            )

            save_experiment_log(self.history, self.output_dir / "experiment_log.json")
            records = [
                ExperimentRecord(r.experiment, r.score, r.kept, r.label)
                for r in self.history
            ]
            plot_progress(records, self.output_dir / "progress.png")

            if self._consecutive_non_improvements >= self.early_stop_n and i > 1:
                print(f"Early stop: {self.early_stop_n} consecutive non-improvements")
                break

        print("-" * 60)
        print(f"Done. Best score: {self.best_score:.4f} ({len(self.history)} experiments)")
        print(f"Results: {self.output_dir}")
        return self.history


def run_demo(output_dir: Path) -> None:
    """Run a demo with synthetic data to show a descending error plot."""
    import random

    random.seed(42)
    history: list[ExperimentRecord] = []
    best = 0.98

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
        "fine-mapping: SuSiE credible set to lead variant",
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
        p_improve = max(0.08, 0.30 - i * 0.003)
        improves = random.random() < p_improve

        if improves and kept_idx < len(labels_kept):
            if i < 20:
                reduction = random.uniform(0.06, 0.14)
            elif i < 45:
                reduction = random.uniform(0.04, 0.09)
            elif i < 65:
                reduction = random.uniform(0.02, 0.05)
            else:
                reduction = random.uniform(0.01, 0.025)
            score = max(best - reduction, 0.008)
            label = labels_kept[kept_idx]
            kept_idx += 1
            best = score
            history.append(ExperimentRecord(i, score, True, label))
        else:
            jitter = random.choice([
                random.uniform(0.01, 0.06),
                random.uniform(0.06, 0.20),
                random.uniform(0.20, 0.55),
                random.uniform(0.30, 0.70),
                random.uniform(-0.005, 0.02),
            ])
            score = max(best + jitter, 0.005)
            label = labels_disc[disc_idx % len(labels_disc)]
            disc_idx += 1
            history.append(ExperimentRecord(i, score, False, label))

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    plot_progress(history, output_dir / "progress.png")

    log = [
        {"experiment": r.experiment, "score": round(r.score, 4), "kept": r.kept, "label": r.label}
        for r in history
    ]
    (output_dir / "experiment_log.json").write_text(json.dumps(log, indent=2))

    n_kept = sum(1 for r in history if r.kept)
    print(f"Demo complete. {len(history)} experiments, {n_kept} kept.")
    print(f"Final score: {best:.4f}")
    print(f"Plot: {output_dir / 'progress.png'}")
    print(f"Log: {output_dir / 'experiment_log.json'}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="clawpathy-autoresearch: iterative skill improvement loop"
    )
    parser.add_argument(
        "--task", type=Path,
        help="Path to workspace directory (containing task.json, ground_truth.json, scorer.py, skill/SKILL.md)"
    )
    parser.add_argument("--output", type=Path, default=Path("/tmp/autoresearch"))
    parser.add_argument("--iterations", type=int, default=None)
    parser.add_argument("--demo", action="store_true", help="Run with synthetic demo data")
    args = parser.parse_args()

    if args.demo:
        run_demo(args.output)
        return

    if not args.task:
        parser.error("--task is required (path to workspace directory, or use --demo)")

    loop = AutoResearchLoop(
        workspace_dir=args.task,
        output_dir=args.output,
        max_iterations=args.iterations,
    )
    loop.run()


if __name__ == "__main__":
    main()
