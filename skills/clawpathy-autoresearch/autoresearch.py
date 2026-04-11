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
import time
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import yaml

# ---------------------------------------------------------------------------
# Ensure project root is importable
# ---------------------------------------------------------------------------
_PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from skills.clawpathy_autoresearch.task import TaskDefinition, load_task
from skills.clawpathy_autoresearch.scorer import (
    automated_score,
    llm_judge_score,
    combine_scores,
)
from skills.clawpathy_autoresearch.plotter import plot_progress, ExperimentRecord
from skills.clawpathy_autoresearch.skill_manager import SkillManager


@dataclass
class ExperimentResult:
    """Result of a single experiment iteration."""

    experiment: int
    score: float
    kept: bool
    label: str
    skill_changes: list[str]
    skills_created: list[str]
    per_paper_scores: dict[str, float]
    automated_score: float
    llm_score: float
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
            "per_paper_scores": self.per_paper_scores,
            "automated_score": self.automated_score,
            "llm_score": self.llm_score,
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
    """The core iterative skill improvement loop."""

    def __init__(
        self,
        task_path: Path,
        output_dir: Path,
        max_iterations: int = 20,
        auto_weight: float = 0.6,
        llm_weight: float = 0.4,
        llm_model: str = "claude-sonnet-4-20250514",
    ) -> None:
        self.task = load_task(task_path)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.max_iterations = max_iterations
        self.auto_weight = auto_weight
        self.llm_weight = llm_weight
        self.llm_model = llm_model
        self.skill_mgr = SkillManager(Path(self.task.skills_dir))
        self.history: list[ExperimentResult] = []
        self.best_score: float = -1.0

    def run_agent_on_paper(self, paper_id: str) -> dict[str, Any]:
        """Run the agent on a single paper using current skills.

        Override this method to plug in different agent backends.
        """
        return {
            "variants_found": [],
            "total_loci_reported": 0,
            "qualitative_summary": "",
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

        paper_scores = {}
        total_auto = 0.0
        total_llm = 0.0

        for paper in self.task.papers:
            agent_output = self.run_agent_on_paper(paper.id)

            auto_result = automated_score(paper.ground_truth or {}, agent_output)
            auto_s = auto_result.automated_total

            gt_text = yaml.dump(paper.ground_truth) if paper.ground_truth else ""
            agent_text = yaml.dump(agent_output)
            llm_s, _ = llm_judge_score(gt_text, agent_text, model=self.llm_model)

            combined = combine_scores(auto_s, llm_s, self.auto_weight, self.llm_weight)
            paper_scores[paper.id] = combined
            total_auto += auto_s
            total_llm += llm_s

        n = max(len(self.task.papers), 1)
        avg_score = sum(paper_scores.values()) / n
        avg_auto = total_auto / n
        avg_llm = total_llm / n

        kept = avg_score > self.best_score
        if kept:
            self.best_score = avg_score
            label = "baseline" if iteration == 1 else "improvement"
        else:
            self.skill_mgr.restore(snapshot)
            label = "reverted"

        diff = self.skill_mgr.diff_from_snapshot(snapshot) if kept else {}
        modified = [k for k, v in diff.items() if v == "modified"]
        created = [k for k, v in diff.items() if v == "created"]

        result = ExperimentResult(
            experiment=iteration,
            score=avg_score,
            kept=kept,
            label=label,
            skill_changes=modified,
            skills_created=created,
            per_paper_scores=paper_scores,
            automated_score=avg_auto,
            llm_score=avg_llm,
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
                f"[{i}/{self.max_iterations}] Score: {result.score:.2f} "
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
        print(f"Done. Best score: {self.best_score:.2f}")
        print(f"Results: {self.output_dir}")
        return self.history


def run_demo(output_dir: Path) -> None:
    """Run a demo with synthetic data to show the plot."""
    import random

    random.seed(42)

    history: list[ExperimentRecord] = []
    best = 2.5
    exp = 0

    labels_kept = [
        "baseline",
        "added gwas-lookup chaining",
        "created fine-mapping skill",
        "expanded PheWAS workflow",
        "added LD cross-reference gotcha",
        "improved variant resolution step",
        "created pathway-enrichment skill",
        "added multi-ancestry check",
        "expanded effect size validation",
        "added cross-database deduplication",
        "refined LLM judge rubric",
        "added confidence interval parsing",
    ]
    labels_disc = [
        "removed validation step",
        "aggressive gotcha pruning",
        "reordered workflow badly",
        "over-specified trigger conditions",
        "merged unrelated skills",
        "stripped safety checks",
        "added redundant API calls",
        "weakened scoring rubric",
    ]

    kept_idx = 0
    disc_idx = 0

    for i in range(1, 51):
        exp += 1
        improves = random.random() < max(0.15, 0.4 - i * 0.005)

        if improves and kept_idx < len(labels_kept):
            improvement = random.uniform(0.3, 1.2) * max(0.5, 1.0 - i * 0.01)
            score = best + improvement
            label = labels_kept[kept_idx]
            kept_idx += 1
            best = score
            history.append(ExperimentRecord(exp, score, True, label))
        else:
            score = best - random.uniform(0.1, 2.0)
            label = labels_disc[disc_idx % len(labels_disc)]
            disc_idx += 1
            history.append(ExperimentRecord(exp, score, False, label))

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    plot_progress(history, output_dir / "progress.png")

    log = [
        {
            "experiment": r.experiment,
            "score": round(r.score, 2),
            "kept": r.kept,
            "label": r.label,
        }
        for r in history
    ]
    (output_dir / "experiment_log.json").write_text(json.dumps(log, indent=2))

    print(f"Demo complete. {len(history)} experiments, {sum(1 for r in history if r.kept)} kept.")
    print(f"Plot: {output_dir / 'progress.png'}")
    print(f"Log: {output_dir / 'experiment_log.json'}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="clawpathy-autoresearch: iterative skill improvement loop"
    )
    parser.add_argument("--task", type=Path, help="Path to task.yaml")
    parser.add_argument("--output", type=Path, default=Path("/tmp/autoresearch"))
    parser.add_argument("--iterations", type=int, default=20)
    parser.add_argument("--demo", action="store_true", help="Run with synthetic demo data")
    parser.add_argument("--auto-weight", type=float, default=0.6)
    parser.add_argument("--llm-weight", type=float, default=0.4)
    parser.add_argument(
        "--llm-model", default="claude-sonnet-4-20250514", help="Model for LLM judge"
    )
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
        auto_weight=args.auto_weight,
        llm_weight=args.llm_weight,
        llm_model=args.llm_model,
    )
    loop.run()


if __name__ == "__main__":
    main()
