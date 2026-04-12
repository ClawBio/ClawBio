"""Stub. The loop engine is rewritten in Task 8."""


def run(*args, **kwargs):
    raise NotImplementedError("rewritten in Task 8")



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
        self.enable_critic: bool = True

    def run_agent_on_skill(
        self, skill_path: Path, sources_dir: Path
    ) -> dict[str, Any]:
        """Run an agent that reads the SKILL.md and produces output.

        MUST be overridden per-workspace. The default raises to prevent
        silent vacuous loops where every iteration scores an empty dict.
        """
        raise NotImplementedError(
            "AutoResearchLoop.run_agent_on_skill must be overridden. "
            "Subclass AutoResearchLoop and implement run_agent_on_skill(skill_path, sources_dir) "
            "to execute the workspace task and return a dict matching ground_truth schema."
        )

    def propose_skill_changes(
        self, score: float, error_breakdown: dict[str, Any], history: list[ExperimentResult]
    ) -> str:
        """Propose a single targeted edit to the SKILL.md via the LLM proposer.

        Default implementation dispatches to `proposer.propose_edit` which
        rewrites SKILL.md in place. Override to plug in a different proposer.
        """
        from skills.clawpathy_autoresearch.proposer import propose_edit

        skill_path = self.workspace.skill_dir / "SKILL.md"
        history_labels = [r.label for r in history]
        result = propose_edit(
            workspace_name=self.workspace.name,
            skill_path=skill_path,
            history_labels=history_labels,
            last_score=score,
        )
        if result.error:
            print(f"  (proposer error: {result.error[:160]})")
        return result.label

    def run_iteration(self, iteration: int) -> ExperimentResult:
        """Execute one iteration of the loop.

        Flow: propose changes (iteration > 1) → snapshot → run → score → keep/revert.
        """
        # Snapshot BEFORE proposing changes so we can revert if the change hurts
        snapshot = self.skill_mgr.snapshot()

        # Propose changes before running (except baseline)
        change_label = ""
        if iteration > 1:
            change_label = self.propose_skill_changes(
                self.best_score,
                self.history[-1].error_breakdown if self.history else {},
                self.history,
            )

            # Reflection critic: reject edits that cheat (memorise ground truth)
            if self.enable_critic:
                from skills.clawpathy_autoresearch import critic as critic_mod

                diff = self.skill_mgr.diff_from_snapshot(snapshot)
                verdict = critic_mod.review_edit(
                    task_description=self.workspace.description or self.workspace.name,
                    ground_truth=self.workspace.ground_truth,
                    diff=diff,
                    proposer_label=change_label,
                )
                if not verdict.approved:
                    self.skill_mgr.restore(snapshot)
                    self._consecutive_non_improvements += 1
                    result = ExperimentResult(
                        experiment=iteration,
                        score=self.best_score,
                        kept=False,
                        label=f"critic-rejected: {change_label} ({verdict.reason})",
                        skill_diff=diff,
                        error_breakdown={},
                    )
                    self.history.append(result)
                    return result

        skill_path = self.workspace.skill_dir / "SKILL.md"
        output = self.run_agent_on_skill(skill_path, self.workspace.sources_dir)
        score_result = self.score_fn(output, self.workspace.ground_truth)
        if isinstance(score_result, tuple):
            score, error_breakdown = score_result
        else:
            score, error_breakdown = score_result, {}

        # Capture the diff BEFORE potentially restoring, so discarded edits are logged too
        diff = self.skill_mgr.diff_from_snapshot(snapshot)

        kept = score < self.best_score
        if kept:
            self.best_score = score
            self._consecutive_non_improvements = 0
            label = change_label or ("baseline" if iteration == 1 else "improvement")
        else:
            self.skill_mgr.restore(snapshot)
            self._consecutive_non_improvements += 1
            label = f"discarded: {change_label}" if change_label else "discarded"

        result = ExperimentResult(
            experiment=iteration,
            score=score,
            kept=kept,
            label=label,
            skill_diff=diff,
            error_breakdown=error_breakdown,
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

        # Final parallel LLM-judge pass (qualitative score + proposed additions)
        if not self._skip_final_judges:
            self._run_final_judges()

        return self.history

    _skip_final_judges = False  # override in tests

    def _run_final_judges(self) -> None:
        """Run the final LLM-judge pass. Safe to call only after loop completes."""
        try:
            from skills.clawpathy_autoresearch.final_judge import run_final_judges
        except ImportError:
            print("(final judges skipped: final_judge module unavailable)")
            return

        skill_path = self.workspace.skill_dir / "SKILL.md"
        skill_text = skill_path.read_text()
        output = self.run_agent_on_skill(skill_path, self.workspace.sources_dir)

        print("-" * 60)
        print("Running final LLM judges in parallel (quality + proposal)...")
        try:
            results = run_final_judges(
                workspace_name=self.workspace.name,
                skill_text=skill_text,
                output=output,
                ground_truth=self.workspace.ground_truth,
                output_dir=self.output_dir,
            )
            q = results["quality"].get("parsed") or {}
            p = results["proposal"].get("parsed") or {}
            q_err = results["quality"].get("error")
            p_err = results["proposal"].get("error")
            if q_err:
                print(f"Quality judge FAILED: {q_err}")
            else:
                print(f"Quality judge: {q.get('verdict', '?')} (score={q.get('score', '?')})")
            if p_err:
                print(f"Proposal judge FAILED: {p_err}")
            elif p.get("priority_change"):
                print(f"Proposal judge priority: {p['priority_change']}")
            print(f"Full judge output: {self.output_dir}/final_judges.json")
        except Exception as exc:
            print(f"(final judges failed: {exc})")


def run_demo(output_dir: Path) -> None:
    """Run a demo with synthetic data to show a descending error plot."""
    import random

    random.seed(42)
    history: list[ExperimentRecord] = []
    best = 0.98

    labels_kept = [
        "baseline",
        "added explicit input-format gotcha",
        "tightened trigger keywords",
        "split workflow step 3 into two substeps",
        "added output-schema example block",
        "reordered gotchas by severity",
        "clarified ambiguous parameter default",
        "added edge-case: empty input handling",
        "added cross-check step against reference",
        "pinned numerical tolerance in scorer",
        "added missing-data imputation note",
        "separated validation from transformation",
        "added provenance logging to workflow",
        "cached intermediate result between steps",
        "replaced prose instruction with table",
        "added explicit type annotations to output",
        "tuned retry count for flaky step",
        "added batching for large inputs",
        "deduplicated overlapping gotchas",
        "surfaced silent failure as loud error",
        "added domain-specific glossary",
        "linked to chaining partner skill",
        "added dry-run mode to workflow",
        "tightened output rounding precision",
        "added structured error-return contract",
        "added rate-limit handling note",
        "explicit ordering of multi-step output",
        "added regression guard for common mistake",
        "pre-flight input schema validation",
        "added progress reporting to long step",
    ]
    labels_disc = [
        "removed validation step",
        "aggressive gotcha pruning",
        "reordered workflow (lost chaining context)",
        "over-specified trigger conditions",
        "merged unrelated steps",
        "stripped safety checks",
        "added redundant API calls",
        "weakened scoring rubric",
        "skipped input normalisation",
        "removed deduplication pass",
        "collapsed gotchas into a single line",
        "hardcoded default that varies by input",
        "dropped edge-case handling",
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


def setup_workspace_interactive(output_dir: Path) -> Path:
    """Scaffold a fresh workspace by asking the user for task details.

    Creates task.json, ground_truth.json (empty stub), scorer.py (exact-match
    stub), skill/SKILL.md (task header only), and sources/ in output_dir.
    Returns the workspace path.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("No workspace found. Let's scaffold one.")
    name = input("Task name (short, e.g. 'IBD GWAS reproduction'): ").strip() or "Untitled task"
    description = input("One-paragraph description of what the agent must do: ").strip()
    max_iters_raw = input("Max iterations [30]: ").strip()
    max_iters = int(max_iters_raw) if max_iters_raw else 30
    early_stop_raw = input("Early-stop after N consecutive non-improvements [5]: ").strip()
    early_stop = int(early_stop_raw) if early_stop_raw else 5

    (output_dir / "task.json").write_text(
        json.dumps(
            {"name": name, "description": description,
             "max_iterations": max_iters, "early_stop_n": early_stop},
            indent=2,
        )
    )
    (output_dir / "ground_truth.json").write_text(
        json.dumps({"_note": "Fill in task-specific ground truth here."}, indent=2)
    )
    (output_dir / "scorer.py").write_text(
        'def score(skill_output: dict, ground_truth: dict):\n'
        '    """Return float (lower is better) or (float, breakdown_dict)."""\n'
        '    raise NotImplementedError("Implement scoring for this task.")\n'
    )
    skill_dir = output_dir / "skill"
    skill_dir.mkdir(exist_ok=True)
    (skill_dir / "SKILL.md").write_text(
        f"---\nname: {name.lower().replace(' ', '-')}\nversion: 0.1.0\n---\n\n"
        f"# {name}\n\n{description}\n\n## Workflow\n\n(to be proposed by autoresearch)\n"
    )
    (output_dir / "sources").mkdir(exist_ok=True)

    print(f"\nWorkspace scaffolded at: {output_dir}")
    print("Next: edit ground_truth.json and scorer.py, then rerun with --task", output_dir)
    return output_dir


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
    parser.add_argument("--setup", type=Path, help="Scaffold a new workspace interactively at this path")
    parser.add_argument("--validate", action="store_true", help="Validate workspace (load task + scorer) and exit")
    args = parser.parse_args()

    if args.demo:
        run_demo(args.output)
        return

    if args.setup:
        setup_workspace_interactive(args.setup)
        return

    if not args.task:
        parser.error("--task is required (or use --demo / --setup)")

    if args.validate:
        from skills.clawpathy_autoresearch.workspace import load_workspace
        ws = load_workspace(args.task)
        load_scorer(ws.scorer_path)
        print(f"OK: {ws.name} — {len(ws.ground_truth)} ground-truth keys, scorer loads.")
        return

    loop = AutoResearchLoop(
        workspace_dir=args.task,
        output_dir=args.output,
        max_iterations=args.iterations,
    )
    loop.run()


if __name__ == "__main__":
    main()
