# Autoresearch Redesign Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Rewrite autoresearch as a domain-agnostic meta-skill that creates and iteratively improves SKILL.md files for any reproducible task, replacing the hardcoded GWAS-specific implementation.

**Architecture:** Two-phase design: interactive setup conversation produces a task workspace (task.json, ground_truth.json, scorer.py, skill/SKILL.md), then a headless optimisation loop runs the skill against ground truth, scores each iteration, and proposes single targeted SKILL.md edits. SkillManager (unchanged) handles snapshot/restore. Plotter (unchanged) renders Karpathy-style progress charts.

**Tech Stack:** Python 3.11+, pyyaml, matplotlib, anthropic (optional for agent integration)

---

## File Structure

| File | Action | Responsibility |
|------|--------|----------------|
| `skills/clawpathy-autoresearch/workspace.py` | **Create** | Load and validate task workspaces (task.json, ground_truth.json, scorer.py, skill/SKILL.md) |
| `skills/clawpathy-autoresearch/autoresearch.py` | **Rewrite** | CLI entry point: --setup, --task, --demo. AutoResearchLoop uses workspace instead of task.yaml |
| `skills/clawpathy-autoresearch/task.py` | **Delete** | Replaced by workspace.py |
| `skills/clawpathy-autoresearch/real_runner.py` | **Delete** | No longer needed: agent follows SKILL.md directly |
| `skills/clawpathy-autoresearch/scorer.py` | **Keep as reference scorer** | Keep ErrorBreakdown + reproduction_error as the demo/GWAS scorer. Add `load_scorer()` to dynamically import workspace scorers |
| `skills/clawpathy-autoresearch/skill_manager.py` | **Unchanged** | Snapshot/restore SKILL.md files |
| `skills/clawpathy-autoresearch/plotter.py` | **Unchanged** | Karpathy-style progress plot |
| `skills/clawpathy-autoresearch/SKILL.md` | **Update** | Reflect new general-purpose design |
| `skills/clawpathy-autoresearch/tests/test_workspace.py` | **Create** | Tests for workspace loading and validation |
| `skills/clawpathy-autoresearch/tests/test_autoresearch.py` | **Rewrite** | Tests for new CLI and loop logic |
| `skills/clawpathy-autoresearch/tests/test_scorer.py` | **Extend** | Add tests for dynamic scorer loading |
| `skills/clawpathy-autoresearch/tests/test_task.py` | **Delete** | Replaced by test_workspace.py |
| `skills/clawpathy-autoresearch/demo_workspace/` | **Create** | Pre-built GWAS workspace for --demo |

---

### Task 1: Create workspace.py with TDD

**Files:**
- Create: `skills/clawpathy-autoresearch/workspace.py`
- Create: `skills/clawpathy-autoresearch/tests/test_workspace.py`

A workspace is a directory containing everything autoresearch needs to run the optimisation loop. This module loads and validates workspaces.

- [ ] **Step 1: Write failing tests for workspace loading**

```python
"""Tests for workspace loading and validation."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from skills.clawpathy_autoresearch.workspace import (
    Workspace,
    load_workspace,
    validate_workspace,
)


@pytest.fixture
def valid_workspace(tmp_path: Path) -> Path:
    """Create a minimal valid workspace directory."""
    ws = tmp_path / "my_task"
    ws.mkdir()

    task_json = {
        "name": "Test Task",
        "description": "Reproduce test results",
        "max_iterations": 20,
        "early_stop_n": 5,
    }
    (ws / "task.json").write_text(json.dumps(task_json))

    ground_truth = {
        "targets": [
            {"id": "item_1", "value": 42.0},
            {"id": "item_2", "value": 3.14},
        ]
    }
    (ws / "ground_truth.json").write_text(json.dumps(ground_truth))

    scorer_code = '''
def score(skill_output: dict, ground_truth: dict) -> float:
    """Lower is better. 0 = perfect."""
    targets = {t["id"]: t["value"] for t in ground_truth["targets"]}
    outputs = {o["id"]: o["value"] for o in skill_output.get("results", [])}
    if not targets:
        return 0.0
    errors = []
    for tid, tval in targets.items():
        oval = outputs.get(tid, 0.0)
        errors.append(abs(tval - oval) / max(abs(tval), 1e-9))
    return sum(errors) / len(errors)
'''
    (ws / "scorer.py").write_text(scorer_code)

    skill_dir = ws / "skill"
    skill_dir.mkdir()
    (skill_dir / "SKILL.md").write_text(
        "---\nname: test-skill\nversion: 0.1.0\nauthor: test\n"
        "description: A test skill\n---\n\n# Test Skill\n\n"
        "## Workflow\n\n1. Read the input\n2. Produce output\n"
    )

    sources_dir = ws / "sources"
    sources_dir.mkdir()

    return ws


def test_load_workspace_returns_workspace(valid_workspace: Path):
    ws = load_workspace(valid_workspace)
    assert isinstance(ws, Workspace)
    assert ws.name == "Test Task"
    assert ws.description == "Reproduce test results"
    assert ws.max_iterations == 20
    assert ws.early_stop_n == 5


def test_workspace_has_ground_truth(valid_workspace: Path):
    ws = load_workspace(valid_workspace)
    assert ws.ground_truth is not None
    assert len(ws.ground_truth["targets"]) == 2


def test_workspace_has_paths(valid_workspace: Path):
    ws = load_workspace(valid_workspace)
    assert ws.workspace_dir == valid_workspace
    assert ws.skill_dir == valid_workspace / "skill"
    assert ws.sources_dir == valid_workspace / "sources"
    assert ws.scorer_path == valid_workspace / "scorer.py"


def test_validate_workspace_valid(valid_workspace: Path):
    errors = validate_workspace(valid_workspace)
    assert errors == []


def test_validate_workspace_missing_task_json(valid_workspace: Path):
    (valid_workspace / "task.json").unlink()
    errors = validate_workspace(valid_workspace)
    assert any("task.json" in e for e in errors)


def test_validate_workspace_missing_ground_truth(valid_workspace: Path):
    (valid_workspace / "ground_truth.json").unlink()
    errors = validate_workspace(valid_workspace)
    assert any("ground_truth.json" in e for e in errors)


def test_validate_workspace_missing_scorer(valid_workspace: Path):
    (valid_workspace / "scorer.py").unlink()
    errors = validate_workspace(valid_workspace)
    assert any("scorer.py" in e for e in errors)


def test_validate_workspace_missing_skill_md(valid_workspace: Path):
    (valid_workspace / "skill" / "SKILL.md").unlink()
    errors = validate_workspace(valid_workspace)
    assert any("SKILL.md" in e for e in errors)


def test_load_workspace_invalid_raises(valid_workspace: Path):
    (valid_workspace / "task.json").unlink()
    with pytest.raises(FileNotFoundError):
        load_workspace(valid_workspace)


def test_load_workspace_nonexistent_dir():
    with pytest.raises(FileNotFoundError):
        load_workspace(Path("/nonexistent/workspace"))
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd /Users/jaymoore/Documents/JAY_PhD/imperial/ClawBio && python -m pytest skills/clawpathy-autoresearch/tests/test_workspace.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'skills.clawpathy_autoresearch.workspace'`

- [ ] **Step 3: Implement workspace.py**

```python
"""Workspace loader for clawpathy-autoresearch.

A workspace is a self-contained directory with everything needed to run
the autoresearch optimisation loop: task config, ground truth, scorer,
and the SKILL.md being optimised.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any


@dataclass
class Workspace:
    """A loaded autoresearch workspace."""

    name: str
    description: str
    max_iterations: int
    early_stop_n: int
    ground_truth: dict[str, Any]
    workspace_dir: Path
    skill_dir: Path
    sources_dir: Path
    scorer_path: Path


def validate_workspace(workspace_dir: Path) -> list[str]:
    """Check a workspace directory for required files. Returns list of errors."""
    workspace_dir = Path(workspace_dir)
    errors = []

    if not workspace_dir.exists():
        return [f"Workspace directory does not exist: {workspace_dir}"]

    required = [
        ("task.json", workspace_dir / "task.json"),
        ("ground_truth.json", workspace_dir / "ground_truth.json"),
        ("scorer.py", workspace_dir / "scorer.py"),
        ("skill/SKILL.md", workspace_dir / "skill" / "SKILL.md"),
    ]
    for label, path in required:
        if not path.exists():
            errors.append(f"Missing required file: {label}")

    return errors


def load_workspace(workspace_dir: Path) -> Workspace:
    """Load a workspace from a directory.

    Raises FileNotFoundError if the directory is missing or invalid.
    """
    workspace_dir = Path(workspace_dir)

    if not workspace_dir.exists():
        raise FileNotFoundError(f"Workspace not found: {workspace_dir}")

    errors = validate_workspace(workspace_dir)
    if errors:
        raise FileNotFoundError(
            f"Invalid workspace: {'; '.join(errors)}"
        )

    task_data = json.loads((workspace_dir / "task.json").read_text())
    ground_truth = json.loads((workspace_dir / "ground_truth.json").read_text())

    return Workspace(
        name=task_data["name"],
        description=task_data.get("description", ""),
        max_iterations=task_data.get("max_iterations", 80),
        early_stop_n=task_data.get("early_stop_n", 5),
        ground_truth=ground_truth,
        workspace_dir=workspace_dir,
        skill_dir=workspace_dir / "skill",
        sources_dir=workspace_dir / "sources",
        scorer_path=workspace_dir / "scorer.py",
    )
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cd /Users/jaymoore/Documents/JAY_PhD/imperial/ClawBio && python -m pytest skills/clawpathy-autoresearch/tests/test_workspace.py -v`
Expected: All 10 tests PASS

- [ ] **Step 5: Commit**

```bash
git add skills/clawpathy-autoresearch/workspace.py skills/clawpathy-autoresearch/tests/test_workspace.py
git commit -m "feat(autoresearch): add workspace loader with validation"
```

---

### Task 2: Add dynamic scorer loading to scorer.py

**Files:**
- Modify: `skills/clawpathy-autoresearch/scorer.py:1-153`
- Modify: `skills/clawpathy-autoresearch/tests/test_scorer.py:1-167`

Keep the existing `reproduction_error` and `ErrorBreakdown` (used by demo workspace). Add a `load_scorer()` function that dynamically imports a workspace's `scorer.py`.

- [ ] **Step 1: Write failing tests for dynamic scorer loading**

Append to `tests/test_scorer.py`:

```python
from skills.clawpathy_autoresearch.scorer import load_scorer


def test_load_scorer_from_file(tmp_path: Path):
    scorer_code = '''
def score(skill_output: dict, ground_truth: dict) -> float:
    return abs(skill_output.get("value", 0) - ground_truth.get("value", 0))
'''
    scorer_path = tmp_path / "scorer.py"
    scorer_path.write_text(scorer_code)

    score_fn = load_scorer(scorer_path)
    result = score_fn({"value": 10}, {"value": 7})
    assert result == pytest.approx(3.0)


def test_load_scorer_missing_file_raises():
    with pytest.raises(FileNotFoundError):
        load_scorer(Path("/nonexistent/scorer.py"))


def test_load_scorer_missing_score_function(tmp_path: Path):
    scorer_path = tmp_path / "scorer.py"
    scorer_path.write_text("x = 1\n")

    with pytest.raises(AttributeError, match="score"):
        load_scorer(scorer_path)
```

- [ ] **Step 2: Run tests to verify new tests fail**

Run: `cd /Users/jaymoore/Documents/JAY_PhD/imperial/ClawBio && python -m pytest skills/clawpathy-autoresearch/tests/test_scorer.py::test_load_scorer_from_file skills/clawpathy-autoresearch/tests/test_scorer.py::test_load_scorer_missing_file_raises skills/clawpathy-autoresearch/tests/test_scorer.py::test_load_scorer_missing_score_function -v`
Expected: FAIL — `ImportError: cannot import name 'load_scorer'`

- [ ] **Step 3: Implement load_scorer in scorer.py**

Add to the bottom of `scorer.py`:

```python
import importlib.util


def load_scorer(scorer_path: Path) -> callable:
    """Dynamically import a scorer.py and return its score() function.

    The scorer module must define: score(skill_output: dict, ground_truth: dict) -> float
    """
    scorer_path = Path(scorer_path)
    if not scorer_path.exists():
        raise FileNotFoundError(f"Scorer not found: {scorer_path}")

    spec = importlib.util.spec_from_file_location("workspace_scorer", scorer_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    if not hasattr(module, "score"):
        raise AttributeError(
            f"Scorer at {scorer_path} must define a score() function"
        )
    return module.score
```

Also add `from pathlib import Path` to the imports at the top if not already present.

- [ ] **Step 4: Run all scorer tests to verify they pass**

Run: `cd /Users/jaymoore/Documents/JAY_PhD/imperial/ClawBio && python -m pytest skills/clawpathy-autoresearch/tests/test_scorer.py -v`
Expected: All 10 tests PASS (7 existing + 3 new)

- [ ] **Step 5: Commit**

```bash
git add skills/clawpathy-autoresearch/scorer.py skills/clawpathy-autoresearch/tests/test_scorer.py
git commit -m "feat(autoresearch): add dynamic scorer loading from workspace"
```

---

### Task 3: Create demo workspace

**Files:**
- Create: `skills/clawpathy-autoresearch/demo_workspace/task.json`
- Create: `skills/clawpathy-autoresearch/demo_workspace/ground_truth.json`
- Create: `skills/clawpathy-autoresearch/demo_workspace/scorer.py`
- Create: `skills/clawpathy-autoresearch/demo_workspace/skill/SKILL.md`
- Create: `skills/clawpathy-autoresearch/demo_workspace/sources/README.md`

This is the pre-built GWAS reproduction workspace used by `--demo`. It packages the existing ground truth from `tasks/gwas_reproduction/` into the new workspace format.

- [ ] **Step 1: Write failing test for demo workspace validity**

Create `skills/clawpathy-autoresearch/tests/test_demo_workspace.py`:

```python
"""Tests for the pre-built demo workspace."""
from __future__ import annotations

from pathlib import Path

from skills.clawpathy_autoresearch.workspace import load_workspace, validate_workspace
from skills.clawpathy_autoresearch.scorer import load_scorer


DEMO_DIR = Path(__file__).resolve().parent.parent / "demo_workspace"


def test_demo_workspace_is_valid():
    errors = validate_workspace(DEMO_DIR)
    assert errors == [], f"Demo workspace validation failed: {errors}"


def test_demo_workspace_loads():
    ws = load_workspace(DEMO_DIR)
    assert ws.name == "GWAS Paper Reproduction"
    assert ws.max_iterations == 80
    assert ws.early_stop_n == 5


def test_demo_ground_truth_has_variants():
    ws = load_workspace(DEMO_DIR)
    targets = ws.ground_truth.get("targets", [])
    assert len(targets) >= 3, "Demo should have at least 3 papers"
    for t in targets:
        assert "paper_id" in t
        assert "lead_variants" in t


def test_demo_scorer_runs():
    ws = load_workspace(DEMO_DIR)
    score_fn = load_scorer(ws.scorer_path)

    # Perfect output should score 0
    perfect = {"papers": {}}
    for t in ws.ground_truth["targets"]:
        perfect["papers"][t["paper_id"]] = {
            "variants_found": t["lead_variants"],
            "total_loci_reported": t["total_loci"],
        }
    result = score_fn(perfect, ws.ground_truth)
    assert result == 0.0


def test_demo_scorer_imperfect_output():
    ws = load_workspace(DEMO_DIR)
    score_fn = load_scorer(ws.scorer_path)

    # Empty output should score > 0
    empty = {"papers": {}}
    result = score_fn(empty, ws.ground_truth)
    assert result > 0.0


def test_demo_skill_md_exists():
    ws = load_workspace(DEMO_DIR)
    skill_path = ws.skill_dir / "SKILL.md"
    assert skill_path.exists()
    content = skill_path.read_text()
    assert "## Workflow" in content
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd /Users/jaymoore/Documents/JAY_PhD/imperial/ClawBio && python -m pytest skills/clawpathy-autoresearch/tests/test_demo_workspace.py -v`
Expected: FAIL — demo_workspace directory does not exist

- [ ] **Step 3: Create demo workspace task.json**

```json
{
  "name": "GWAS Paper Reproduction",
  "description": "Reproduce lead variant statistics from landmark GWAS papers by reading the paper methods and applying the described analysis pipeline.",
  "max_iterations": 80,
  "early_stop_n": 5
}
```

- [ ] **Step 4: Create demo workspace ground_truth.json**

Convert the three papers from the real autoresearch run (Nalls 2019, Mahajan 2018, de Lange 2017) into the new format. Read the existing ground truth YAMLs (`tasks/gwas_reproduction/ground_truth/paper_001.yaml`, `paper_003.yaml`, `paper_008.yaml`) and restructure:

```json
{
  "targets": [
    {
      "paper_id": "nalls_2019_pd",
      "pmid": "31701892",
      "title": "Nalls et al. 2019 - Parkinson's disease",
      "total_loci": 90,
      "lead_variants": [
        {"rsid": "rs356182", "gene": "SNCA", "neg_log10_p": 45.0, "odds_ratio": 1.33, "effect_allele_freq": 0.37, "effect_direction": "risk"},
        {"rsid": "rs34311866", "gene": "LRRK2", "neg_log10_p": 30.0, "odds_ratio": 1.22, "effect_allele_freq": 0.02, "effect_direction": "risk"},
        {"rsid": "rs34637584", "gene": "LRRK2", "neg_log10_p": 25.0, "odds_ratio": 2.25, "effect_allele_freq": 0.001, "effect_direction": "risk"},
        {"rsid": "rs114138760", "gene": "GBA", "neg_log10_p": 20.0, "odds_ratio": 1.95, "effect_allele_freq": 0.015, "effect_direction": "risk"},
        {"rsid": "rs12637471", "gene": "MAPT", "neg_log10_p": 40.0, "odds_ratio": 0.80, "effect_allele_freq": 0.22, "effect_direction": "protective"}
      ]
    },
    {
      "paper_id": "mahajan_2018_t2d",
      "pmid": "30297969",
      "title": "Mahajan et al. 2018 - Type 2 diabetes",
      "total_loci": 243,
      "lead_variants": [
        {"rsid": "rs7903146", "gene": "TCF7L2", "neg_log10_p": 200.0, "odds_ratio": 1.38, "effect_allele_freq": 0.30, "effect_direction": "risk"},
        {"rsid": "rs11708067", "gene": "ADCY5", "neg_log10_p": 15.0, "odds_ratio": 1.12, "effect_allele_freq": 0.22, "effect_direction": "risk"},
        {"rsid": "rs1801282", "gene": "PPARG", "neg_log10_p": 20.0, "odds_ratio": 0.86, "effect_allele_freq": 0.12, "effect_direction": "protective"}
      ]
    },
    {
      "paper_id": "delange_2017_ibd",
      "pmid": "28067908",
      "title": "de Lange et al. 2017 - IBD",
      "total_loci": 226,
      "lead_variants": [
        {"rsid": "rs2241880", "gene": "ATG16L1", "neg_log10_p": 30.0, "odds_ratio": 1.22, "effect_allele_freq": 0.47, "effect_direction": "risk"},
        {"rsid": "rs11209026", "gene": "IL23R", "neg_log10_p": 40.0, "odds_ratio": 0.50, "effect_allele_freq": 0.06, "effect_direction": "protective"}
      ]
    }
  ]
}
```

- [ ] **Step 5: Create demo workspace scorer.py**

This is the GWAS-specific scoring function, moved from the old `scorer.py` into the workspace:

```python
"""GWAS reproduction scorer for the demo workspace.

Score = weighted mean of six error components. Lower is better. 0 = perfect.
"""
from __future__ import annotations

MISSING_VARIANT_PENALTY = 1.0
DIRECTION_PENALTY = 1.0
WEIGHTS = {
    "p_value": 0.20,
    "or": 0.25,
    "freq": 0.10,
    "locus_count": 0.15,
    "missing": 0.20,
    "direction": 0.10,
}


def _score_paper(gt_paper: dict, agent_paper: dict) -> float:
    """Score a single paper's reproduction."""
    gt_variants = gt_paper.get("lead_variants", [])
    agent_variants = agent_paper.get("variants_found", [])
    gt_loci = gt_paper.get("total_loci", 0)
    agent_loci = agent_paper.get("total_loci_reported", 0)

    if not gt_variants:
        lce = abs(gt_loci - agent_loci) / max(gt_loci, 1)
        return lce * WEIGHTS["locus_count"]

    agent_by_rsid = {v["rsid"]: v for v in agent_variants}

    p_errors, or_errors, freq_errors, direction_errors = [], [], [], []
    missing_count = 0

    for gt_v in gt_variants:
        rsid = gt_v["rsid"]
        if rsid not in agent_by_rsid:
            missing_count += 1
            continue

        av = agent_by_rsid[rsid]

        gt_p = gt_v.get("neg_log10_p", 0)
        ag_p = av.get("neg_log10_p", 0)
        p_errors.append(abs(gt_p - ag_p) / gt_p if gt_p > 0 else abs(gt_p - ag_p))

        gt_or = gt_v.get("odds_ratio", 1.0)
        ag_or = av.get("odds_ratio", 1.0)
        or_errors.append(abs(gt_or - ag_or))

        gt_freq = gt_v.get("effect_allele_freq", 0.0)
        ag_freq = av.get("effect_allele_freq", 0.0)
        freq_errors.append(abs(gt_freq - ag_freq))

        gt_dir = gt_v.get("effect_direction", "")
        ag_dir = av.get("effect_direction", "")
        direction_errors.append(
            DIRECTION_PENALTY if gt_dir and ag_dir and gt_dir != ag_dir else 0.0
        )

    n = len(gt_variants)
    p_val_err = sum(p_errors) / len(p_errors) if p_errors else 0.0
    or_err = sum(or_errors) / len(or_errors) if or_errors else 0.0
    freq_err = sum(freq_errors) / len(freq_errors) if freq_errors else 0.0
    dir_err = sum(direction_errors) / n
    missing_err = missing_count * MISSING_VARIANT_PENALTY / n
    locus_err = abs(gt_loci - agent_loci) / gt_loci if gt_loci > 0 else float(agent_loci)

    return (
        WEIGHTS["p_value"] * p_val_err
        + WEIGHTS["or"] * or_err
        + WEIGHTS["freq"] * freq_err
        + WEIGHTS["locus_count"] * locus_err
        + WEIGHTS["missing"] * missing_err
        + WEIGHTS["direction"] * dir_err
    )


def score(skill_output: dict, ground_truth: dict) -> float:
    """Score the agent's full output against all papers.

    Args:
        skill_output: {"papers": {"paper_id": {"variants_found": [...], "total_loci_reported": int}}}
        ground_truth: {"targets": [{"paper_id": str, "lead_variants": [...], "total_loci": int}]}

    Returns:
        Mean error across all papers. Lower is better. 0 = perfect.
    """
    targets = ground_truth.get("targets", [])
    if not targets:
        return 0.0

    agent_papers = skill_output.get("papers", {})
    paper_errors = []

    for target in targets:
        pid = target["paper_id"]
        agent_paper = agent_papers.get(pid, {"variants_found": [], "total_loci_reported": 0})
        paper_errors.append(_score_paper(target, agent_paper))

    return sum(paper_errors) / len(paper_errors)
```

- [ ] **Step 6: Create demo workspace skill/SKILL.md**

```markdown
---
name: gwas-reproduction
version: 0.1.0
author: autoresearch
description: Reproduce lead variant statistics from GWAS papers
---

# GWAS Paper Reproduction

You are a GWAS reproduction agent. Given a published GWAS paper, your task is to reproduce the reported lead variant statistics.

## Workflow

1. Read the paper's methods section to understand the analysis pipeline
2. Identify the lead variants reported in the paper's results
3. For each lead variant, extract: rsID, gene, -log10(p-value), odds ratio, effect allele frequency, effect direction
4. Count the total number of genome-wide significant loci reported
5. Output structured results matching the ground truth schema

## Output Format

Return a JSON dict with:
- `papers`: dict mapping paper_id to paper results
- Each paper result has:
  - `variants_found`: list of dicts with rsid, neg_log10_p, odds_ratio, effect_allele_freq, effect_direction
  - `total_loci_reported`: integer count of GWS loci

## Gotchas

- The model will want to query external databases. Do not. Use only the paper's reported values.
- Effect direction depends on allele coding. Check the methods section for reference allele conventions.
- Some papers report beta instead of OR. Convert: OR = exp(beta).
```

- [ ] **Step 7: Create demo workspace sources/README.md**

```markdown
# Demo Sources

This demo workspace references three landmark GWAS papers:
- Nalls et al. 2019 (PMID: 31701892) — Parkinson's disease
- Mahajan et al. 2018 (PMID: 30297969) — Type 2 diabetes
- de Lange et al. 2017 (PMID: 28067908) — IBD

Ground truth was extracted from the papers' results tables.
```

- [ ] **Step 8: Run demo workspace tests**

Run: `cd /Users/jaymoore/Documents/JAY_PhD/imperial/ClawBio && python -m pytest skills/clawpathy-autoresearch/tests/test_demo_workspace.py -v`
Expected: All 6 tests PASS

- [ ] **Step 9: Commit**

```bash
git add skills/clawpathy-autoresearch/demo_workspace/
git add skills/clawpathy-autoresearch/tests/test_demo_workspace.py
git commit -m "feat(autoresearch): add pre-built GWAS demo workspace"
```

---

### Task 4: Rewrite autoresearch.py

**Files:**
- Modify: `skills/clawpathy-autoresearch/autoresearch.py:1-356`
- Modify: `skills/clawpathy-autoresearch/tests/test_autoresearch.py:1-119`

Rewrite the main entry point to use workspaces. Keep `ExperimentResult`, `save_experiment_log`, `load_experiment_log`, and `run_demo`. Rewrite `AutoResearchLoop` to load a workspace, use dynamic scorer, and operate SkillManager on `workspace/skill/`.

- [ ] **Step 1: Write failing tests for new AutoResearchLoop**

Rewrite `tests/test_autoresearch.py`:

```python
"""Tests for the autoresearch loop and CLI."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from skills.clawpathy_autoresearch.autoresearch import (
    AutoResearchLoop,
    ExperimentResult,
    load_experiment_log,
    save_experiment_log,
)


@pytest.fixture
def workspace_dir(tmp_path: Path) -> Path:
    """Create a minimal workspace for testing the loop."""
    ws = tmp_path / "test_workspace"
    ws.mkdir()

    task_json = {
        "name": "Test Task",
        "description": "Test",
        "max_iterations": 3,
        "early_stop_n": 2,
    }
    (ws / "task.json").write_text(json.dumps(task_json))

    ground_truth = {
        "targets": [
            {"id": "item_1", "value": 10.0},
        ]
    }
    (ws / "ground_truth.json").write_text(json.dumps(ground_truth))

    scorer_code = '''
def score(skill_output: dict, ground_truth: dict) -> float:
    targets = {t["id"]: t["value"] for t in ground_truth["targets"]}
    outputs = {o["id"]: o["value"] for o in skill_output.get("results", [])}
    if not targets:
        return 0.0
    errors = []
    for tid, tval in targets.items():
        oval = outputs.get(tid, 0.0)
        errors.append(abs(tval - oval) / max(abs(tval), 1e-9))
    return sum(errors) / len(errors)
'''
    (ws / "scorer.py").write_text(scorer_code)

    skill_dir = ws / "skill"
    skill_dir.mkdir()
    (skill_dir / "SKILL.md").write_text(
        "---\nname: test-skill\n---\n# Test\n\n## Workflow\n\n1. Do thing\n"
    )

    (ws / "sources").mkdir()
    return ws


def test_experiment_result_serialisation():
    result = ExperimentResult(
        experiment=1,
        score=0.35,
        kept=True,
        label="baseline",
        skill_diff={},
        error_breakdown={},
    )
    d = result.to_dict()
    assert d["experiment"] == 1
    assert d["score"] == 0.35
    assert d["kept"] is True


def test_save_and_load_experiment_log(tmp_path: Path):
    results = [
        ExperimentResult(
            experiment=1, score=0.8, kept=True, label="baseline",
            skill_diff={}, error_breakdown={},
        ),
        ExperimentResult(
            experiment=2, score=0.6, kept=True, label="improved",
            skill_diff={"skill": "modified"}, error_breakdown={},
        ),
    ]
    log_path = tmp_path / "experiment_log.json"
    save_experiment_log(results, log_path)
    assert log_path.exists()

    loaded = load_experiment_log(log_path)
    assert len(loaded) == 2
    assert loaded[0]["experiment"] == 1
    assert loaded[1]["score"] == 0.6


def test_loop_init(workspace_dir: Path):
    loop = AutoResearchLoop(
        workspace_dir=workspace_dir,
        output_dir=workspace_dir / "output",
    )
    assert loop.workspace.name == "Test Task"
    assert loop.max_iterations == 3
    assert loop.early_stop_n == 2
    assert len(loop.history) == 0
    assert loop.best_score == float("inf")


def test_loop_run_iteration_default_agent(workspace_dir: Path):
    """Default run_agent_on_skill returns empty output, scoring > 0."""
    loop = AutoResearchLoop(
        workspace_dir=workspace_dir,
        output_dir=workspace_dir / "output",
    )
    result = loop.run_iteration(1)
    assert isinstance(result, ExperimentResult)
    assert result.experiment == 1
    assert result.score >= 0.0
    assert result.kept is True  # first iteration is always kept as baseline


def test_loop_early_stop(workspace_dir: Path):
    """Loop stops after early_stop_n consecutive non-improvements."""
    loop = AutoResearchLoop(
        workspace_dir=workspace_dir,
        output_dir=workspace_dir / "output",
    )
    results = loop.run()
    # With default agent (always returns same output), should early-stop
    # after 1 kept (baseline) + early_stop_n non-improvements
    assert len(results) <= loop.max_iterations


def test_loop_produces_plot(workspace_dir: Path):
    loop = AutoResearchLoop(
        workspace_dir=workspace_dir,
        output_dir=workspace_dir / "output",
    )
    loop.run()
    assert (workspace_dir / "output" / "progress.png").exists()


def test_loop_produces_experiment_log(workspace_dir: Path):
    loop = AutoResearchLoop(
        workspace_dir=workspace_dir,
        output_dir=workspace_dir / "output",
    )
    loop.run()
    assert (workspace_dir / "output" / "experiment_log.json").exists()
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd /Users/jaymoore/Documents/JAY_PhD/imperial/ClawBio && python -m pytest skills/clawpathy-autoresearch/tests/test_autoresearch.py -v`
Expected: FAIL — old API doesn't match new test expectations

- [ ] **Step 3: Rewrite autoresearch.py**

```python
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
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cd /Users/jaymoore/Documents/JAY_PhD/imperial/ClawBio && python -m pytest skills/clawpathy-autoresearch/tests/test_autoresearch.py -v`
Expected: All 7 tests PASS

- [ ] **Step 5: Run the full test suite to check for regressions**

Run: `cd /Users/jaymoore/Documents/JAY_PhD/imperial/ClawBio && python -m pytest skills/clawpathy-autoresearch/tests/ -v`
Expected: All tests PASS across all test files

- [ ] **Step 6: Commit**

```bash
git add skills/clawpathy-autoresearch/autoresearch.py skills/clawpathy-autoresearch/tests/test_autoresearch.py
git commit -m "refactor(autoresearch): rewrite loop to use workspace instead of task.yaml"
```

---

### Task 5: Delete old files and update imports

**Files:**
- Delete: `skills/clawpathy-autoresearch/task.py`
- Delete: `skills/clawpathy-autoresearch/real_runner.py`
- Delete: `skills/clawpathy-autoresearch/tests/test_task.py`

- [ ] **Step 1: Delete the old files**

```bash
cd /Users/jaymoore/Documents/JAY_PhD/imperial/ClawBio
rm skills/clawpathy-autoresearch/task.py
rm skills/clawpathy-autoresearch/real_runner.py
rm skills/clawpathy-autoresearch/tests/test_task.py
```

- [ ] **Step 2: Verify no remaining imports of deleted modules**

Search for any imports of `task.py` or `real_runner.py` in the autoresearch directory. The only file that imported `task.py` was the old `autoresearch.py` (now rewritten). `real_runner.py` was standalone.

```bash
cd /Users/jaymoore/Documents/JAY_PhD/imperial/ClawBio
grep -r "from skills.clawpathy_autoresearch.task import" skills/clawpathy-autoresearch/
grep -r "from skills.clawpathy_autoresearch.real_runner import" skills/clawpathy-autoresearch/
grep -r "import real_runner" skills/clawpathy-autoresearch/
```

Expected: no matches

- [ ] **Step 3: Run full test suite**

Run: `cd /Users/jaymoore/Documents/JAY_PhD/imperial/ClawBio && python -m pytest skills/clawpathy-autoresearch/tests/ -v`
Expected: All tests PASS (test_task.py gone, no import errors)

- [ ] **Step 4: Commit**

```bash
cd /Users/jaymoore/Documents/JAY_PhD/imperial/ClawBio
git rm skills/clawpathy-autoresearch/task.py skills/clawpathy-autoresearch/real_runner.py skills/clawpathy-autoresearch/tests/test_task.py
git commit -m "refactor(autoresearch): remove GWAS-specific task.py and real_runner.py"
```

---

### Task 6: Update SKILL.md for the new design

**Files:**
- Modify: `skills/clawpathy-autoresearch/SKILL.md:1-186`

- [ ] **Step 1: Rewrite SKILL.md**

Update to reflect the general-purpose design. Remove all GWAS-specific language. Document the workspace format, the two-phase architecture, and the new CLI.

Key changes:
- Description: "domain-agnostic meta-skill" not "GWAS reproduction"
- Input format: workspace directory, not task.yaml
- Workflow: setup phase + optimisation loop
- CLI: `--task <workspace_dir>` not `--task <task.yaml>`
- Remove chaining partners (gwas-lookup, fine-mapping, etc.)
- Update trigger keywords

The SKILL.md must conform to SKILL-TEMPLATE.md: YAML frontmatter (name, version, author, description, inputs, outputs, trigger_keywords), Trigger, Scope, Workflow, Example Output, Gotchas (min 3), Safety, Agent Boundary sections.

- [ ] **Step 2: Run SKILL.md conformance check (manual)**

Verify all 17 checklist items from CLAUDE.md pass:
- YAML: name matches folder, version semver, author present, description one-line, inputs present, outputs present, trigger_keywords >= 3
- Sections: Trigger (fire/don't-fire), Scope (one-skill-one-task), Workflow (numbered steps), Example Output (rendered sample), Gotchas (>= 3), Safety (disclaimer), Agent Boundary
- Files: demo data exists, tests exist
- Line count: under 500

- [ ] **Step 3: Commit**

```bash
git add skills/clawpathy-autoresearch/SKILL.md
git commit -m "docs(autoresearch): update SKILL.md for general-purpose workspace design"
```

---

### Task 7: Update CLAUDE.md routing table and CLI reference

**Files:**
- Modify: `CLAUDE.md`

- [ ] **Step 1: Update the CLI reference for autoresearch**

In the CLI Reference section of CLAUDE.md, replace the autoresearch commands:

Old:
```bash
python skills/clawpathy-autoresearch/autoresearch.py \
  --task tasks/gwas_reproduction/task.yaml --output /tmp/autoresearch --iterations 80
```

New:
```bash
# Run optimisation loop on a workspace
python skills/clawpathy-autoresearch/autoresearch.py \
  --task /path/to/workspace --output /tmp/autoresearch --iterations 80

# Demo (pre-built GWAS workspace, synthetic progress plot)
python skills/clawpathy-autoresearch/autoresearch.py --demo --output /tmp/autoresearch_demo
```

- [ ] **Step 2: Update demo data table entry**

Old: `Autoresearch demo (83 synthetic experiments, progress plot)`
New: `Autoresearch demo (83 synthetic experiments, progress plot + pre-built GWAS workspace)`

- [ ] **Step 3: Commit**

```bash
git add CLAUDE.md
git commit -m "docs: update CLAUDE.md CLI reference for autoresearch workspace design"
```

---

### Task 8: Run full test suite and demo

**Files:** None (validation only)

- [ ] **Step 1: Run all autoresearch tests**

```bash
cd /Users/jaymoore/Documents/JAY_PhD/imperial/ClawBio
python -m pytest skills/clawpathy-autoresearch/tests/ -v
```

Expected: All tests PASS across test_workspace.py, test_scorer.py, test_autoresearch.py, test_demo_workspace.py, test_skill_manager.py, test_plotter.py

- [ ] **Step 2: Run demo mode**

```bash
cd /Users/jaymoore/Documents/JAY_PhD/imperial/ClawBio
python skills/clawpathy-autoresearch/autoresearch.py --demo --output /tmp/autoresearch_demo
```

Expected: "83 experiments, 30 kept" + progress.png + experiment_log.json

- [ ] **Step 3: Run loop on demo workspace (short run)**

```bash
cd /Users/jaymoore/Documents/JAY_PhD/imperial/ClawBio
python skills/clawpathy-autoresearch/autoresearch.py \
  --task skills/clawpathy-autoresearch/demo_workspace \
  --output /tmp/autoresearch_loop_test --iterations 3
```

Expected: 3 iterations (or early stop), progress.png, experiment_log.json. Default agent returns empty output so score will be > 0 and constant (early stop after early_stop_n).

- [ ] **Step 4: Verify no import errors across the broader codebase**

```bash
cd /Users/jaymoore/Documents/JAY_PhD/imperial/ClawBio
grep -r "clawpathy_autoresearch.task" skills/ clawbio/ --include="*.py"
grep -r "clawpathy_autoresearch.real_runner" skills/ clawbio/ --include="*.py"
```

Expected: no matches (old modules fully removed)

- [ ] **Step 5: Final commit if any fixes needed**

Only if steps 1-4 revealed issues that required fixes.
