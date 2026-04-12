import json
from pathlib import Path

from skills.clawpathy_autoresearch.dispatcher import CallableDispatcher, DispatchRequest
from skills.clawpathy_autoresearch.loop import run_loop


def _build_task(tmp_path: Path, target_text: str) -> Path:
    ws = tmp_path / "ws"
    (ws / "skill").mkdir(parents=True)
    (ws / "skill" / "SKILL.md").write_text("")
    (ws / "task.json").write_text(json.dumps({
        "name": "hello",
        "max_iterations": 3,
        "early_stop_n": 2,
        "target_score": 0.0,
    }))
    (ws / "ground_truth.json").write_text(json.dumps({"target": target_text}))
    (ws / "sandbox.yaml").write_text("allowed_tools: [Read]\n")
    (ws / "scorer.py").write_text(
        "import json\n"
        "from pathlib import Path\n"
        "def score(output, task_dir):\n"
        "    gt = json.loads((Path(task_dir)/'ground_truth.json').read_text())\n"
        "    return (0.0 if output.get('text')==gt['target'] else 1.0, {})\n"
    )
    return ws


def test_loop_converges_on_trivial_task(tmp_path):
    ws = _build_task(tmp_path, target_text="hello")

    def fake(req: DispatchRequest) -> str:
        if req.role == "proposer":
            return '```markdown\n# SKILL\nReturn {"text": "hello"}\n```'
        return '```json\n{"text": "hello"}\n```'

    d = CallableDispatcher(fake)
    result = run_loop(workspace_dir=ws, dispatcher=d)
    assert result["best_score"] == 0.0
    assert result["iterations_run"] >= 1
    history = [json.loads(l) for l in (ws / "history.jsonl").read_text().splitlines()]
    assert history[0]["kept"] is True
    assert (ws / "skill" / "SKILL.md").read_text().startswith("# SKILL")


def test_loop_reverts_when_score_worsens(tmp_path):
    ws = _build_task(tmp_path, target_text="hello")
    # Drop target_score so the loop doesn't stop on first success; cap iters
    task_json = ws / "task.json"
    data = json.loads(task_json.read_text())
    data.pop("target_score", None)
    data["max_iterations"] = 2
    task_json.write_text(json.dumps(data))
    proposals = iter([
        '```markdown\n# good\nreturn hello\n```',
        '```markdown\n# bad\nreturn wrong\n```',
    ])
    outputs = iter([
        '{"text": "hello"}',
        '{"text": "goodbye"}',
    ])

    def fake(req):
        return next(proposals) if req.role == "proposer" else next(outputs)

    d = CallableDispatcher(fake)
    result = run_loop(workspace_dir=ws, dispatcher=d)
    assert result["best_score"] == 0.0
    assert "good" in (ws / "skill" / "SKILL.md").read_text()
    history = [json.loads(l) for l in (ws / "history.jsonl").read_text().splitlines()]
    assert history[0]["kept"] is True
    assert history[1]["kept"] is False


def test_loop_stops_on_early_stop(tmp_path):
    ws = _build_task(tmp_path, target_text="hello")
    def fake(req):
        if req.role == "proposer":
            return '```markdown\nbad\n```'
        return '{"text": "wrong"}'
    d = CallableDispatcher(fake)
    result = run_loop(workspace_dir=ws, dispatcher=d)
    assert result["iterations_run"] <= 3
    assert result["best_score"] == 1.0
