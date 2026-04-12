import json
import subprocess
import sys
from pathlib import Path


def test_cli_dry_run_on_hello_task(tmp_path):
    ws = tmp_path / "ws"
    (ws / "skill").mkdir(parents=True)
    (ws / "skill" / "SKILL.md").write_text("")
    (ws / "task.json").write_text(json.dumps({
        "name": "hello", "max_iterations": 2, "early_stop_n": 2, "target_score": 0.0,
    }))
    (ws / "ground_truth.json").write_text(json.dumps({"target": "hello"}))
    (ws / "sandbox.yaml").write_text("allowed_tools: [Read]\n")
    (ws / "scorer.py").write_text(
        "import json\nfrom pathlib import Path\n"
        "def score(o, d):\n"
        "    gt = json.loads((Path(d)/'ground_truth.json').read_text())\n"
        "    return (0.0 if o.get('text')==gt['target'] else 1.0, {})\n"
    )
    result = subprocess.run(
        [sys.executable, "-m", "skills.clawpathy_autoresearch.autoresearch",
         "--workspace", str(ws), "--dry-run"],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    assert "best_score" in result.stdout
