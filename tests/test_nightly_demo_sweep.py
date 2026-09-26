import importlib.util
import sys
from pathlib import Path


def _load_sweep_module():
    module_path = Path(__file__).resolve().parents[1] / "scripts" / "nightly_demo_sweep.py"
    spec = importlib.util.spec_from_file_location("nightly_demo_sweep", module_path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_run_demo_does_not_interpret_demo_command_as_shell(tmp_path):
    sweep = _load_sweep_module()
    marker = tmp_path / "marker"
    skill = {"name": "evil-skill", "demo_command": f"python -c pass $(touch {marker})"}

    sweep.run_demo(skill, timeout=30, output_dir=None)

    assert not marker.exists()


def test_run_demo_runs_argv_with_current_interpreter(monkeypatch):
    sweep = _load_sweep_module()
    calls = []

    def fake_run(args, **kwargs):
        calls.append((args, kwargs))
        return sweep.subprocess.CompletedProcess(args, 0, "", "")

    monkeypatch.setattr(sweep.subprocess, "run", fake_run)
    skill = {"name": "pharmgx-reporter", "demo_command": "python clawbio.py run pharmgx --demo"}

    result = sweep.run_demo(skill, timeout=30, output_dir=None)

    args, kwargs = calls[0]
    assert args == [sys.executable, "clawbio.py", "run", "pharmgx", "--demo"]
    assert not kwargs.get("shell")
    assert result["passed"]


def test_run_demo_fails_closed_on_unparseable_command():
    sweep = _load_sweep_module()
    skill = {"name": "broken-skill", "demo_command": "python 'unterminated"}

    result = sweep.run_demo(skill, timeout=30, output_dir=None)

    assert not result["passed"]
