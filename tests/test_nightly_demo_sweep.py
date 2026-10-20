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
    assert args[:5] == [sys.executable, "clawbio.py", "run", "pharmgx", "--demo"]
    assert not kwargs.get("shell")
    assert result["passed"]


def test_run_demo_fails_closed_on_unparseable_command():
    sweep = _load_sweep_module()
    skill = {"name": "broken-skill", "demo_command": "python 'unterminated"}

    result = sweep.run_demo(skill, timeout=30, output_dir=None)

    assert not result["passed"]


# --- every demo gets an explicit --output -----------------------------------

from hypothesis import given, settings
from hypothesis import strategies as st


def _capture_runs(monkeypatch, sweep):
    calls = []

    def fake_run(args, **kwargs):
        calls.append((list(args), kwargs))
        return sweep.subprocess.CompletedProcess(args, 0, "", "")

    monkeypatch.setattr(sweep.subprocess, "run", fake_run)
    return calls


def test_run_demo_adds_output_when_the_command_has_none(monkeypatch, tmp_path):
    sweep = _load_sweep_module()
    calls = _capture_runs(monkeypatch, sweep)

    sweep.run_demo({"name": "vcf-annotator", "demo_command": "python skills/vcf-annotator/vcf_annotator.py --demo"},
                   timeout=30, output_dir=tmp_path)

    args, kwargs = calls[0]
    assert args[-2:] == ["--output", str(tmp_path / "vcf-annotator")]
    assert kwargs["cwd"] == str(sweep.PROJECT_ROOT)  # stays in the repo, so stray writes show in git status


def test_run_demo_replaces_an_existing_output(monkeypatch, tmp_path):
    sweep = _load_sweep_module()
    calls = _capture_runs(monkeypatch, sweep)

    sweep.run_demo({"name": "x", "demo_command": "python x.py --demo --output /tmp/x_demo"},
                   timeout=30, output_dir=tmp_path)

    args, _ = calls[0]
    assert args.count("--output") == 1
    assert args[args.index("--output") + 1] == str(tmp_path / "x")


def test_run_demo_without_output_dir_uses_a_temp_dir_and_removes_it(monkeypatch):
    sweep = _load_sweep_module()
    seen = []

    def fake_run(args, **kwargs):
        out = Path(args[args.index("--output") + 1])
        assert out.is_dir()
        seen.append(out)
        return sweep.subprocess.CompletedProcess(args, 0, "", "")

    monkeypatch.setattr(sweep.subprocess, "run", fake_run)

    sweep.run_demo({"name": "x", "demo_command": "python x.py --demo"}, timeout=30, output_dir=None)

    assert not seen[0].exists()


_token = st.sampled_from(["--demo", "--output", "/tmp/a", "x.py", "--input", "f.vcf", "--output=/tmp/b"])


@settings(deadline=None)
@given(st.lists(_token, max_size=8))
def test_with_output_leaves_exactly_one_output_pointing_at_the_skill_dir(argv):
    sweep = _load_sweep_module()
    target = Path("/out/skill")

    result = sweep.with_output(argv, target)

    outputs = [a for a in result if a == "--output" or a.startswith("--output=")]
    assert outputs == ["--output"]
    assert result[result.index("--output") + 1] == str(target)
