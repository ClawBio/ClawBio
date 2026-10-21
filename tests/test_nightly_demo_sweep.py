import importlib.util
import sys
from pathlib import Path

import pytest


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



def _write_bundle(out: Path, files=("commands.sh", "environment.yml", "checksums.sha256")):
    (out / "reproducibility").mkdir(parents=True)
    for name in files:
        (out / "reproducibility" / name).write_text("x\n")


def test_run_demo_records_which_reproducibility_files_the_demo_left_out(monkeypatch, tmp_path):
    sweep = _load_sweep_module()

    def fake_run(args, **kwargs):
        _write_bundle(Path(args[args.index("--output") + 1]), files=("commands.sh",))
        return sweep.subprocess.CompletedProcess(args, 0, "", "")

    monkeypatch.setattr(sweep.subprocess, "run", fake_run)

    result = sweep.run_demo({"name": "x", "demo_command": "python x.py --demo"}, timeout=30, output_dir=tmp_path)

    assert result["contract_gaps"] == {
        "repro_files": "reproducibility/ is missing environment.yml, checksums.sha256"}


def test_run_demo_checks_the_contract_before_removing_a_temp_output_dir(monkeypatch):
    sweep = _load_sweep_module()

    def fake_run(args, **kwargs):
        _write_bundle(Path(args[args.index("--output") + 1]))
        return sweep.subprocess.CompletedProcess(args, 0, "", "")

    monkeypatch.setattr(sweep.subprocess, "run", fake_run)

    result = sweep.run_demo({"name": "x", "demo_command": "python x.py --demo"}, timeout=30, output_dir=None)

    assert result["contract_gaps"] == {}


def _result(name, gaps, passed=True):
    return {"name": name, "passed": passed, "contract_gaps": gaps}


def test_a_gap_not_listed_for_the_skill_fails():
    sweep = _load_sweep_module()

    problems = sweep.contract_regressions(
        [_result("new-skill", {"repro_files": "reproducibility/ is missing checksums.sha256"})], known={})

    assert problems == ["new-skill: reproducibility/ is missing checksums.sha256"]


def test_a_listed_gap_the_skill_now_passes_must_leave_the_table():
    sweep = _load_sweep_module()

    problems = sweep.contract_regressions([_result("fixed-skill", {})], known={"fixed-skill": {"repro_files"}})

    assert problems == ["fixed-skill: now passes repro_files; remove it from CONTRACT_KNOWN_GAPS"]


def test_listed_gaps_and_failed_demos_are_not_reported():
    sweep = _load_sweep_module()
    results = [_result("known-gap", {"repro_files": "missing"}), _result("crashed", {}, passed=False),
               _result("good", {})]

    assert sweep.contract_regressions(results, known={"known-gap": {"repro_files"}}) == []



def test_gaps_for_checks_the_sweep_does_not_run_are_left_to_their_own_tests():
    sweep = _load_sweep_module()

    problems = sweep.contract_regressions([_result("skill", {})], known={"skill": {"demo_flag"}})

    assert problems == []

def test_the_sweep_exits_nonzero_on_a_contract_gap_even_when_every_demo_passed(monkeypatch, tmp_path):
    sweep = _load_sweep_module()
    monkeypatch.setattr(sweep, "load_demo_skills", lambda: [{"name": "new-skill", "demo_command": "python x.py --demo"}])
    monkeypatch.setattr(sweep, "run_demo", lambda skill, timeout, output_dir: {
        "name": "new-skill", "command": "", "exit_code": 0, "elapsed": 0.0, "stderr_tail": "", "passed": True,
        "contract_gaps": {"repro_files": "reproducibility/ is missing commands.sh"}})
    monkeypatch.delenv("GITHUB_STEP_SUMMARY", raising=False)
    monkeypatch.setattr(sys, "argv", ["nightly_demo_sweep.py", "--output", str(tmp_path)])

    with pytest.raises(SystemExit) as exit_info:
        sweep.main()

    assert exit_info.value.code == 1
    assert "new-skill: reproducibility/ is missing commands.sh" in (tmp_path / "sweep_summary.md").read_text()
