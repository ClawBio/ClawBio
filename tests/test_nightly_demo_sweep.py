import importlib.util
import sys
from pathlib import Path

import pytest

# Running a demo (argv handling, --output, failure modes) is tested in
# clawbio/common/tests/test_demo_runner.py; this file covers what the sweep adds.


def _load_sweep_module():
    module_path = Path(__file__).resolve().parents[1] / "scripts" / "nightly_demo_sweep.py"
    spec = importlib.util.spec_from_file_location("nightly_demo_sweep", module_path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def _capture_demo_runs(monkeypatch, sweep):
    seen = []

    def fake_run_demo(name, output_dir, *, timeout):
        assert Path(output_dir).is_dir()
        seen.append((name, Path(output_dir), timeout))
        return {"name": name, "passed": True}

    monkeypatch.setattr(sweep.demo_runner, "run_demo", fake_run_demo)
    return seen


def test_run_demo_gives_each_skill_its_own_dir_under_output(monkeypatch, tmp_path):
    sweep = _load_sweep_module()
    seen = _capture_demo_runs(monkeypatch, sweep)

    sweep.run_demo({"name": "vcf-annotator"}, timeout=30, output_dir=tmp_path)

    assert seen == [("vcf-annotator", tmp_path / "vcf-annotator", 30)]


def test_run_demo_without_output_dir_uses_a_temp_dir_and_removes_it(monkeypatch):
    sweep = _load_sweep_module()
    seen = _capture_demo_runs(monkeypatch, sweep)

    sweep.run_demo({"name": "x"}, timeout=30, output_dir=None)

    assert not seen[0][1].exists()


def _write_bundle(out: Path, files=("commands.sh", "environment.yml", "checksums.sha256")):
    (out / "reproducibility").mkdir(parents=True)
    for name in files:
        (out / "reproducibility" / name).write_text("x\n")


def test_run_demo_records_which_reproducibility_files_the_demo_left_out(monkeypatch, tmp_path):
    sweep = _load_sweep_module()

    def fake_run_demo(name, output_dir, *, timeout):
        _write_bundle(output_dir, files=("commands.sh",))
        return {"name": name, "passed": True}

    monkeypatch.setattr(sweep.demo_runner, "run_demo", fake_run_demo)

    result = sweep.run_demo({"name": "x"}, timeout=30, output_dir=tmp_path)

    assert result["contract_gaps"] == {
        "repro_files": "reproducibility/ is missing environment.yml, checksums.sha256"}


def test_run_demo_checks_the_contract_before_removing_a_temp_output_dir(monkeypatch):
    sweep = _load_sweep_module()

    def fake_run_demo(name, output_dir, *, timeout):
        _write_bundle(output_dir)
        return {"name": name, "passed": True}

    monkeypatch.setattr(sweep.demo_runner, "run_demo", fake_run_demo)

    result = sweep.run_demo({"name": "x"}, timeout=30, output_dir=None)

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
