import importlib.util
from pathlib import Path

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
