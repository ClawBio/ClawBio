import json
import sys
from pathlib import Path

import pytest
from hypothesis import example, given, settings
from hypothesis import strategies as st

from clawbio.common import demo_runner


@pytest.fixture
def catalog(tmp_path, monkeypatch):
    """Point the runner at a throwaway catalog; returns a function that adds a demo to it."""
    path = tmp_path / "catalog.json"
    skills = []

    def add(name: str, demo_command: str) -> str:
        skills.append({"name": name, "demo_command": demo_command})
        path.write_text(json.dumps({"skills": skills}))
        return name

    monkeypatch.setattr(demo_runner, "CATALOG", path)
    return add


@pytest.fixture
def captured_runs(monkeypatch):
    calls = []

    def fake_run(args, **kwargs):
        calls.append((list(args), kwargs))
        return demo_runner.subprocess.CompletedProcess(args, 0, "", "")

    monkeypatch.setattr(demo_runner.subprocess, "run", fake_run)
    return calls


# --- running a catalog demo ---------------------------------------------------

def test_run_demo_does_not_interpret_demo_command_as_shell(catalog, tmp_path):
    marker = tmp_path / "marker"
    name = catalog("evil-skill", f"python -c pass $(touch {marker})")

    demo_runner.run_demo(name, tmp_path / "out", timeout=30)

    assert not marker.exists()


def test_run_demo_runs_argv_with_current_interpreter(catalog, captured_runs, tmp_path):
    name = catalog("pharmgx-reporter", "python clawbio.py run pharmgx --demo")

    result = demo_runner.run_demo(name, tmp_path, timeout=30)

    args, kwargs = captured_runs[0]
    assert args[:5] == [sys.executable, "clawbio.py", "run", "pharmgx", "--demo"]
    assert not kwargs.get("shell")
    assert result["passed"]


def test_run_demo_fails_closed_on_unparseable_command(catalog, tmp_path):
    name = catalog("broken-skill", "python 'unterminated")

    result = demo_runner.run_demo(name, tmp_path, timeout=30)

    assert not result["passed"]


def test_run_demo_fails_closed_on_a_skill_missing_from_the_catalog(catalog, tmp_path):
    catalog("some-other-skill", "python x.py --demo")

    result = demo_runner.run_demo("not-in-catalog", tmp_path, timeout=30)

    assert not result["passed"]
    assert "not-in-catalog" in result["stderr_tail"]


def test_run_demo_adds_output_when_the_command_has_none(catalog, captured_runs, tmp_path):
    name = catalog("vcf-annotator", "python skills/vcf-annotator/vcf_annotator.py --demo")

    demo_runner.run_demo(name, tmp_path, timeout=30)

    args, kwargs = captured_runs[0]
    assert args[-2:] == ["--output", str(tmp_path)]
    assert kwargs["cwd"] == str(demo_runner.PROJECT_ROOT)  # stays in the repo, so stray writes show in git status


def test_run_demo_replaces_an_existing_output(catalog, captured_runs, tmp_path):
    name = catalog("x", "python x.py --demo --output /tmp/x_demo")

    demo_runner.run_demo(name, tmp_path, timeout=30)

    args, _ = captured_runs[0]
    assert args.count("--output") == 1
    assert args[args.index("--output") + 1] == str(tmp_path)


def test_run_demo_keeps_the_stderr_of_a_failing_demo(catalog, tmp_path):
    name = catalog("crashes", "python -c \"import sys; sys.exit('demo exploded')\"")

    result = demo_runner.run_demo(name, tmp_path, timeout=30)

    assert not result["passed"]
    assert "demo exploded" in result["stderr_tail"]


_token = st.sampled_from(["--demo", "--output", "/tmp/a", "x.py", "--input", "f.vcf", "--output=/tmp/b"])


@settings(deadline=None)
@given(st.lists(_token, max_size=8))
def test_with_output_leaves_exactly_one_output_pointing_at_the_target(argv):
    target = Path("/out/skill")

    result = demo_runner.with_output(argv, target)

    outputs = [a for a in result if a == "--output" or a.startswith("--output=")]
    assert outputs == ["--output"]
    assert result[result.index("--output") + 1] == str(target)


# --- masking volatile fields --------------------------------------------------

OUT = Path("/tmp/pytest-of-ci/equity_demo0")
_text = st.lists(
    st.sampled_from([str(OUT), "/tmp", "equity_demo0", "2026-11-05T23:11:04.123456+00:00",
                     "2026-11-05 23:11 UTC", "76.2", ",", "\n", " ", "AFR"]),
    max_size=12,
).map("".join)


@given(_text)
@example(str(OUT) + "/tables/x.csv")
def test_mask_removes_the_output_dir(text):
    assert str(OUT) not in demo_runner.mask(text, OUT)


@given(_text)
@example("completed_at: 2026-11-05T23:11:04.123456+00:00")
def test_mask_is_idempotent(text):
    once = demo_runner.mask(text, OUT)
    assert demo_runner.mask(once, OUT) == once


@given(_text)
def test_two_runs_differing_only_in_time_and_dir_mask_equal(text):
    other = Path("/tmp/pytest-of-ci/equity_demo7")
    later = text.replace(str(OUT), str(other)).replace("23:11", "23:58")

    assert demo_runner.mask(text, OUT) == demo_runner.mask(later, other)


def test_mask_keeps_numbers_that_are_not_timestamps():
    assert demo_runner.mask("heim_score,76.2\nn_samples,50\n", OUT) == "heim_score,76.2\nn_samples,50\n"


# --- reading a demo's outputs for a snapshot ------------------------------------

_WRITER = '''
import sys
from pathlib import Path
out = Path(sys.argv[sys.argv.index("--output") + 1])
(out / "tables").mkdir(parents=True)
(out / "reproducibility").mkdir()
(out / "report.md").write_text(f"Date: 2026-11-05T23:11:04+00:00\\nOutput: {out}\\nRun: --output {out.name}\\n")
(out / "tables" / "scores.csv").write_text("pop,score\\nAFR,0.35\\n")
(out / "figures.png").write_bytes(b"\\x89PNG\\r\\n\\x1a\\n\\x00\\xff")
(out / "run.log").write_text("took 0.13s\\n")
(out / "reproducibility" / "commands.sh").write_text(f"python x.py --output {out}\\n")
(out / "big.tsv").write_text("x" * 60_000)
'''


@pytest.fixture
def writer_demo(catalog, tmp_path):
    script = tmp_path / "writer.py"
    script.write_text(_WRITER)
    return catalog("writer", f"python {script}")


def test_demo_outputs_returns_masked_text_outputs_by_relative_path(writer_demo, tmp_path):
    out = tmp_path / "out"

    outputs = demo_runner.demo_outputs(writer_demo, out)

    assert outputs == {
        "report.md": "Date: <TIMESTAMP>\nOutput: <OUTPUT_DIR>\nRun: --output demo_output\n",
        "tables/scores.csv": "pop,score\nAFR,0.35\n",
    }


def test_demo_outputs_raises_with_the_demo_stderr_when_it_fails(catalog, tmp_path):
    name = catalog("crashes", "python -c \"import sys; sys.exit('demo exploded')\"")

    with pytest.raises(RuntimeError, match="demo exploded"):
        demo_runner.demo_outputs(name, tmp_path / "out")


def test_demo_outputs_match_across_runs_into_different_dirs(writer_demo, tmp_path):
    # The writer echoes the output folder's bare name, as equity's report.md does.
    assert demo_runner.demo_outputs(writer_demo, tmp_path / "a") == demo_runner.demo_outputs(writer_demo, tmp_path / "b")
