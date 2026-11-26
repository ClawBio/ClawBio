import json
import os
import subprocess
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


# Points the runner at the catalog named by the first argument and prints the
# run_demo result for the skill it holds.
_RUN_DEMO_FROM_CATALOG = """
import json, sys
from pathlib import Path
from clawbio.common import demo_runner
demo_runner.CATALOG = Path(sys.argv[1])
print(json.dumps(demo_runner.run_demo("tcga-query", Path(sys.argv[2]), timeout=30)))
"""


def test_run_demo_reads_the_catalog_as_utf8_whatever_the_locale(tmp_path):
    # generate_catalog.py writes the catalog as UTF-8 and keeps non-ASCII text,
    # such as xena-tcga-gene-query's Chinese trigger keywords. 预后 is E5 90 8E
    # in UTF-8, and neither ASCII nor cp1252 can decode 0x90.
    catalog = tmp_path / "catalog.json"
    skill = {"name": "tcga-query", "trigger_keywords": ["生存分析", "预后"], "demo_command": "python -c pass"}
    catalog.write_text(json.dumps({"skills": [skill]}, ensure_ascii=False), encoding="utf-8")
    # With locale coercion and UTF-8 mode off, the C locale gives Python an
    # ASCII locale encoding on Linux and macOS. Windows ignores LC_ALL and uses
    # its ANSI code page, such as cp1252.
    env = {**os.environ, "LC_ALL": "C", "PYTHONCOERCECLOCALE": "0", "PYTHONUTF8": "0"}
    env.pop("PYTHONIOENCODING", None)

    proc = subprocess.run(
        [sys.executable, "-c", _RUN_DEMO_FROM_CATALOG, str(catalog), str(tmp_path / "out")],
        cwd=demo_runner.PROJECT_ROOT,
        env=env,
        capture_output=True,
        text=True,
        timeout=120,
    )

    assert proc.returncode == 0, proc.stderr
    result = json.loads(proc.stdout)
    assert result["passed"], result["stderr_tail"]


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


# --- rounding away floating-point noise ----------------------------------------
# A library upgrade can move a computed float in its last bit (scipy 1.18's
# moment code moved drug-repurposing-screen's bimodality_coef from
# 0.41691191789245413 to 0.4169119178924542). round_floats lets a snapshot pin
# the first N significant digits instead, so it fails on real changes only.

def test_round_floats_rounds_only_long_decimals():
    text = "BRD-0015,other,0.5,0.635056282655958,0.41691191789245413,1.0,10\n"
    assert demo_runner.round_floats(text, 12) == "BRD-0015,other,0.5,0.635056282656,0.416911917892,1.0,10\n"


def test_round_floats_absorbs_a_last_bit_change():
    before, after = "bc 0.41691191789245413", "bc 0.4169119178924542"
    assert demo_runner.round_floats(before, 12) == demo_runner.round_floats(after, 12)


def test_round_floats_still_sees_a_real_change():
    assert demo_runner.round_floats("bc 0.4169119178", 6) != demo_runner.round_floats("bc 0.4169919178", 6)


@pytest.mark.parametrize("text", [
    "numpy 2.5.4", "10.0.0.1", "rs3798220-C", "GCST90002370", "1.00e-45", "p=5.00e-22",
    "x1.234567890123456", "1.234567890123456y", "v1.2345678901234567.8",
])
def test_round_floats_leaves_versions_ids_and_short_numbers_alone(text):
    assert demo_runner.round_floats(text, 12) == text


@given(st.floats(allow_nan=False, allow_infinity=False, width=64))
@example(0.41691191789245413)
@example(-0.00045888077361219204)
@example(1e-300)
def test_property_round_floats_keeps_n_significant_digits(x):
    text = repr(x)
    out = demo_runner.round_floats(text, 12)
    assert demo_runner.round_floats(out, 12) == out  # idempotent
    assert float(out) == pytest.approx(x, rel=5e-12, abs=0.0) or float(out) == x


def test_demo_outputs_rounds_floats_only_when_asked(catalog, tmp_path):
    script = tmp_path / "floats.py"
    script.write_text(
        "import sys\nfrom pathlib import Path\n"
        "out = Path(sys.argv[sys.argv.index('--output') + 1]); out.mkdir(parents=True)\n"
        "(out / 'r.csv').write_text('a,0.41691191789245413\\n')\n"
    )
    name = catalog("floats", f"python {script}")

    assert demo_runner.demo_outputs(name, tmp_path / "a") == {"r.csv": "a,0.41691191789245413\n"}
    assert demo_runner.demo_outputs(name, tmp_path / "b", float_digits=12) == {"r.csv": "a,0.416911917892\n"}
