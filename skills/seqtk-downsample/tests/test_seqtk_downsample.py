"""Tests for the seqtk-downsample skill (red/green TDD).

Unit tests build tiny FASTQs in tmp_path and need the ``seqtk`` binary
(on PATH, or named by SEQTK_BIN). They skip when it is missing.

The real-data test runs the skill on one 10x PBMC FASTQ and is skipped unless
SEQTK_TEST_FASTQ_DIR points at a directory holding such files.
"""

import gzip
import hashlib
import importlib.util
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
from hypothesis import example, given, settings
from hypothesis import strategies as st

SKILL_DIR = Path(__file__).resolve().parent.parent
SCRIPT = SKILL_DIR / "seqtk_downsample.py"

SEQTK = os.environ.get("SEQTK_BIN") or shutil.which("seqtk")
needs_seqtk = pytest.mark.skipif(
    SEQTK is None, reason="needs seqtk (conda install -c bioconda seqtk)"
)


def _load_module():
    spec = importlib.util.spec_from_file_location("seqtk_downsample", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["seqtk_downsample"] = mod
    spec.loader.exec_module(mod)
    return mod


def write_fastq(path: Path, n: int, mate: int = 1, order=None) -> Path:
    """Write n synthetic reads named read0..read{n-1}; gzip if path ends .gz."""
    ids = list(order) if order is not None else list(range(n))
    text = "".join(
        f"@read{i}/{mate} synthetic\nACGTACGTAC\n+\nIIIIIIIIII\n" for i in ids
    )
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "wt") as fh:
        fh.write(text)
    return path


def read_names(path: Path) -> list[str]:
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        lines = fh.read().splitlines()
    return [ln[1:].split()[0] for ln in lines[0::4]]


def run_cli(args, expect_code=0, cwd=None, env=None):
    full_env = {**os.environ, **(env or {})}
    result = subprocess.run(
        [sys.executable, str(SCRIPT), *map(str, args)],
        capture_output=True, text=True, cwd=cwd, env=full_env,
    )
    if expect_code is not None:
        assert result.returncode == expect_code, (
            f"Expected exit {expect_code}, got {result.returncode}\n"
            f"stdout: {result.stdout}\nstderr: {result.stderr}"
        )
    return result


def outputs_of(out: Path) -> list[Path]:
    return sorted((out / "downsampled").glob("*.fastq.gz"))


# ──────────────────────────────────────────────
# Demo mode
# ──────────────────────────────────────────────

@needs_seqtk
def test_demo_runs_and_writes_report(tmp_path):
    out = tmp_path / "demo"
    run_cli(["--demo", "--output", out])
    assert (out / "report.md").exists()
    assert (out / "result.json").exists()
    assert len(outputs_of(out)) == 2, "demo is paired-end: expect R1 and R2"


@needs_seqtk
def test_report_carries_the_disclaimer(tmp_path):
    out = tmp_path / "demo"
    run_cli(["--demo", "--output", out])
    assert (
        "ClawBio is a research and educational tool. It is not a medical "
        "device and does not provide clinical diagnoses."
    ) in (out / "report.md").read_text()


# ──────────────────────────────────────────────
# Single-end downsampling
# ──────────────────────────────────────────────

@needs_seqtk
def test_n_reads_gives_exactly_n_reads(tmp_path):
    src = write_fastq(tmp_path / "a.fastq", 200)
    out = tmp_path / "out"
    run_cli(["--input", src, "--output", out, "--n-reads", 25])
    (got,) = outputs_of(out)
    assert len(read_names(got)) == 25


@needs_seqtk
def test_output_reads_come_from_the_input(tmp_path):
    src = write_fastq(tmp_path / "a.fastq", 200)
    out = tmp_path / "out"
    run_cli(["--input", src, "--output", out, "--n-reads", 25])
    (got,) = outputs_of(out)
    assert set(read_names(got)) <= set(read_names(src))


@needs_seqtk
def test_asking_for_more_reads_than_exist_returns_all_and_says_so(tmp_path):
    src = write_fastq(tmp_path / "a.fastq", 30)
    out = tmp_path / "out"
    run_cli(["--input", src, "--output", out, "--n-reads", 1000])
    (got,) = outputs_of(out)
    assert len(read_names(got)) == 30
    data = json.loads((out / "result.json").read_text())
    assert data["fewer_than_requested"] is True
    assert "fewer" in (out / "report.md").read_text().lower()


@needs_seqtk
def test_fraction_mode_is_roughly_the_fraction(tmp_path):
    src = write_fastq(tmp_path / "a.fastq", 2000)
    out = tmp_path / "out"
    run_cli(["--input", src, "--output", out, "--fraction", 0.25])
    (got,) = outputs_of(out)
    assert 400 <= len(read_names(got)) <= 600


@needs_seqtk
def test_gzip_in_and_gzip_out(tmp_path):
    src = write_fastq(tmp_path / "a.fastq.gz", 100)
    out = tmp_path / "out"
    run_cli(["--input", src, "--output", out, "--n-reads", 10])
    (got,) = outputs_of(out)
    assert got.name.endswith(".fastq.gz")
    with gzip.open(got, "rt") as fh:  # raises if not real gzip
        assert fh.read().startswith("@read")


# ──────────────────────────────────────────────
# Seed
# ──────────────────────────────────────────────

@needs_seqtk
def test_default_seed_is_12345_and_is_recorded(tmp_path):
    src = write_fastq(tmp_path / "a.fastq", 100)
    out = tmp_path / "out"
    run_cli(["--input", src, "--output", out, "--n-reads", 10])
    assert json.loads((out / "result.json").read_text())["seed"] == 12345


@needs_seqtk
def test_default_n_reads_is_10000(tmp_path):
    src = write_fastq(tmp_path / "a.fastq", 50)
    out = tmp_path / "out"
    run_cli(["--input", src, "--output", out])
    data = json.loads((out / "result.json").read_text())
    assert data["mode"] == "n_reads"
    assert data["requested"] == 10000


@needs_seqtk
def test_same_seed_is_reproducible_and_different_seed_differs(tmp_path):
    src = write_fastq(tmp_path / "a.fastq", 500)
    names = {}
    for label, seed in (("a", 12345), ("b", 12345), ("c", 999)):
        out = tmp_path / label
        run_cli(["--input", src, "--output", out, "--n-reads", 20,
                 "--seed", seed])
        names[label] = read_names(outputs_of(out)[0])
    assert names["a"] == names["b"]
    assert names["a"] != names["c"]


@needs_seqtk
def test_same_seed_gives_byte_identical_gzip_output(tmp_path):
    """gzip headers carry a timestamp by default; checksums.sha256 must still
    match when the same command is re-run, so the header must be fixed."""
    import time

    src = write_fastq(tmp_path / "a.fastq", 200)
    digests = []
    for label in ("a", "b"):
        out = tmp_path / label
        run_cli(["--input", src, "--output", out, "--n-reads", 20])
        digests.append(hashlib.sha256(outputs_of(out)[0].read_bytes()).hexdigest())
        time.sleep(1.1)  # gzip mtime has 1-second resolution
    assert digests[0] == digests[1]


# ──────────────────────────────────────────────
# Paired-end
# ──────────────────────────────────────────────

@needs_seqtk
def test_paired_end_keeps_mates_in_sync(tmp_path):
    r1 = write_fastq(tmp_path / "s_R1.fastq.gz", 300, mate=1)
    r2 = write_fastq(tmp_path / "s_R2.fastq.gz", 300, mate=2)
    out = tmp_path / "out"
    run_cli(["--input", r1, "--input2", r2, "--output", out, "--n-reads", 40])
    o1, o2 = outputs_of(out)
    n1, n2 = read_names(o1), read_names(o2)
    assert len(n1) == len(n2) == 40
    strip = lambda names: [re.sub(r"/[12]$", "", n) for n in names]  # noqa: E731
    assert strip(n1) == strip(n2)
    assert json.loads((out / "result.json").read_text())["paired"] is True


@needs_seqtk
def test_paired_end_with_mismatched_read_order_fails_loudly(tmp_path):
    r1 = write_fastq(tmp_path / "s_R1.fastq", 300, mate=1)
    shuffled = list(range(300))[::-1]
    r2 = write_fastq(tmp_path / "s_R2.fastq", 300, mate=2, order=shuffled)
    out = tmp_path / "out"
    res = run_cli(["--input", r1, "--input2", r2, "--output", out,
                   "--n-reads", 40], expect_code=None)
    assert res.returncode != 0
    assert "order" in res.stderr.lower() or "mismatch" in res.stderr.lower()


# ──────────────────────────────────────────────
# Argument validation
# ──────────────────────────────────────────────

@pytest.mark.parametrize("bad", [
    ["--n-reads", 0],
    ["--n-reads", -5],
    ["--fraction", 0],
    # seqtk reads any amount >= 1 as a read COUNT, so --fraction 1.0 would
    # silently return 1 read instead of all of them. Must be refused.
    ["--fraction", 1.0],
    ["--fraction", 1.5],
    ["--fraction", "nan"],
    ["--fraction", "inf"],
    ["--n-reads", 10, "--fraction", 0.5],
])
def test_invalid_amounts_are_rejected(tmp_path, bad):
    src = write_fastq(tmp_path / "a.fastq", 50)
    res = run_cli(["--input", src, "--output", tmp_path / "out", *bad],
                  expect_code=None)
    assert res.returncode != 0
    assert not (tmp_path / "out" / "downsampled").exists()


def test_missing_input_is_an_error(tmp_path):
    res = run_cli(["--input", tmp_path / "nope.fastq",
                   "--output", tmp_path / "out"], expect_code=None)
    assert res.returncode != 0
    assert "nope.fastq" in res.stderr


def test_missing_seqtk_gives_a_clear_error(tmp_path):
    src = write_fastq(tmp_path / "a.fastq", 50)
    res = run_cli(["--input", src, "--output", tmp_path / "out"],
                  expect_code=None,
                  env={"SEQTK_BIN": str(tmp_path / "no-such-seqtk")})
    assert res.returncode != 0
    assert "seqtk" in res.stderr.lower()


# ──────────────────────────────────────────────
# Output hygiene
# ──────────────────────────────────────────────

@needs_seqtk
def test_existing_output_is_not_overwritten(tmp_path):
    src = write_fastq(tmp_path / "a.fastq", 100)
    out = tmp_path / "out"
    run_cli(["--input", src, "--output", out, "--n-reads", 10])
    before = outputs_of(out)[0].read_bytes()
    res = run_cli(["--input", src, "--output", out, "--n-reads", 20],
                  expect_code=None)
    assert res.returncode != 0
    assert outputs_of(out)[0].read_bytes() == before


@needs_seqtk
def test_writes_only_under_output_dir(tmp_path):
    src = write_fastq(tmp_path / "a.fastq", 100)
    work = tmp_path / "cwd"
    work.mkdir()
    run_cli(["--input", src, "--output", tmp_path / "out", "--n-reads", 10],
            cwd=work)
    assert list(work.iterdir()) == []


@needs_seqtk
def test_reproducibility_bundle_checksums_verify(tmp_path):
    src = write_fastq(tmp_path / "a.fastq", 100)
    out = tmp_path / "out"
    run_cli(["--input", src, "--output", out, "--n-reads", 10])
    repro = out / "reproducibility"
    for name in ("commands.sh", "environment.yml", "checksums.sha256"):
        assert (repro / name).exists(), name
    assert "seqtk" in (repro / "environment.yml").read_text()
    assert "--seed" in (repro / "commands.sh").read_text()
    for line in (repro / "checksums.sha256").read_text().splitlines():
        digest, label = line.split("  ", 1)
        actual = hashlib.sha256((out / label).read_bytes()).hexdigest()
        assert actual == digest, label


def _parse_output_contract(skill_md: Path) -> list[str]:
    text = skill_md.read_text()
    m = re.search(r"## Output Structure\s+```[^\n]*\n(.*?)```", text, re.S)
    if not m:
        return []
    files, stack = [], []
    for line in m.group(1).splitlines()[1:]:
        if not line.strip() or "(optional)" in line:
            continue
        # each tree level is a 4-char prefix ("├── " / "│   "); top level is 1
        depth = (len(line) - len(line.lstrip("│├└─ "))) // 4 - 1
        name = re.split(r"\s+#", line.lstrip("│├└─ "))[0].strip()
        stack = stack[:depth] + [name.rstrip("/")]
        if not name.endswith("/"):
            files.append("/".join(stack))
    return files


class TestOutputContract:
    """Every artifact promised in SKILL.md '## Output Structure' is produced."""

    @needs_seqtk
    def test_documented_outputs_are_produced(self, tmp_path):
        promised = _parse_output_contract(SKILL_DIR / "SKILL.md")
        assert promised, "SKILL.md has no parseable Output Structure"
        run_cli(["--demo", "--output", tmp_path])
        missing = [p for p in promised if not (tmp_path / p).exists()]
        assert not missing, f"promised but not produced: {missing}"


# ──────────────────────────────────────────────
# Properties
# ──────────────────────────────────────────────

@needs_seqtk
@settings(deadline=None, max_examples=25)
@given(
    total=st.integers(min_value=1, max_value=60),
    n=st.integers(min_value=1, max_value=80),
    seed=st.integers(min_value=0, max_value=2**31 - 1),
)
@example(total=1, n=80, seed=0)
@example(total=60, n=60, seed=12345)
def test_property_count_is_min_of_n_and_total_and_is_a_subset(
        tmp_path_factory, total, n, seed):
    tmp = tmp_path_factory.mktemp("prop")
    src = write_fastq(tmp / "a.fastq", total)
    out = tmp / "out"
    run_cli(["--input", src, "--output", out, "--n-reads", n,
             "--seed", seed])
    got = read_names(outputs_of(out)[0])
    assert len(got) == min(n, total)
    assert len(set(got)) == len(got), "a read was emitted twice"
    assert set(got) <= set(read_names(src))


@settings(deadline=None, max_examples=100)
@given(
    base=st.text(alphabet="ABCDEFGH0123456789:_-", min_size=1, max_size=20),
    suffix=st.sampled_from(["", "/1", "/2"]),
    comment=st.sampled_from(["", " 1:N:0:ACGT", " extra comment"]),
)
def test_property_mate_names_normalise_to_the_same_key(base, suffix, comment):
    mod = _load_module()
    assert mod.normalise_read_name(f"@{base}{suffix}{comment}") == base
    # idempotent
    once = mod.normalise_read_name(f"@{base}{suffix}{comment}")
    assert mod.normalise_read_name(once) == once


# ──────────────────────────────────────────────
# Real data (opt-in)
# ──────────────────────────────────────────────

REAL_DIR = os.environ.get("SEQTK_TEST_FASTQ_DIR")


@pytest.mark.slow
@pytest.mark.integration
@needs_seqtk
@pytest.mark.skipif(
    not REAL_DIR or not Path(REAL_DIR or ".").is_dir(),
    reason="set SEQTK_TEST_FASTQ_DIR to a directory of 10x PBMC FASTQs",
)
def test_real_pbmc_fastq_downsamples_to_10000_reads(tmp_path):
    src = sorted(Path(REAL_DIR).glob("read-RA_*lane-001*.fastq.gz"))[0]
    out = tmp_path / "out"
    run_cli(["--input", src, "--output", out, "--n-reads", 10000,
             "--seed", 12345])
    (got,) = outputs_of(out)
    assert len(read_names(got)) == 10000
    data = json.loads((out / "result.json").read_text())
    assert data["files"][0]["reads_out"] == 10000
