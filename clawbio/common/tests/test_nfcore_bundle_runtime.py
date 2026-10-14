"""Relocation mechanics and emitted single-file standard-library bundles."""

import csv
import hashlib
import importlib
import subprocess
import sys
from pathlib import Path

import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st


def test_csv_preview_backup_and_prefix_boundary(tmp_path):
    api = importlib.import_module("clawbio.common.nfcore.bundle_runtime")
    path = tmp_path / "samplesheet.csv"
    original = "sample,fastq_1\nA,/data/a.fq\nB,/data2/b.fq\nC,s3://bucket/c.fq\n"
    path.write_text(original)
    changes = api.remap_csv(path, "/data", "/moved", columns=("fastq_1",), dry_run=True)
    assert changes == [("fastq_1", "/data/a.fq", "/moved/a.fq")]
    assert path.read_text() == original
    assert not (tmp_path / "samplesheet.csv.bak").exists()
    api.remap_csv(path, "/data", "/moved", columns=("fastq_1",), dry_run=False)
    assert (tmp_path / "samplesheet.csv.bak").read_text() == original
    assert path.read_text() == original.replace("/data/a.fq", "/moved/a.fq")


@settings(deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture])
@given(st.text(alphabet="abc é-_012", min_size=1, max_size=20))
def test_remap_round_trip(tmp_path, suffix):
    api = importlib.import_module("clawbio.common.nfcore.bundle_runtime")
    path = tmp_path / "input.csv"
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["fastq_1"])
        writer.writerow(["/old/" + suffix])
    original = path.read_bytes()
    api.remap_csv(path, "/old", "/new", columns=["fastq_1"], dry_run=False)
    api.remap_csv(path, "/new", "/old", columns=["fastq_1"], dry_run=False)
    assert path.read_bytes() == original


def test_all_actions_run_and_first_failure_wins():
    api = importlib.import_module("clawbio.common.nfcore.bundle_runtime")
    seen = []

    def action(code):
        def run():
            seen.append(code)
            return code

        return run

    assert api.run_actions([action(0), action(2), action(1), action(0)]) == 2
    assert seen == [0, 2, 1, 0]


@pytest.mark.parametrize("pipeline", ["rnaseq", "sarek", "scrnaseq"])
def test_emitted_remapper_runs_without_site_packages(pipeline, tmp_path):
    api = importlib.import_module("clawbio.common.nfcore.bundle_emit")
    root = Path(__file__).resolve().parents[3]
    adapter = root / "skills" / f"nfcore-{pipeline}-wrapper/remap_paths.py"
    dest = tmp_path / "remap_paths.py"
    api.write_standalone_remapper(dest, adapter_path=adapter)
    result = subprocess.run(
        [sys.executable, "-I", "-S", str(dest), "--help"],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    assert "--repair-bundle" in result.stdout
    path = tmp_path / "samplesheet.valid.csv"
    path.write_text("sample,fastq_1\nA,/old/a.fq.gz\n")
    result = subprocess.run(
        [
            sys.executable,
            "-I",
            "-S",
            str(dest),
            "--old",
            "/old",
            "--new",
            "/moved",
            "--dry-run",
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    assert path.read_text() == "sample,fastq_1\nA,/old/a.fq.gz\n"


_BLOCK = "# >>> SHARED NFCORE RUNTIME >>>\npass\n# <<< SHARED NFCORE RUNTIME <<<\n"


def test_emitter_rejects_missing_and_duplicate_markers():
    api = importlib.import_module("clawbio.common.nfcore.bundle_emit")
    with pytest.raises(ValueError):
        api.render_standalone_remapper("print(1)")
    with pytest.raises(ValueError):
        api.render_standalone_remapper(_BLOCK + _BLOCK)


@pytest.mark.parametrize(
    "adapter",
    [
        # After the block: the runtime's remap_csv would call this, not its own helper.
        _BLOCK + "def prefix_matches(path, prefix):\n    return False\n",
        _BLOCK + "escape_double_quoted = str\n",
        # Before the block: the pasted runtime would silently replace the adapter's.
        "def repair_bundle():\n    return 0\n" + _BLOCK,
    ],
)
def test_emitter_refuses_names_the_pasted_runtime_would_shadow(adapter):
    """The runtime is pasted into the adapter's own namespace, unlike the in-repo
    import. A colliding top-level name would make the shipped script behave
    differently from the tested one, so it must fail at generation time."""
    api = importlib.import_module("clawbio.common.nfcore.bundle_emit")
    with pytest.raises(ValueError, match="shadow"):
        api.render_standalone_remapper(adapter)


def test_emitter_allows_identical_stdlib_imports_and_an_adapter_entry_point():
    api = importlib.import_module("clawbio.common.nfcore.bundle_emit")
    adapter = (
        "import csv\nfrom pathlib import Path\n" + _BLOCK
        + "import re\n\ndef remap_csv(path):\n    return _nfcore_remap_csv(path)\n"
    )
    assert "def remap_csv(path):" in api.render_standalone_remapper(adapter)


@pytest.mark.parametrize("pipeline", ["rnaseq", "sarek", "scrnaseq"])
def test_standalone_remapper_mutates_verifies_and_repairs_bundle(pipeline, tmp_path):
    api = importlib.import_module("clawbio.common.nfcore.bundle_emit")
    root = Path(__file__).resolve().parents[3]
    adapter = root / "skills" / f"nfcore-{pipeline}-wrapper/remap_paths.py"
    bundle = tmp_path / "relocated é bundle" / "reproducibility"
    dest = bundle / "remap_paths.py"
    api.write_standalone_remapper(dest, adapter_path=adapter)
    reads = tmp_path / "moved é reads"
    reads.mkdir()
    (reads / "a.fq.gz").write_bytes(b"public test read")
    samplesheet = bundle / "samplesheet.valid.csv"
    original = b"sample,fastq_1\nA,/old/a.fq.gz\n"
    samplesheet.write_bytes(original)
    (bundle / "commands.sh").write_text("#!/bin/bash\ntrue\n")
    (bundle / "params.yaml").write_text(
        "input: reproducibility/samplesheet.valid.csv\n"
    )
    (bundle.parent / "report.md").write_text("Synthetic test report\n")

    def run(*args):
        result = subprocess.run(
            [sys.executable, "-I", "-S", str(dest), *args],
            cwd=tmp_path,
            capture_output=True,
            text=True,
            check=False,
        )
        assert result.returncode == 0, result.stdout + result.stderr

    run("--old", "/old", "--new", str(reads))
    backup = (
        samplesheet.with_suffix(".bak")
        if pipeline == "scrnaseq"
        else samplesheet.with_name(samplesheet.name + ".bak")
    )
    assert backup.read_bytes() == original
    assert list(csv.DictReader(samplesheet.read_text().splitlines())) == [
        {"sample": "A", "fastq_1": str(reads / "a.fq.gz")}
    ]
    run("--repair-bundle", "--verify")
    required = ("manifest.json", "environment.yml", "checksums.sha256")
    repaired = {name: (bundle / name).read_bytes() for name in required}
    assert b'"regenerated_post_hoc": true' in repaired["manifest.json"]
    for line in repaired["checksums.sha256"].decode().splitlines():
        digest, relative = line.split("  ", 1)
        assert (
            hashlib.sha256((bundle.parent / relative).read_bytes()).hexdigest()
            == digest
        )
    assert repaired["checksums.sha256"]
    run("--repair-bundle", "--verify")
    assert {name: (bundle / name).read_bytes() for name in required} == repaired


def test_repair_bundle_recomputes_checksum_after_new_metadata(tmp_path):
    api = importlib.import_module("clawbio.common.nfcore.bundle_runtime")
    (tmp_path / "commands.sh").write_text("#!/bin/bash\ntrue\n")
    (tmp_path / "checksums.sha256").write_text("old\n")

    def manifest(directory, commands):
        (directory / "manifest.json").write_text('{"regenerated_post_hoc":true}')

    def environment(directory):
        (directory / "environment.yml").write_text("name: recipe\n")

    def checksums(directory):
        assert (directory / "manifest.json").exists()
        assert (directory / "environment.yml").exists()
        (directory / "checksums.sha256").write_text("updated\n")

    assert (
        api.repair_bundle(
            tmp_path,
            required_files=("manifest.json", "environment.yml", "checksums.sha256"),
            commands_finder=lambda **kw: tmp_path / "commands.sh",
            regenerate_manifest=manifest,
            regenerate_environment=environment,
            regenerate_checksums=checksums,
        )
        == 0
    )
    assert (tmp_path / "checksums.sha256").read_text() == "updated\n"
    assert (
        api.repair_bundle(
            tmp_path,
            required_files=("manifest.json", "environment.yml", "checksums.sha256"),
            commands_finder=lambda **kw: tmp_path / "commands.sh",
            regenerate_manifest=manifest,
            regenerate_environment=environment,
            regenerate_checksums=lambda _: pytest.fail(
                "complete bundle must be untouched"
            ),
        )
        == 0
    )


def test_portable_shell_value_preserves_literal_metacharacters():
    api = importlib.import_module("clawbio.common.nfcore.bundle_runtime")
    value = '/a b/é/$value/`literal`/"quoted"'
    rendered = api.format_output_value(value, old_value='"old"')
    result = subprocess.run(
        ["bash", "-c", "printf '%s' " + rendered],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout == value
