"""Emitted single-file, standard-library remappers shipped in every bundle."""

import csv
import hashlib
import importlib
import subprocess
import sys
from pathlib import Path

import pytest


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
    (bundle.parent / "upstream" / "results").mkdir(parents=True)
    (bundle.parent / "upstream" / "results" / "calls.vcf").write_text("synthetic\n")

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
