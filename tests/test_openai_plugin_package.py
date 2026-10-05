"""Exercise the distributed archive, including its actual isolated demo runs."""
import importlib.util
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import zipfile

import pytest

ROOT = Path(__file__).resolve().parents[1]


def builder():
    path = ROOT / "scripts/build_openai_plugin.py"
    spec = importlib.util.spec_from_file_location("openai_plugin_builder", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def package(tmp_path_factory):
    directory = tmp_path_factory.mktemp("plugin with spaces")
    archive = builder().build(ROOT, directory / "clawbio.zip")
    with zipfile.ZipFile(archive) as source:
        source.extractall(directory / "installed")
    return directory / "installed"


def test_archive_is_reproducible_and_has_only_curated_skills(tmp_path):
    module = builder()
    first = module.build(ROOT, tmp_path / "first.zip")
    second = module.build(ROOT, tmp_path / "second.zip")
    assert first.read_bytes() == second.read_bytes()
    with zipfile.ZipFile(first) as archive:
        names = archive.namelist()
        manifest = json.loads(archive.read("plugin.json"))
        assert manifest["name"] == "clawbio-research"
        assert manifest["extensions"]["com.openai"]["interface"]["displayName"]
        skills = {n.split("/")[1] for n in names if n.startswith("skills/")}
        assert skills == {"pharmgx-reporter", "clinical-variant-reporter", "equity-scorer"}
        assert all(".." not in Path(n).parts and not n.startswith("/") for n in names)
        assert not any("mcp.json" in n or "/evals/" in n or "/tests/" in n for n in names)
        assert b"Manuel Corpas's publicly available" not in archive.read(
            "runtime/skills/pharmgx-reporter/demo_patient.txt")
        assert not any(n.startswith(("profiles/", "corpas-30x/", "GENOMEBOOK/")) for n in names)


def invoke(package, work, *args):
    env = dict(os.environ)
    # Deny actual network connections in the launcher and its subprocess.
    hook = work / "network_guard"
    hook.mkdir(exist_ok=True)
    (hook / "sitecustomize.py").write_text(
        "import socket\n"
        "def deny(*args, **kwargs):\n"
        "    raise RuntimeError('network forbidden during bundled demos')\n"
        "socket.socket.connect = deny\n"
        "socket.create_connection = deny\n"
    )
    env["PYTHONPATH"] = str(hook)
    env["CLAWBIO_AUDIT_LOG"] = str(work / "audit.jsonl")
    # No optional telemetry exporter should be invoked during package tests.
    env.pop("CLAWBIO_OTLP_ENDPOINT", None)
    return subprocess.run(
        [sys.executable, str(package / "scripts/run_research_demo.py"), *args],
        cwd=work, env=env, capture_output=True, text=True, timeout=90,
    )


@pytest.mark.parametrize("workflow", ["pharmgx", "acmg", "equity"])
def test_extracted_demo_runs_outside_repo_and_writes_provenance(package, tmp_path, workflow):
    output = tmp_path / workflow
    result = invoke(package, tmp_path, workflow, "--output", str(output))
    assert result.returncode == 0, result.stdout + result.stderr
    assert (output / "report.md").is_file()
    assert (output / "result.json").is_file()
    receipt = json.loads((output / "plugin_run.json").read_text())
    assert receipt["mode"] == "bundled-demo"
    assert receipt["workflow"] == workflow
    assert receipt["private_input_supported"] is False
    assert receipt["package_verified"] is True
    assert "not a medical device" in (output / "report.md").read_text()
    assert (output / "reproducibility").is_dir()
    for line in (output / "reproducibility/checksums.sha256").read_text().splitlines():
        checksum, name = line.split("  ", 1)
        assert hashlib.sha256((output / name).read_bytes()).hexdigest() == checksum
    if workflow == "pharmgx":
        data = json.loads((output / "result.json").read_text())
        assert "Indeterminate" in data["data"]["gene_profiles"]["CYP2D6"]["phenotype"]


def test_launcher_rejects_private_input_and_existing_output(package, tmp_path):
    result = invoke(package, tmp_path, "pharmgx", "--input", "private.txt", "--output", str(tmp_path / "out"))
    assert result.returncode != 0
    assert not (tmp_path / "out").exists()
    existing = tmp_path / "existing"
    existing.mkdir()
    marker = existing / "keep.txt"
    marker.write_text("keep")
    result = invoke(package, tmp_path, "equity", "--output", str(existing))
    assert result.returncode != 0
    assert marker.read_text() == "keep"


def test_tampered_runtime_abstains_before_execution(package, tmp_path):
    import shutil
    copy = tmp_path / "tampered"
    shutil.copytree(package, copy)
    code = copy / "runtime/skills/pharmgx-reporter/pharmgx_reporter.py"
    code.write_text("raise RuntimeError('should never execute')")
    output = tmp_path / "out"
    result = invoke(copy, tmp_path, "pharmgx", "--output", str(output))
    assert result.returncode != 0
    assert "PACKAGE_INTEGRITY_FAILED" in result.stderr
    assert not output.exists()
