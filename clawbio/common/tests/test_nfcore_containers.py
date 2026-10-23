"""Container provenance: which images a finished run used, pinned to digests."""

from __future__ import annotations

import json
import re
import subprocess

import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from clawbio.common.nfcore import containers

FASTQC = "quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"
FASTQC_DIGEST = "sha256:" + "a" * 64
BCFTOOLS = "community.wave.seqera.io/library/bcftools:1.21--4335bec1d7b44d11"
BCFTOOLS_DIGEST = "sha256:" + "b" * 64
LOG = (
    f"NFCORE_SAREK:SAREK:FASTQC\t{FASTQC}\tCOMPLETED\n"
    f"NFCORE_SAREK:SAREK:FASTQC\t{FASTQC}\tCACHED\n"
    f"NFCORE_SAREK:SAREK:VCF_STATS:BCFTOOLS_STATS\t{BCFTOOLS}\tCOMPLETED\n"
    "NFCORE_SAREK:SAREK:FAILED_ONE\tquay.io/x:1\tFAILED\n"
    "NFCORE_SAREK:SAREK:NATIVE\t\tCOMPLETED\n"
)


def _proc(stdout: str = "", returncode: int = 0) -> subprocess.CompletedProcess[str]:
    return subprocess.CompletedProcess(args=[], returncode=returncode, stdout=stdout, stderr="")


def _inspect(image: str, digest: str) -> str:
    repo = image.rsplit(":", 1)[0]
    return json.dumps({
        "Id": "sha256:" + "f" * 64, "RepoDigests": ["other.io/mirror@sha256:" + "9" * 64, f"{repo}@{digest}"],
        "Os": "linux", "Architecture": "amd64", "Variant": "",
    })


def _run_dir(path):
    history = path / ".nextflow" / "history"
    history.parent.mkdir(parents=True, exist_ok=True)
    history.write_text("2026-10-16 run\n")
    return path


class FakeEngine:
    """Answers the commands the recorder runs; ``local`` maps image -> digest."""

    def __init__(self, *, log: str = LOG, local: dict[str, str] | None = None, registry: dict[str, object] | None = None,
                 log_returncode: int = 0, unreachable: bool = False):
        self.log, self.local = log, local or {}
        self.registry = {} if registry is None else registry
        self.log_returncode = log_returncode
        if unreachable:
            self.registry = None
        self.calls: list[list[str]] = []

    def __call__(self, command, *, cwd=None, env=None, timeout_seconds):
        self.calls.append(list(command))
        if command[:2] == ["nextflow", "log"]:
            return _proc(self.log, self.log_returncode)
        if command[1:3] == ["image", "inspect"]:
            image = command[-1]
            return _proc(_inspect(image, self.local[image])) if image in self.local else _proc("", 1)
        if command[1:3] == ["manifest", "inspect"]:
            if self.registry is None:  # unreachable: the runner timed out
                return None
            image = command[-1]
            return _proc(json.dumps(self.registry[image])) if image in self.registry else _proc("", 1)
        if command[1:2] == ["version"]:
            return _proc("29.1.3 linux/amd64\n")
        raise AssertionError(command)


@pytest.mark.parametrize(("profile", "engine"), [
    ("docker", "docker"), ("test,docker", "docker"), ("docker,arm64", "docker"), ("podman", "podman"),
    ("singularity", "singularity"), ("test,apptainer", "apptainer"), ("conda", None), ("mamba", None), ("", None),
])
def test_engine_comes_from_the_profile(profile, engine):
    assert containers.container_engine(profile) == engine


def test_record_lists_each_image_with_its_digest_and_processes(tmp_path):
    engine = FakeEngine(local={FASTQC: FASTQC_DIGEST, BCFTOOLS: BCFTOOLS_DIGEST})
    record = containers.collect_container_record(_run_dir(tmp_path), profile="test,docker", runner=engine)
    assert record["complete"] is True and record["warnings"] == []
    assert record["engine"] == "docker" and record["engine_version"] == "29.1.3"
    assert record["source"] == "nextflow log" and record["task_count"] == 4
    assert [image["image"] for image in record["images"]] == sorted([FASTQC, BCFTOOLS])
    fastqc = next(image for image in record["images"] if image["image"] == FASTQC)
    assert fastqc == {
        "image": FASTQC,
        "digest": FASTQC_DIGEST,
        "pinned_reference": f"quay.io/biocontainers/fastqc@{FASTQC_DIGEST}",
        "digest_source": "local",
        "image_id": "sha256:" + "f" * 64,
        "platform": "linux/amd64",
        "processes": ["NFCORE_SAREK:SAREK:FASTQC"],
    }
    assert ["nextflow", "log", "last", "-f", "process,container,status"] in engine.calls


def test_a_pruned_image_is_resolved_from_the_registry(tmp_path):
    """Another job may remove an image before the run is recorded."""
    single = {"Descriptor": {"digest": BCFTOOLS_DIGEST, "platform": {"os": "linux", "architecture": "amd64"}}}
    engine = FakeEngine(local={FASTQC: FASTQC_DIGEST}, registry={BCFTOOLS: single})
    record = containers.collect_container_record(_run_dir(tmp_path), profile="docker", runner=engine)
    bcftools = next(image for image in record["images"] if image["image"] == BCFTOOLS)
    assert bcftools["digest"] == BCFTOOLS_DIGEST and bcftools["digest_source"] == "registry"
    assert bcftools["image_id"] is None and record["complete"] is True


def test_a_multi_platform_registry_answer_picks_the_engine_platform(tmp_path):
    arm = {"Descriptor": {"digest": "sha256:" + "c" * 64, "platform": {"os": "linux", "architecture": "arm64"}}}
    amd = {"Descriptor": {"digest": BCFTOOLS_DIGEST, "platform": {"os": "linux", "architecture": "amd64"}}}
    engine = FakeEngine(local={FASTQC: FASTQC_DIGEST}, registry={BCFTOOLS: [arm, amd]})
    record = containers.collect_container_record(_run_dir(tmp_path), profile="docker", runner=engine)
    bcftools = next(image for image in record["images"] if image["image"] == BCFTOOLS)
    assert bcftools["digest"] == BCFTOOLS_DIGEST and bcftools["platform"] == "linux/amd64"


def test_an_unresolved_image_makes_the_record_incomplete_and_withholds_the_pin(tmp_path):
    engine = FakeEngine(local={FASTQC: FASTQC_DIGEST})
    record = containers.collect_container_record(_run_dir(tmp_path), profile="docker", runner=engine)
    assert record["complete"] is False
    assert any(BCFTOOLS in warning for warning in record["warnings"])
    assert containers.pin_config(record) is None
    written = containers.write_container_provenance(tmp_path, record)
    assert [path.name for path in written] == ["containers.json"]


def test_a_digest_reference_is_already_pinned(tmp_path):
    pinned = f"quay.io/biocontainers/fastqc@{FASTQC_DIGEST}"
    engine = FakeEngine(log=f"P:FASTQC\t{pinned}\tCOMPLETED\n")
    record = containers.collect_container_record(_run_dir(tmp_path), profile="docker", runner=engine)
    (image,) = record["images"]
    assert image["digest"] == FASTQC_DIGEST and image["digest_source"] == "reference"
    assert image["pinned_reference"] == pinned


def test_conda_runs_record_no_images(tmp_path):
    engine = FakeEngine(log="P:FASTQC\t\tCOMPLETED\n")
    record = containers.collect_container_record(_run_dir(tmp_path), profile="conda", runner=engine)
    assert record["engine"] is None and record["images"] == [] and record["complete"] is True
    assert containers.pin_config(record) is None


def test_singularity_images_are_named_but_not_pinned(tmp_path):
    engine = FakeEngine(log=f"P:FASTQC\t{FASTQC}\tCOMPLETED\n")
    record = containers.collect_container_record(_run_dir(tmp_path), profile="singularity", runner=engine)
    (image,) = record["images"]
    assert image["image"] == FASTQC and image["digest"] is None and record["complete"] is False
    assert containers.pin_config(record) is None


def test_an_unreadable_task_log_is_reported(tmp_path):
    engine = FakeEngine(log="", log_returncode=1)
    record = containers.collect_container_record(_run_dir(tmp_path), profile="docker", runner=engine)
    assert record["complete"] is False and record["images"] == [] and record["warnings"]


def test_pin_config_pins_every_process_to_its_digest(tmp_path):
    engine = FakeEngine(local={FASTQC: FASTQC_DIGEST, BCFTOOLS: BCFTOOLS_DIGEST})
    record = containers.collect_container_record(_run_dir(tmp_path), profile="docker", runner=engine)
    config = containers.pin_config(record)
    assert config is not None
    pins = dict(re.findall(r"withName: '([^']+)' \{\n\s+container = '([^']+)'", config))
    assert pins == {
        "NFCORE_SAREK:SAREK:FASTQC": f"quay.io/biocontainers/fastqc@{FASTQC_DIGEST}",
        "NFCORE_SAREK:SAREK:VCF_STATS:BCFTOOLS_STATS": f"community.wave.seqera.io/library/bcftools@{BCFTOOLS_DIGEST}",
    }
    written = containers.write_container_provenance(tmp_path, record)
    assert [path.name for path in written] == ["containers.json", "containers.config"]
    assert json.loads((tmp_path / "containers.json").read_text()) == record
    assert (tmp_path / "containers.config").read_bytes().count(b"\r") == 0


@settings(deadline=None, max_examples=60, suppress_health_check=[HealthCheck.function_scoped_fixture])
@given(st.permutations(LOG.strip("\n").split("\n")))
def test_the_record_does_not_depend_on_task_order(tmp_path, lines):
    run_dir = _run_dir(tmp_path)
    engine = FakeEngine(log="\n".join(lines) + "\n", local={FASTQC: FASTQC_DIGEST, BCFTOOLS: BCFTOOLS_DIGEST})
    reference = FakeEngine(local={FASTQC: FASTQC_DIGEST, BCFTOOLS: BCFTOOLS_DIGEST})
    shuffled = containers.collect_container_record(run_dir, profile="docker", runner=engine)
    ordered = containers.collect_container_record(run_dir, profile="docker", runner=reference)
    assert shuffled == ordered
    assert containers.pin_config(shuffled) == containers.pin_config(ordered)


def test_a_process_that_ran_two_images_cannot_be_pinned(tmp_path):
    other = "quay.io/biocontainers/fastqc:0.11.9--0"
    engine = FakeEngine(
        log=f"P:FASTQC\t{FASTQC}\tCOMPLETED\nP:FASTQC\t{other}\tCOMPLETED\n",
        local={FASTQC: FASTQC_DIGEST, other: "sha256:" + "d" * 64},
    )
    record = containers.collect_container_record(_run_dir(tmp_path), profile="docker", runner=engine)
    assert record["complete"] is False and any("P:FASTQC" in warning for warning in record["warnings"])
    assert containers.pin_config(record) is None


def test_manifest_summary_points_at_the_bundle_files(tmp_path):
    engine = FakeEngine(local={FASTQC: FASTQC_DIGEST, BCFTOOLS: BCFTOOLS_DIGEST})
    record = containers.record_run_containers(_run_dir(tmp_path), profile="docker", runner=engine)
    assert (tmp_path / "reproducibility" / "containers.json").is_file()
    assert containers.manifest_summary(record) == {
        "engine": "docker", "engine_version": "29.1.3", "images": 2, "complete": True,
        "record": "reproducibility/containers.json", "pin_config": "reproducibility/containers.config",
    }
    engine = FakeEngine(local={FASTQC: FASTQC_DIGEST})
    record = containers.record_run_containers(_run_dir(tmp_path), profile="docker", runner=engine)
    assert not (tmp_path / "reproducibility" / "containers.config").exists()
    assert containers.manifest_summary(record)["pin_config"] is None


def _run_guard(tmp_path, lines: str, *, resume: str, pinned: bool) -> str:
    repro = tmp_path / "reproducibility"
    repro.mkdir(exist_ok=True)
    pin = repro / "containers.config"
    if pinned:
        pin.write_text("process {}\n")
    else:
        pin.unlink(missing_ok=True)
    script = repro / "probe.sh"
    script.write_text(
        "#!/bin/bash\nset -euo pipefail\n"
        'SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"\n'
        f'RESUME="{resume}"\n' + lines + 'printf "[%s]" $CONTAINER_PIN\n'
    )
    return subprocess.run(["/bin/bash", str(script)], capture_output=True, text=True, check=True).stdout


def test_replay_pins_images_only_on_a_fresh_replay(tmp_path):
    lines = containers.replay_pin_guard(resume_variable="RESUME")
    assert _run_guard(tmp_path, lines, resume="", pinned=True) == "[-c][reproducibility/containers.config]"
    assert _run_guard(tmp_path, lines, resume="-resume", pinned=True) == "[]"
    assert _run_guard(tmp_path, lines, resume="", pinned=False) == "[]"


def test_a_directory_without_a_nextflow_run_starts_no_process(tmp_path):
    """Bundle writers also run in unit tests; without a Nextflow history there is
    nothing to ask, so the record says so without launching Nextflow or Docker."""
    engine = FakeEngine()
    record = containers.collect_container_record(tmp_path, profile="docker", runner=engine)
    assert engine.calls == [] and record["images"] == [] and record["complete"] is False
    assert "history" in record["warnings"][0]


def _manifest_calls(engine):
    return [call for call in engine.calls if call[1:3] == ["manifest", "inspect"]]


@pytest.mark.parametrize("value", ["true", "1", "YES"])
def test_offline_hosts_never_query_a_registry(tmp_path, monkeypatch, value):
    """Under NXF_OFFLINE every image without a local digest waited out a 120 s
    registry query, one after another, after the run had finished."""
    monkeypatch.setenv("NXF_OFFLINE", value)
    engine = FakeEngine(local={}, registry={FASTQC: {"Descriptor": {"digest": FASTQC_DIGEST}}})
    record = containers.collect_container_record(_run_dir(tmp_path), profile="docker", runner=engine)
    assert _manifest_calls(engine) == []
    assert record["complete"] is False
    assert any("NXF_OFFLINE" in warning for warning in record["warnings"])


def test_an_unreachable_registry_is_queried_once(tmp_path, monkeypatch):
    monkeypatch.delenv("NXF_OFFLINE", raising=False)
    engine = FakeEngine(local={}, unreachable=True)
    record = containers.collect_container_record(_run_dir(tmp_path), profile="docker", runner=engine)
    assert len(_manifest_calls(engine)) == 1
    assert record["complete"] is False
    assert sum("registry" in warning and "unreachable" in warning for warning in record["warnings"]) == 1
