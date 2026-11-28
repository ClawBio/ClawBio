"""Tests for environment_check.py. Every probe is fed fixed inputs, so no test depends on the
machine it runs on."""

import importlib.util
import json
import sys
from pathlib import Path

import pytest

SKILL_DIR = Path(__file__).resolve().parent.parent
MODULE = SKILL_DIR / "environment_check.py"


def load():
    spec = importlib.util.spec_from_file_location("environment_check", MODULE)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["environment_check"] = mod
    spec.loader.exec_module(mod)
    return mod


ec = load()

MOUNTS_9P = "none / 9p rw 0 0\nnone /home/u 9p rw 0 0\ntmpfs /tmp tmpfs rw 0 0\n"
MOUNTS_EXT4 = "/dev/sda1 / ext4 rw 0 0\n/dev/sdb1 /data ext4 rw 0 0\n"
MEM_128 = "MemTotal:       134217728 kB\n"


def tools(*names):
    present = {n: f"/usr/bin/{n}" for n in names}
    return lambda n: present.get(n)


def probes(**kw):
    base = {
        "kernel": "6.1.0", "version_text": "Linux version 6.1.0", "meminfo": MEM_128,
        "mounts": MOUNTS_EXT4, "cpus": 32, "dev_fuse": True, "disk_free_gb": 900.0,
        "run_version": lambda cmd, pat, timeout=30: "1.0",
        "apptainer_conf": "", "docker_ok": None,
        "which": tools("java", "nextflow", "uv"),
    }
    base.update(kw)
    return base


def test_docker_with_a_live_daemon_gets_docker(tmp_path):
    rep = ec.detect(tmp_path, probes=probes(which=tools("java", "nextflow", "docker"),
                                            docker_ok=True))
    assert rep.profile == "docker"


def test_docker_with_a_dead_daemon_falls_to_apptainer(tmp_path):
    rep = ec.detect(tmp_path, probes=probes(
        which=tools("java", "nextflow", "docker", "apptainer", "fusermount3"), docker_ok=False))
    assert rep.profile == "singularity-sif"


def test_plain_linux_host_with_apptainer_gets_sif(tmp_path):
    rep = ec.detect(tmp_path, probes=probes(
        which=tools("java", "nextflow", "apptainer", "fusermount3")))
    assert rep.profile == "singularity-sif"
    assert any("--profile singularity" in s for s in rep.next_steps)


def test_gvisor_host_gets_sandbox_and_paths_use_the_callers_home(tmp_path):
    rep = ec.detect(tmp_path, probes=probes(
        kernel="4.19.0-gvisor", mounts=MOUNTS_9P, dev_fuse=False,
        which=tools("java", "nextflow", "apptainer")))
    assert rep.profile == "singularity-sandbox"
    joined = " ".join(rep.next_steps)
    assert str(Path.home()) in joined and "--extra-config" in joined
    assert "/home/ccloud" not in joined or str(Path.home()) == "/home/ccloud"


def test_overlay_turned_off_is_a_sandbox_signal(tmp_path):
    rep = ec.detect(tmp_path, probes=probes(
        apptainer_conf="enable overlay = no\n",
        which=tools("java", "nextflow", "apptainer", "fusermount3")))
    assert rep.profile == "singularity-sandbox"


def test_no_runtime_gives_none_with_install_advice(tmp_path):
    rep = ec.detect(tmp_path, probes=probes())
    assert rep.profile == "none"
    assert "Install Docker" in rep.next_steps[0]


def test_successful_sif_probe_overrides_weak_signals(tmp_path):
    img = tmp_path / "x.sif"
    img.write_text("")
    rep = ec.detect(tmp_path, image_to_probe=img, probes=probes(
        dev_fuse=False, which=tools("java", "nextflow", "apptainer"),
        probe_image=lambda rt, im: (True, "ok")))
    assert rep.profile == "singularity-sif"


def test_failed_probe_forces_sandbox(tmp_path):
    img = tmp_path / "x.sif"
    img.write_text("")
    rep = ec.detect(tmp_path, image_to_probe=img, probes=probes(
        which=tools("java", "nextflow", "apptainer", "fusermount3"),
        probe_image=lambda rt, im: (False, "mount failed")))
    assert rep.profile == "singularity-sandbox"
    assert any(f.status == "problem" and "probe" in f.name for f in rep.findings)


def test_low_disk_and_missing_java_are_reported(tmp_path):
    rep = ec.detect(tmp_path, probes=probes(disk_free_gb=50.0, which=tools("nextflow")))
    by = {f.name: f for f in rep.findings}
    assert by["free disk at output"].status == "warn"
    assert by["java"].status == "problem"
    assert rep.worst() == "problem"


def test_min_free_gb_is_adjustable(tmp_path):
    rep = ec.detect(tmp_path, min_free_gb=20, probes=probes(disk_free_gb=50.0))
    assert {f.name: f for f in rep.findings}["free disk at output"].status == "ok"


def test_filesystem_lookup_uses_longest_mount():
    assert ec.filesystem_of(Path("/data/run1"), MOUNTS_EXT4) == "ext4"
    assert ec.filesystem_of(Path("/home/u/x"), MOUNTS_9P) == "9p"
    assert ec.filesystem_of(Path("/elsewhere"), "") == "unknown"


def test_9p_output_warns(tmp_path):
    rep = ec.detect(tmp_path, probes=probes(mounts=f"none / 9p rw 0 0\n"))
    assert {f.name: f for f in rep.findings}["filesystem at output"].status == "warn"


def test_overlay_setting_parser():
    assert ec.apptainer_overlay_setting("# c\nenable overlay = try\n") == "try"
    assert ec.apptainer_overlay_setting("") is None


def test_memory_parser():
    assert ec.mem_total_gb(MEM_128) == 128.0
    assert ec.mem_total_gb("") is None


def test_report_text_and_json_are_consistent(tmp_path):
    rep = ec.detect(tmp_path, probes=probes())
    text = ec.format_report(rep)
    assert "Recommended profile: none" in text
    data = ec.to_json(rep)
    assert data["profile"] == "none" and json.dumps(data)


def test_detect_changes_nothing(tmp_path):
    before = sorted(p.name for p in tmp_path.iterdir())
    ec.detect(tmp_path / "not-yet-created", probes=probes())
    assert sorted(p.name for p in tmp_path.iterdir()) == before
    assert not (tmp_path / "not-yet-created").exists()
