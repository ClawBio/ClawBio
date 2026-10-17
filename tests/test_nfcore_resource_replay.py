from __future__ import annotations

import os
import subprocess

import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st


@pytest.mark.parametrize("system", ["Darwin", "Linux"])
def test_resource_replay_bounds_capacity_and_preserves_snapshot(tmp_path, system):
    from clawbio.common.nfcore import resources as nfcore

    repro = tmp_path / "reproducibility"
    repro.mkdir()
    source = repro / "resource_limits.config"
    original = "executor.cpus = 64\nexecutor.memory = '64.GB'\nprocess {\n    resourceLimits = [\n        cpus: 64,\n        memory: '64.GB',\n        time: '24.h'\n    ]\n}\n"
    source.write_text(original)
    binaries = tmp_path / "bin"
    binaries.mkdir()
    (binaries / "uname").write_text(f'#!/bin/sh\nprintf "{system}\\n"\n')
    (binaries / "docker").write_text('#!/bin/sh\nprintf "2 2147483648\\n"\n')
    for path in binaries.iterdir():
        path.chmod(0o755)
    script = tmp_path / "test.sh"
    script.write_text(
        "#!/usr/bin/env bash\nset -euo pipefail\n"
        + nfcore.build_runtime_resource_guard()
        + "\ncp reproducibility/resource_limits.replay.config captured.config\n"
    )
    result = subprocess.run(
        ["bash", str(script)],
        cwd=tmp_path,
        env={**os.environ, "PATH": str(binaries) + os.pathsep + os.environ["PATH"]},
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    captured = (tmp_path / "captured.config").read_text()
    assert "executor.cpus = 2" in captured
    assert "executor.memory = '1.GB'" in captured
    assert "memory: '1.GB'" in captured
    assert source.read_text() == original
    assert not (repro / "resource_limits.replay.config").exists()


def _replay_config(tmp_path, original: str, *, runtime_cpus: int, runtime_bytes: int):
    """Run the replay guard against ``original`` with a fake Docker; return
    (exit code, stderr, bounded replay config or None)."""
    from clawbio.common.nfcore import resources as nfcore

    repro = tmp_path / "reproducibility"
    repro.mkdir(exist_ok=True)
    (repro / "resource_limits.config").write_text(original)
    binaries = tmp_path / "bin"
    binaries.mkdir(exist_ok=True)
    (binaries / "docker").write_text(f'#!/bin/sh\nprintf "{runtime_cpus} {runtime_bytes}\\n"\n')
    (binaries / "docker").chmod(0o755)
    script = tmp_path / "guard.sh"
    script.write_text(
        "#!/usr/bin/env bash\nset -euo pipefail\n"
        + nfcore.build_runtime_resource_guard()
        + "\ncp reproducibility/resource_limits.replay.config captured.config\n"
    )
    captured = tmp_path / "captured.config"
    captured.unlink(missing_ok=True)
    result = subprocess.run(
        ["bash", str(script)], cwd=tmp_path, capture_output=True, text=True, check=False,
        env={**os.environ, "PATH": str(binaries) + os.pathsep + os.environ["PATH"]},
    )
    return result.returncode, result.stderr, captured.read_text() if captured.exists() else None


_GB = 1024**3


def test_replay_reads_fractional_memory_without_dropping_the_point(tmp_path):
    """A hand-edited 15.5 GB used to be read as 155 GB: every non-digit was stripped."""
    original = (
        "executor.cpus = 8\nexecutor.memory = '15.5.GB'\n"
        "process {\n    resourceLimits = [\n        cpus: 8,\n        memory: '15.5.GB',\n        time: '24.h'\n    ]\n}\n"
    )
    rc, err, replay = _replay_config(tmp_path, original, runtime_cpus=16, runtime_bytes=64 * _GB)
    assert rc == 0, err
    assert "executor.memory = '15872.MB'" in replay and "memory: '15872.MB'" in replay
    rc, err, replay = _replay_config(tmp_path, original, runtime_cpus=16, runtime_bytes=10 * _GB)
    assert rc == 0, err
    assert "executor.memory = '9.GB'" in replay


def test_replay_reads_recorded_limits_from_the_executor_lines(tmp_path):
    """Other lines that mention cpus or memory must not be taken for the limits."""
    original = (
        "// cpus: 128 and memory: '512.GB' are the cluster's, not this run's\n"
        "executor.cpus = 4\nexecutor.memory = '6.GB'\n"
        "process {\n    resourceLimits = [\n        cpus: 4,\n        memory: '6.GB',\n        time: '24.h'\n    ]\n}\n"
    )
    rc, err, replay = _replay_config(tmp_path, original, runtime_cpus=32, runtime_bytes=128 * _GB)
    assert rc == 0, err
    assert "executor.cpus = 4" in replay and "executor.memory = '6.GB'" in replay
    assert "        cpus: 4," in replay and "        memory: '6.GB'," in replay


def test_replay_refuses_unreadable_recorded_memory(tmp_path):
    original = "executor.cpus = 4\nexecutor.memory = 'lots'\n"
    rc, err, replay = _replay_config(tmp_path, original, runtime_cpus=4, runtime_bytes=8 * _GB)
    assert rc != 0 and replay is None and "resource" in err.lower()



@settings(deadline=None, max_examples=40, suppress_health_check=[HealthCheck.function_scoped_fixture])
@given(
    recorded_cpus=st.integers(1, 256),
    recorded_tenths=st.integers(1, 4096),
    runtime_cpus=st.integers(1, 256),
    runtime_gb=st.integers(1, 1024),
)
def test_replay_never_exceeds_recorded_or_runtime_limits(tmp_path, recorded_cpus, recorded_tenths, runtime_cpus, runtime_gb):
    recorded_gb = recorded_tenths / 10
    text = f"{recorded_gb:g}"
    original = f"executor.cpus = {recorded_cpus}\nexecutor.memory = '{text}.GB'\n"
    rc, err, replay = _replay_config(tmp_path, original, runtime_cpus=runtime_cpus, runtime_bytes=runtime_gb * _GB)
    assert rc == 0, err
    cpus = int(replay.split("executor.cpus = ")[1].split("\n")[0])
    memory = replay.split("executor.memory = '")[1].split("'")[0]
    number, unit = memory.rsplit(".", 1)
    gb = float(number) / (1024 if unit == "MB" else 1)
    assert cpus == min(recorded_cpus, runtime_cpus)
    limit = max(1, int(runtime_gb * 0.9))
    assert abs(gb - min(recorded_gb, limit)) <= 1 / 1024


def test_replay_bounds_aggregate_and_per_process_limits_separately(tmp_path):
    """A demo records the host's aggregate executor budget but the test profile's
    smaller per-process limits; replay must bound each by the runtime, never copy the
    aggregate into resourceLimits (which would lift a demo's 4 CPUs to the host's)."""
    original = (
        "executor.cpus = 28\nexecutor.memory = '58.GB'\n"
        "process {\n    resourceLimits = [\n        cpus: 4,\n        memory: '15.GB',\n        time: '24.h'\n    ]\n}\n"
    )
    rc, err, replay = _replay_config(tmp_path, original, runtime_cpus=16, runtime_bytes=64 * _GB)
    assert rc == 0, err
    assert "executor.cpus = 16" in replay and "executor.memory = '57.GB'" in replay
    assert "        cpus: 4," in replay and "        memory: '15.GB'," in replay
    rc, err, replay = _replay_config(tmp_path, original, runtime_cpus=2, runtime_bytes=8 * _GB)
    assert rc == 0, err
    assert "executor.cpus = 2" in replay and "        cpus: 2," in replay
    assert "executor.memory = '7.GB'" in replay and "        memory: '7.GB'," in replay


@settings(deadline=None, max_examples=40, suppress_health_check=[HealthCheck.function_scoped_fixture])
@given(
    aggregate=st.integers(1, 256), per_process=st.integers(1, 256), memory_tenths=st.integers(1, 4096),
    limit_tenths=st.integers(1, 4096), runtime_cpus=st.integers(1, 256), runtime_gb=st.integers(1, 1024),
)
def test_replay_per_process_limits_never_exceed_recorded_or_runtime(
    tmp_path, aggregate, per_process, memory_tenths, limit_tenths, runtime_cpus, runtime_gb
):
    original = (
        f"executor.cpus = {aggregate}\nexecutor.memory = '{memory_tenths / 10:g}.GB'\n"
        f"process {{\n    resourceLimits = [\n        cpus: {per_process},\n"
        f"        memory: '{limit_tenths / 10:g}.GB',\n        time: '24.h'\n    ]\n}}\n"
    )
    rc, err, replay = _replay_config(tmp_path, original, runtime_cpus=runtime_cpus, runtime_bytes=runtime_gb * _GB)
    assert rc == 0, err
    limit = max(1, int(runtime_gb * 0.9))

    def gb(text):
        number, unit = text.rsplit(".", 1)
        return float(number) / (1024 if unit == "MB" else 1)

    cpus = int(replay.split("        cpus: ")[1].split(",")[0])
    memory = gb(replay.split("        memory: '")[1].split("'")[0])
    assert cpus == min(per_process, runtime_cpus)
    assert abs(memory - min(limit_tenths / 10, limit)) <= 1 / 1024
