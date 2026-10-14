from __future__ import annotations

import os
import subprocess

import pytest


@pytest.mark.parametrize("system", ["Darwin", "Linux"])
def test_resource_replay_bounds_capacity_and_preserves_snapshot(tmp_path, system):
    from clawbio.common import nfcore

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
