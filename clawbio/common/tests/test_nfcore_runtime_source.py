"""Runtime facts preserve failure distinctions for local policy adapters, and the
shared preflight checks."""

import importlib
import subprocess
import sys

import pytest
from hypothesis import given
from hypothesis import strategies as st

from clawbio.common.nfcore import runtime
from clawbio.common.nfcore.errors import SkillError


@pytest.mark.parametrize(
    "text,want",
    [
        ('openjdk version "17.0.20.1"', (17, 0, 20)),
        ("Nextflow version 25.04.3 build 1", (25, 4, 3)),
        ("version 17", (17,)),
        ("unavailable", ()),
    ],
)
def test_version_parser(text, want):
    api = importlib.import_module("clawbio.common.nfcore.runtime")
    assert api.parse_version_tuple(text) == want
    assert api.pad_version(want) == want + (0,) * max(0, 3 - len(want))


@given(st.integers(0, 99), st.integers(0, 99), st.integers(0, 99))
def test_version_parser_preserves_three_components(a, b, c):
    api = importlib.import_module("clawbio.common.nfcore.runtime")
    assert api.parse_version_tuple(f"version {a}.{b}.{c}") == (a, b, c)


def test_probe_preserves_nonzero_result_and_stderr():
    api = importlib.import_module("clawbio.common.nfcore.runtime")
    result = api.run_probe(
        [
            sys.executable,
            "-c",
            'import sys; print("version 17",file=sys.stderr); sys.exit(3)',
        ],
        timeout_seconds=3,
    )
    assert result.returncode == 3
    assert result.stdout == ""
    assert result.stderr == "version 17\n"
    assert api.run_probe(["/nonexistent/nfcore-executable"], timeout_seconds=1) is None


def test_probe_timeout_returns_unknown():
    api = importlib.import_module("clawbio.common.nfcore.runtime")
    assert (
        api.run_probe(
            [sys.executable, "-c", "import time; time.sleep(60)"], timeout_seconds=0.01
        )
        is None
    )


def test_host_probe_reports_unknown_invalid_signals(monkeypatch):
    api = importlib.import_module("clawbio.common.nfcore.runtime")
    monkeypatch.setattr(api.os, "sysconf", lambda key: 0)
    assert api.host_memory_gb() is None
    monkeypatch.setattr(
        api.os,
        "sysconf",
        lambda key: {"SC_PHYS_PAGES": 4194304, "SC_PAGE_SIZE": 4096}[key],
    )
    assert api.host_memory_gb() == 16


def test_git_probe_and_required_files_are_explicit(tmp_path):
    api = importlib.import_module("clawbio.common.nfcore.source")
    assert api.missing_pipeline_files(tmp_path, ["main.nf", "nextflow.config"]) == [
        "main.nf",
        "nextflow.config",
    ]
    (tmp_path / "main.nf").touch()
    assert api.missing_pipeline_files(tmp_path, ["main.nf", "nextflow.config"]) == [
        "nextflow.config"
    ]
    assert api.git_stdout(tmp_path, ["rev-parse", "HEAD"]) == ""
    subprocess.run(["git", "init", "-q", str(tmp_path)], check=True)
    assert api.git_stdout(tmp_path, ["status", "--porcelain"]) == "?? main.nf"


def test_raw_host_probe_distinguishes_zero_from_unknown(monkeypatch):
    api = importlib.import_module("clawbio.common.nfcore.runtime")
    monkeypatch.setattr(api.os, "sysconf", lambda key: 0)
    assert api.host_memory_bytes() == 0

    def unavailable(key):
        raise OSError("unknown")

    monkeypatch.setattr(api.os, "sysconf", unavailable)
    assert api.host_memory_bytes() is None


@pytest.mark.parametrize("value", ["true", "1", "YES", " on "])
def test_demo_refuses_to_start_offline(value, monkeypatch):
    monkeypatch.setenv("NXF_OFFLINE", value)
    with pytest.raises(SkillError) as failure:
        runtime.require_network_for_demo()
    assert failure.value.error_code == "DEMO_REQUIRES_NETWORK"
    assert failure.value.details == {"NXF_OFFLINE": value}


@pytest.mark.parametrize("value", [None, "", "false", "0"])
def test_demo_runs_online(value, monkeypatch):
    if value is None:
        monkeypatch.delenv("NXF_OFFLINE", raising=False)
    else:
        monkeypatch.setenv("NXF_OFFLINE", value)
    runtime.require_network_for_demo()


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        ("  N E X T F L O W\n  version 26.04.3 build 5928", "26.04.3"),
        ('openjdk version "17.0.10" 2024-01-16', "17.0.10"),
        ("version 25.10", "25.10"),
        ("no version here", ""),
    ],
)
def test_reported_version_keeps_the_raw_token(text, expected):
    # 26.04.3 must not become 26.4.3: that is not a release NXF_VER can fetch.
    assert runtime.reported_version(text) == expected
