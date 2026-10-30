"""Run-control values every nf-core wrapper parses from its command line."""

import argparse
import math
from pathlib import Path

import pytest
from hypothesis import example, given
from hypothesis import strategies as st

from clawbio.common.nfcore import cli


@given(
    st.floats(allow_nan=True, allow_infinity=True).map(repr)
    | st.text(alphabet="0123456789.-+einfa", max_size=8)
)
@example("nan")
@example("inf")
@example("-inf")
@example("1e309")
@example("1e306")
@example("-0.0")
def test_accepted_timeouts_always_convert_to_a_usable_wait(value):
    """NaN and infinity never pass as a timeout: int(nan * 3600) raises after
    preflight, so the wrapper died with UNEXPECTED_ERROR instead of a usage error."""
    try:
        hours = cli.timeout_hours(value)
    except argparse.ArgumentTypeError:
        return
    assert math.isfinite(hours)
    assert hours >= 0
    seconds = cli.timeout_seconds(hours, default_hours=12)
    if hours == 0:
        assert seconds is None
    else:
        assert isinstance(seconds, int)
        assert seconds >= 1


@given(st.floats(min_value=0, max_value=1e6))
def test_longer_timeouts_never_wait_less(hours):
    shorter = cli.timeout_seconds(hours, default_hours=12)
    longer = cli.timeout_seconds(hours * 2, default_hours=12)
    if shorter is not None:
        assert longer >= shorter


def test_missing_timeout_falls_back_to_the_default_cap():
    assert cli.timeout_seconds(None, default_hours=12) == 12 * 3600


def test_sub_second_timeout_still_waits():
    # proc.wait(timeout=0) would kill Nextflow on launch.
    assert cli.timeout_seconds(1e-9, default_hours=12) == 1


def test_work_dir_defaults_under_the_output_directory(tmp_path):
    assert cli.resolve_work_dir(None, output_dir=tmp_path) == tmp_path / "upstream" / "work"
    assert cli.resolve_work_dir("", output_dir=tmp_path) == tmp_path / "upstream" / "work"


def test_work_dir_keeps_object_store_uris_verbatim(tmp_path):
    assert cli.resolve_work_dir(" s3://bucket/work ", output_dir=tmp_path) == "s3://bucket/work"


def test_work_dir_resolves_local_paths(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    assert cli.resolve_work_dir("scratch", output_dir=Path("/elsewhere")) == (tmp_path / "scratch").resolve()


@pytest.mark.parametrize("value", ["", "abc", "-1"])
def test_timeout_usage_errors_name_the_flag(value):
    with pytest.raises(argparse.ArgumentTypeError, match="--timeout-hours"):
        cli.timeout_hours(value)
