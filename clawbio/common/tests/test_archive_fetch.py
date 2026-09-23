"""Tests for the shared archive-skill machinery.

Run with: pytest clawbio/common/tests/test_archive_fetch.py -v
"""
from __future__ import annotations

import argparse
import types

import pytest

from clawbio.common import archive_fetch as af


def _fake_api(*subcommands: str) -> types.SimpleNamespace:
    """A stand-in for a vendored CLI module, with a real argparse surface."""
    def main(argv=None):
        p = argparse.ArgumentParser(prog="vendored")
        sub = p.add_subparsers(dest="cmd", required=True)
        for name in subcommands:
            sp = sub.add_parser(name)
            sp.add_argument("accession")
        args = p.parse_args(argv)
        print(f"ran {args.cmd} {args.accession}")

    return types.SimpleNamespace(main=main)


class TestRunUpstream:
    def test_it_returns_what_the_vendored_cli_printed(self, tmp_path):
        out, err = af.run_upstream(_fake_api("metadata"), ["metadata", "X1"], tmp_path)
        assert "ran metadata X1" in out
        assert err == ""

    def test_an_unknown_subcommand_surfaces_the_argparse_error(self, tmp_path):
        """The failure this guards against was silent.

        `run_upstream` redirects stderr into a StringIO, so when argparse
        raised SystemExit the `return` was never reached and the captured
        message was discarded -- the caller got exit 2 and no explanation at
        all, while a previous run's report.md stayed on disk looking current.
        """
        with pytest.raises(SystemExit) as excinfo:
            af.run_upstream(_fake_api("metadata"), ["download-script", "X1"], tmp_path)
        assert "invalid choice" in str(excinfo.value), (
            "argparse's message was swallowed by the stderr redirect")

    def test_a_missing_argument_also_surfaces(self, tmp_path):
        with pytest.raises(SystemExit) as excinfo:
            af.run_upstream(_fake_api("metadata"), ["metadata"], tmp_path)
        assert "required" in str(excinfo.value)

    def test_a_clean_exit_is_not_turned_into_an_error(self, tmp_path):
        """SystemExit(0) means the vendored CLI finished, not that it failed."""
        def main(argv=None):
            print("done")
            raise SystemExit(0)

        out, _ = af.run_upstream(types.SimpleNamespace(main=main), ["x"], tmp_path)
        assert "done" in out

    def test_output_paths_are_still_scrubbed(self, tmp_path):
        """Scrubbing must survive the error handling -- reports are compared
        byte-for-byte across output directories by TestDeterminism."""
        def main(argv=None):
            print(f"wrote {tmp_path}/tables/metadata.tsv")

        out, _ = af.run_upstream(types.SimpleNamespace(main=main), ["x"], tmp_path)
        assert str(tmp_path) not in out
