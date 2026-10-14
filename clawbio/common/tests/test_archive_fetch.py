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


class TestSafeJoin:
    """A file path supplied by an archive API must not escape the output dir."""

    def test_a_nested_relative_path_stays_under_the_base(self, tmp_path):
        assert af.safe_join(tmp_path, "sub/ok.txt") == (tmp_path / "sub" / "ok.txt").resolve()

    @pytest.mark.parametrize("hostile", ["../evil", "a/../../evil", "/etc/passwd", ".."])
    def test_an_escaping_path_is_refused(self, tmp_path, hostile):
        with pytest.raises(SystemExit, match="outside"):
            af.safe_join(tmp_path, hostile)

    def test_a_dot_dot_that_stays_inside_is_allowed(self, tmp_path):
        assert af.safe_join(tmp_path, "a/../b.txt") == (tmp_path / "b.txt").resolve()

    def test_a_symlink_out_of_the_base_is_refused(self, tmp_path):
        outside = tmp_path / "outside"
        outside.mkdir()
        base = tmp_path / "base"
        base.mkdir()
        (base / "link").symlink_to(outside)
        with pytest.raises(SystemExit, match="outside"):
            af.safe_join(base, "link/x.txt")


class TestTrustedEbiBase:
    """`httpLink` from /studies/{acc}/info decides where bytes come from, so
    only an https URL on an ebi.ac.uk host is accepted."""

    @pytest.mark.parametrize("link", [
        "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/030/E-MTAB-10030",
        "https://www.ebi.ac.uk/biostudies/files/S-B1",
        "https://ebi.ac.uk/x",
    ])
    def test_accepts_https_on_ebi(self, link):
        assert af.trusted_ebi_base(link) == link

    @pytest.mark.parametrize("link", [
        None, "", "https://evil.example/x", "https://ebi.ac.uk.evil.com/x",
        "https://evilebi.ac.uk/x", "http://ftp.ebi.ac.uk/x", "file:///etc",
        "ftp://ftp.ebi.ac.uk/x", "https://user@evil.example/x?ebi.ac.uk",
        "-o/etc/passwd", 42,
    ])
    def test_rejects_everything_else(self, link):
        assert af.trusted_ebi_base(link) is None


class TestReplayArgv:
    """commands.sh must hold the command that was run, not a summary of it."""

    def test_every_flag_is_kept_and_output_is_re_anchored(self, tmp_path):
        argv = ["--command", "download", "--accession", "GSE1", "--suppl",
                "--output", "rel/out"]
        assert af.replay_argv(argv, tmp_path) == [
            "--command", "download", "--accession", "GSE1", "--suppl",
            "--output", str(tmp_path)]

    def test_the_equals_form_of_output_is_replaced_too(self, tmp_path):
        assert af.replay_argv(["--demo", "--output=x"], tmp_path) == [
            "--demo", "--output", str(tmp_path)]

    def test_out_is_not_mistaken_for_output(self, tmp_path):
        argv = ["--command", "samplesheet", "--out", "s.csv", "--output", "o"]
        assert af.replay_argv(argv, tmp_path) == [
            "--command", "samplesheet", "--out", "s.csv", "--output", str(tmp_path)]

    def test_write_bundle_records_the_replay_command(self, tmp_path):
        script = tmp_path / "skill.py"
        af.write_bundle(tmp_path, script=script, env_name="e",
                        argv=["--command", "runs", "--accession", "PRJ 1", "--output", "o"],
                        written=[])
        line = (tmp_path / "reproducibility" / "commands.sh").read_text().splitlines()[1]
        assert line == (f"python {script} --command runs --accession 'PRJ 1' "
                        f"--output {tmp_path}")
