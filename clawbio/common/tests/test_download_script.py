"""Tests for the shared FASTQ download-script emitter.

Run with: pytest clawbio/common/tests/test_download_script.py -v

The emitter produces a bash script; it never downloads anything itself.
"""

from __future__ import annotations

import csv
import stat
import subprocess

import pytest

from clawbio.common.download_script import (
    SlurmOptions,
    build_download_script,
    read_md5_sidecar,
    urls_from_samplesheet,
    urls_from_url_list,
    write_download_script,
)
from clawbio.common.tests.shims import run_with_shims


def _samplesheet(tmp_path, rows, header=("sample", "fastq_1", "fastq_2")):
    path = tmp_path / "samplesheet.csv"
    with path.open("w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)
    return path


class TestSamplesheetParsing:
    def test_reads_fastq_columns_in_row_order(self, tmp_path):
        sheet = _samplesheet(tmp_path, [
            ["s1", "https://x/1_1.fq.gz", "https://x/1_2.fq.gz"],
            ["s2", "https://x/2_1.fq.gz", ""],
        ])
        assert urls_from_samplesheet(sheet) == [
            ("s1", ["https://x/1_1.fq.gz", "https://x/1_2.fq.gz"]),
            ("s2", ["https://x/2_1.fq.gz"]),
        ]

    def test_rejects_a_sheet_with_no_fastq_columns(self, tmp_path):
        sheet = _samplesheet(tmp_path, [["s1", "a", "b"]], header=("sample", "x", "y"))
        with pytest.raises(SystemExit, match="fastq"):
            urls_from_samplesheet(sheet)

    def test_skips_local_paths_and_keeps_urls(self, tmp_path):
        sheet = _samplesheet(tmp_path, [["s1", "/local/1.fq.gz", "https://x/1_2.fq.gz"]])
        assert urls_from_samplesheet(sheet) == [("s1", ["https://x/1_2.fq.gz"])]

    def test_skips_a_url_carrying_a_quote_or_control_character(self, tmp_path):
        """Neither belongs in a URL; a control character would also split the line."""
        sheet = _samplesheet(tmp_path, [["s1", 'https://x/a".fq.gz', "https://x/ok.fq.gz"]])
        assert urls_from_samplesheet(sheet) == [("s1", ["https://x/ok.fq.gz"])]

    def test_falls_back_to_a_row_label_when_sample_is_missing(self, tmp_path):
        sheet = _samplesheet(tmp_path, [["", "https://x/1.fq.gz"]],
                             header=("sample", "fastq_1"))
        assert urls_from_samplesheet(sheet)[0][0] == "row1"

    def test_url_list_ignores_comments_and_blanks(self, tmp_path):
        path = tmp_path / "urls.txt"
        path.write_text("# a comment\n\nhttps://x/1.fq.gz\nnot-a-url\n")
        assert urls_from_url_list(path) == [("file3", ["https://x/1.fq.gz"])]


class TestScriptBody:
    def test_wget_and_curl_are_both_quiet_and_retrying(self):
        groups = [("s1", ["https://x/1.fq.gz"])]
        wget, _ = build_download_script(groups, tool="wget", outdir="fastq")
        curl, _ = build_download_script(groups, tool="curl", outdir="fastq")
        assert "wget -q --tries=5" in wget
        assert "curl -fsSL --retry 5" in curl

    def test_curl_is_the_default_tool(self):
        """curl resumes and retries 403; wget is the opt-in alternative."""
        body, _ = build_download_script([("s1", ["https://x/1.fq.gz"])])
        assert "curl -fsSL" in body
        assert "wget" not in body

    def test_wget_resumes_a_partial_transfer(self):
        """-c verified live 2026-09-22: 206 Partial Content, only the remainder.

        Not -nc: that means --no-clobber and would SKIP an existing partial,
        stranding it truncated rather than finishing it.
        """
        body, _ = build_download_script([("s1", ["https://x/1.fq.gz"])], tool="wget")
        assert " -c " in body
        assert "-nc" not in body

    def test_curl_resumes_a_partial_transfer(self):
        """Archive FASTQs are multi-GB; a drop at 90% must not restart at 0."""
        body, _ = build_download_script([("s1", ["https://x/1.fq.gz"])], tool="curl")
        assert "-C -" in body

    def test_both_tools_bound_the_connect_time(self):
        """An unattended job must not hang forever on a dead host."""
        groups = [("s1", ["https://x/1.fq.gz"])]
        curl, _ = build_download_script(groups, tool="curl")
        wget, _ = build_download_script(groups, tool="wget")
        assert "--connect-timeout 30" in curl
        assert "--timeout=30" in wget

    def test_curl_aborts_a_stalled_transfer(self):
        """A transfer trickling at ~0 B/s would otherwise hang until walltime."""
        body, _ = build_download_script([("s1", ["https://x/1.fq.gz"])], tool="curl")
        assert "--speed-limit 1024" in body
        assert "--speed-time 120" in body

    def test_curl_waits_between_retries(self):
        body, _ = build_download_script([("s1", ["https://x/1.fq.gz"])], tool="curl")
        assert "--retry-delay 10" in body

    def test_curl_probes_for_retry_all_errors(self):
        """--retry-all-errors is curl >= 7.71; probe rather than assume.

        Without it curl retries only transient network failures, not an HTTP
        error -- and archives do answer 403/429/503 under load. Hardcoding the
        flag would break the script outright on an older curl (RHEL 7 ships
        7.29), so the script detects support at run time.
        """
        body, _ = build_download_script([("s1", ["https://x/1.fq.gz"])], tool="curl")
        assert "--retry-all-errors" in body
        assert "RETRY_ALL" in body

    def test_wget_script_has_no_curl_probe(self):
        body, _ = build_download_script([("s1", ["https://x/1.fq.gz"])], tool="wget")
        assert "RETRY_ALL" not in body

    def test_script_is_strict_bash(self):
        body, _ = build_download_script([("s1", ["https://x/1.fq.gz"])], tool="wget")
        assert "set -euo pipefail" in body

    def test_counts_every_file(self):
        groups = [("s1", ["https://x/1.fq.gz", "https://x/2.fq.gz"]), ("s2", ["https://x/3.fq.gz"])]
        _, n = build_download_script(groups, tool="wget")
        assert n == 3

    def test_no_slurm_header_by_default_request(self):
        body, _ = build_download_script([("s1", ["https://x/1.fq.gz"])], tool="wget", slurm=None)
        assert "#SBATCH" not in body
        assert body.startswith("#!/bin/bash")


class TestSlurmHeader:
    def test_partition_is_commented_out_when_unset(self):
        """Slurm has no standard default partition name; emitting a guess makes
        sbatch fail with 'invalid partition specified'. Omitting the directive
        lets the controller pick the site default, which is correct everywhere."""
        body, _ = build_download_script(
            [("s1", ["https://x/1.fq.gz"])], tool="wget", slurm=SlurmOptions())
        assert "# #SBATCH --partition=<your_partition>" in body
        assert "\n#SBATCH --partition=" not in body

    def test_partition_is_emitted_when_given(self):
        body, _ = build_download_script(
            [("s1", ["https://x/1.fq.gz"])], tool="wget",
            slurm=SlurmOptions(partition="batch"))
        assert "\n#SBATCH --partition=batch" in body

    def test_account_and_email_follow_the_same_rule(self):
        body, _ = build_download_script(
            [("s1", ["https://x/1.fq.gz"])], tool="wget", slurm=SlurmOptions())
        assert "# #SBATCH --account=<your_account>" in body
        assert "# #SBATCH --mail-user=<you@example.org>" in body

    def test_carries_no_site_specific_defaults(self):
        """No site names may reach a shipped file."""
        body, _ = build_download_script(
            [("s1", ["https://x/1.fq.gz"])], tool="wget", slurm=SlurmOptions())
        lowered = body.lower()
        for token in ("nfsdata", "ukdri", "htc"):
            assert token not in lowered

    @pytest.mark.parametrize("field", ["job_name", "partition", "account", "cpus",
                                       "mem", "time", "email"])
    def test_values_cannot_break_out_of_their_line(self, field):
        """A newline in a header value would end the #SBATCH line and start a
        command that runs when the job does."""
        opts = SlurmOptions(job_name="j", partition="p", account="a",
                            email="me@example.org")
        setattr(opts, field, "x\nrm -rf ~\r\n#SBATCH --wrap=evil")
        body, _ = build_download_script(
            [("s1", ["https://x/1.fq.gz"])], tool="wget", slurm=opts)
        lines = body.splitlines()
        assert not any(line.startswith(("rm", "#SBATCH --wrap")) for line in lines)


class TestWrite:
    def test_written_script_is_executable(self, tmp_path):
        out = tmp_path / "download_fastq.sh"
        write_download_script([("s1", ["https://x/1.fq.gz"])], out, tool="wget")
        assert out.exists()
        assert out.stat().st_mode & stat.S_IXUSR

    def test_refuses_to_write_an_empty_script(self, tmp_path):
        with pytest.raises(SystemExit, match="No FASTQ URLs"):
            write_download_script([], tmp_path / "x.sh", tool="wget")

    def test_refuses_when_every_url_was_skipped(self, tmp_path):
        with pytest.raises(SystemExit, match="No FASTQ URLs"):
            write_download_script([("s1", ["https://x/.."])], tmp_path / "x.sh")


# Values bash would expand inside double quotes. SDRF Comment[FASTQ_URI] cells
# are submitter free text and reach these scripts via ena-fetch, so a URL is
# not trusted input.
HOSTILE_URLS = [
    "https://x/a$(touch PWNED_SUBST).fq.gz",
    "https://x/b`touch PWNED_TICK`.fq.gz",
    "https://x/c${HOME}.fq.gz",
    "https://x/d\\$HOME.fq.gz",
    "https://x/e.fq.gz?x=1&y=a+b;touch PWNED_SEMI",
]


def _transfers(calls):
    """The download calls, without curl's --help probe."""
    return [c for c in calls if "--help" not in c]


@pytest.mark.parametrize("tool", ["curl", "wget"])
class TestShellInjection:
    """A generated script must pass every value to the tool literally."""

    def _run(self, tmp_path, urls, tool, outdir="fastq"):
        script = tmp_path / "dl.sh"
        write_download_script([("s1", urls)], script, tool=tool, outdir=outdir)
        return run_with_shims(script, tmp_path)

    def test_no_url_can_run_a_command(self, tmp_path, tool):
        result, _ = self._run(tmp_path, HOSTILE_URLS, tool)
        assert result.returncode == 0, result.stderr
        assert not list(tmp_path.glob("PWNED*"))

    def test_each_url_reaches_the_tool_byte_for_byte(self, tmp_path, tool):
        _, calls = self._run(tmp_path, HOSTILE_URLS, tool)
        transfers = _transfers(calls)
        assert [c[-1] for c in transfers] == HOSTILE_URLS

    def test_destination_names_are_literal(self, tmp_path, tool):
        self._run(tmp_path, HOSTILE_URLS[:2], tool)
        names = sorted(p.name for p in (tmp_path / "fastq").iterdir())
        assert names == ["a$(touch PWNED_SUBST).fq.gz", "b`touch PWNED_TICK`.fq.gz"]

    def test_outdir_is_quoted_too(self, tmp_path, tool):
        outdir = "my dir $(touch PWNED_OUT)"
        result, _ = self._run(tmp_path, ["https://x/1.fq.gz"], tool, outdir=outdir)
        assert result.returncode == 0, result.stderr
        assert (tmp_path / outdir / "1.fq.gz").exists()
        assert not list(tmp_path.glob("PWNED*"))

    def test_script_parses(self, tmp_path, tool):
        script = tmp_path / "dl.sh"
        write_download_script([("s1", HOSTILE_URLS)], script, tool=tool)
        assert subprocess.run(["bash", "-n", str(script)]).returncode == 0

    def test_a_name_that_is_a_dot_segment_is_skipped(self, tmp_path, tool):
        """`https://x/..` would make the destination `$OUTDIR/..` itself."""
        body, n = build_download_script(
            [("s1", ["https://x/..", "https://x/.", "https://x/ok.fq.gz"])], tool=tool)
        assert n == 1
        assert "ok.fq.gz" in body


DATA_MD5 = "8d777f385d3dfec8815d20f7496026dc"  # md5 of b"data", the shim's default content


@pytest.mark.parametrize("tool", ["curl", "wget"])
class TestMd5Verification:
    """ENA publishes fastq_md5; `curl -C -` / `wget -c` can exit 0 on a truncated
    or stale file, so the checksum is the only proof a file is complete."""

    def _script(self, tmp_path, tool, md5):
        script = tmp_path / "dl.sh"
        write_download_script([("s1", ["https://x/1.fq.gz"])], script, tool=tool, md5=md5)
        return script

    def test_a_matching_file_passes(self, tmp_path, tool):
        script = self._script(tmp_path, tool, {"https://x/1.fq.gz": DATA_MD5})
        result, _ = run_with_shims(script, tmp_path)
        assert result.returncode == 0, result.stderr

    def test_a_truncated_file_fails_the_script(self, tmp_path, tool):
        script = self._script(tmp_path, tool, {"https://x/1.fq.gz": DATA_MD5})
        result, _ = run_with_shims(script, tmp_path, content="dat")
        assert result.returncode != 0
        assert "MD5 mismatch" in result.stderr
        assert "delete it and re-run" in result.stderr

    def test_a_malformed_checksum_is_never_embedded(self, tmp_path, tool):
        body, _ = build_download_script(
            [("s1", ["https://x/1.fq.gz"])], tool=tool,
            md5={"https://x/1.fq.gz": "$(touch PWNED)"})
        assert "PWNED" not in body
        assert "verify_md5 " not in body.replace("verify_md5() ", "")

    def test_only_files_with_a_checksum_are_verified(self, tmp_path, tool):
        body, _ = build_download_script(
            [("s1", ["https://x/1.fq.gz", "https://x/2.fq.gz"])], tool=tool,
            md5={"https://x/1.fq.gz": DATA_MD5.upper()})
        checks = [ln for ln in body.splitlines() if ln.startswith("verify_md5 ")]
        assert checks == [f'verify_md5 {DATA_MD5} "$OUTDIR"/1.fq.gz']

    def test_no_checksums_means_no_verification_block(self, tmp_path, tool):
        body, _ = build_download_script([("s1", ["https://x/1.fq.gz"])], tool=tool)
        assert "md5sum" not in body


class TestMd5Sidecar:
    def test_reads_url_tab_md5_rows(self, tmp_path):
        path = tmp_path / "fastq_md5.tsv"
        path.write_text(f"url\tmd5\nhttps://x/1.fq.gz\t{DATA_MD5}\n")
        assert read_md5_sidecar(path) == {"https://x/1.fq.gz": DATA_MD5}

    def test_drops_malformed_rows(self, tmp_path):
        path = tmp_path / "fastq_md5.tsv"
        path.write_text("url\tmd5\nhttps://x/1.fq.gz\tnot-a-hash\nonly-one-field\n")
        assert read_md5_sidecar(path) == {}

    def test_a_missing_sidecar_is_no_checksums(self, tmp_path):
        assert read_md5_sidecar(tmp_path / "absent.tsv") == {}
