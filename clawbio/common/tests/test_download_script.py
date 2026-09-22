"""Tests for the shared FASTQ download-script emitter.

Run with: pytest clawbio/common/tests/test_download_script.py -v

The emitter produces a bash script; it never downloads anything itself.
"""

from __future__ import annotations

import csv
import stat

import pytest

from clawbio.common.download_script import (
    SlurmOptions,
    build_download_script,
    urls_from_samplesheet,
    urls_from_url_list,
    write_download_script,
)


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
        """The URL is interpolated into double quotes in bash; a quote escapes them."""
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


class TestWrite:
    def test_written_script_is_executable(self, tmp_path):
        out = tmp_path / "download_fastq.sh"
        write_download_script([("s1", ["https://x/1.fq.gz"])], out, tool="wget")
        assert out.exists()
        assert out.stat().st_mode & stat.S_IXUSR

    def test_refuses_to_write_an_empty_script(self, tmp_path):
        with pytest.raises(SystemExit, match="No FASTQ URLs"):
            write_download_script([], tmp_path / "x.sh", tool="wget")
