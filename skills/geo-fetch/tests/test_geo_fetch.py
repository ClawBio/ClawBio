"""Tests for geo-fetch.

Run with: pytest skills/geo-fetch/tests/test_geo_fetch.py -v

No network required: the vendored client's single HTTP entry point replays a
recorded URL->body map captured from one real GSE30720 run.
"""

import csv
import json
import re
import subprocess
import sys
from pathlib import Path
from unittest.mock import patch

import pytest
from inline_snapshot import snapshot
from clawbio.common import demo_runner

SKILL_DIR = Path(__file__).resolve().parents[1]
SCRIPT = SKILL_DIR / "geo_fetch.py"
DEMO = "GSE30720"

sys.path.insert(0, str(SKILL_DIR))


class TestCredentialGate:
    """'Presence of an env var is not consent.' Upstream's rule, kept verbatim."""

    def test_credentials_are_not_sent_by_default(self, monkeypatch):
        import geo_fetch_api as api

        monkeypatch.setattr(api, "NCBI_EMAIL", "me@example.org")
        monkeypatch.setattr(api, "NCBI_API_KEY", "secret")
        monkeypatch.setattr(api, "_USE_CREDENTIALS", False)
        params = api._eutil_params({"db": "gds"})
        assert "email" not in params
        assert "api_key" not in params
        assert "tool" not in params

    def test_credentials_are_sent_only_when_opted_in(self, monkeypatch):
        import geo_fetch_api as api

        monkeypatch.setattr(api, "NCBI_EMAIL", "me@example.org")
        monkeypatch.setattr(api, "NCBI_API_KEY", "secret")
        monkeypatch.setattr(api, "_USE_CREDENTIALS", True)
        params = api._eutil_params({"db": "gds"})
        assert params["email"] == "me@example.org"
        assert params["api_key"] == "secret"
        assert params["tool"] == "geo-skill"

    def test_opting_in_with_no_env_vars_sends_nothing(self, monkeypatch):
        import geo_fetch_api as api

        monkeypatch.setattr(api, "NCBI_EMAIL", "")
        monkeypatch.setattr(api, "NCBI_API_KEY", "")
        monkeypatch.setattr(api, "_USE_CREDENTIALS", True)
        assert api._eutil_params({"db": "gds"}) == {"db": "gds"}

    def test_demo_reports_that_no_credentials_were_sent(self, tmp_path):
        import geo_fetch as app

        app.main(["--demo", "--output", str(tmp_path)])
        result = json.loads((tmp_path / "result.json").read_text())
        assert result["summary"]["ncbi_credentials_sent"] is False


TO_HTTPS_CASES = [
    ("ftp://ftp.sra.ebi.ac.uk/vol1/x.fastq.gz", "https://ftp.sra.ebi.ac.uk/vol1/x.fastq.gz"),
    ("https://ftp.sra.ebi.ac.uk/vol1/x.fastq.gz", "https://ftp.sra.ebi.ac.uk/vol1/x.fastq.gz"),
    ("http://example.org/x.fastq.gz", "http://example.org/x.fastq.gz"),
    # ENA's fastq_ftp is scheme-less host/path; that is why a scheme is added at all
    ("ftp.sra.ebi.ac.uk/vol1/x.fastq.gz", "https://ftp.sra.ebi.ac.uk/vol1/x.fastq.gz"),
    # not URLs: returned as None, never dressed up as https://<filename>
    ("sample_R1.fastq.gz", None),
    ("/data/fastq/sample_R1.fastq.gz", None),
    ("fastq/sample_R1.fastq.gz", None),
    ("", None),
]


class TestVendoredApi:
    @pytest.mark.parametrize("value, expected", TO_HTTPS_CASES)
    def test_to_https_only_accepts_urls(self, value, expected):
        import geo_fetch_api as api

        assert api.to_https(value) == expected

    def test_ftp_paths_shard_the_accession(self):
        """GEO's FTP layout buckets by accession: GSE30720 -> series/GSE30nnn."""
        import geo_fetch_api as api

        assert api.geo_ftp_dir("GSE30720") == ("series", "GSE30nnn")
        assert api.geo_ftp_dir("GSM12345") == ("samples", "GSM12nnn")
        assert api.geo_ftp_dir("GSE123") == ("series", "GSEnnn")

    def test_clean_val_replaces_control_characters(self):
        import geo_fetch_api as api

        assert api._clean_val("a\tb\nc") == "a_b_c"

    def test_metadata_tsv_is_lf_terminated(self, tmp_path):
        import geo_fetch_api as api

        out = tmp_path / "metadata.tsv"
        api.write_metadata_tsv([api.harmonize_row("GSM1", "1", {})], out)
        assert b"\r\n" not in out.read_bytes()

    def test_read_map_rejects_malformed_values(self):
        import geo_fetch_api as api

        for bad in ("0,1", "x,2", "1", "1,2,3"):
            with pytest.raises(SystemExit):
                api.parse_read_map(bad)


class TestDemo:
    def test_demo_makes_no_network_call(self, tmp_path):
        import geo_fetch as app
        import geo_fetch_api as api

        with patch.object(api.urllib.request, "urlopen") as mocked:
            app.main(["--demo", "--output", str(tmp_path)])
            mocked.assert_not_called()

    def test_demo_writes_report_and_result(self, tmp_path):
        import geo_fetch as app

        app.main(["--demo", "--output", str(tmp_path)])
        assert (tmp_path / "report.md").exists()
        assert json.loads((tmp_path / "result.json").read_text())["skill"] == "geo-fetch"

    def test_demo_report_carries_the_disclaimer(self, tmp_path):
        import geo_fetch as app

        app.main(["--demo", "--output", str(tmp_path)])
        text = (tmp_path / "report.md").read_text()
        assert "research and educational tool" in text
        assert "not a medical device" in text

    def test_demo_writes_the_reproducibility_bundle(self, tmp_path):
        import geo_fetch as app

        app.main(["--demo", "--output", str(tmp_path)])
        for name in ("commands.sh", "environment.yml", "checksums.sha256"):
            assert (tmp_path / "reproducibility" / name).exists()

    def test_checksum_labels_resolve_from_the_output_dir(self, tmp_path):
        import geo_fetch as app

        app.main(["--demo", "--output", str(tmp_path)])
        for line in (tmp_path / "reproducibility" / "checksums.sha256").read_text().splitlines():
            assert (tmp_path / line.split("  ", 1)[1]).exists()

    def test_demo_metadata_table_has_one_row_per_sample(self, tmp_path):
        import geo_fetch as app

        app.main(["--demo", "--output", str(tmp_path)])
        rows = list(csv.DictReader((tmp_path / "tables" / "metadata.tsv").open(), delimiter="\t"))
        assert len(rows) == 42
        assert {r["species"] for r in rows} == {"Arabidopsis thaliana"}

    def test_demo_writes_a_pipeline_ready_samplesheet(self, tmp_path):
        """GEO resolves to its SRA project, then to ENA for the FASTQ links."""
        import geo_fetch as app

        app.main(["--demo", "--output", str(tmp_path)])
        rows = list(csv.DictReader((tmp_path / "samplesheet.csv").open()))
        assert list(rows[0]) == ["sample", "fastq_1", "fastq_2", "strandedness"]
        assert len(rows) == 42
        assert all(r["fastq_1"].startswith("http") for r in rows)

    def test_demo_download_script_downloads_nothing(self, tmp_path):
        import geo_fetch as app

        app.main(["--demo", "--output", str(tmp_path)])
        body = (tmp_path / "download_geo.sh").read_text()
        assert "set -euo pipefail" in body
        assert "# #SBATCH --partition=<your_partition>" in body
        assert not (tmp_path / "fastq").exists()


class TestCLI:
    def test_no_args_exits_nonzero(self):
        result = subprocess.run([sys.executable, str(SCRIPT)], capture_output=True, text=True)
        assert result.returncode != 0

    def test_command_flag_dispatches(self, tmp_path):
        import geo_fetch as app

        app.main(["--demo", "--command", "metadata", "--output", str(tmp_path)])
        assert "## metadata" in (tmp_path / "report.md").read_text()

    def test_upstream_positional_form_still_works(self, tmp_path):
        import geo_fetch as app

        app._install_demo_transport()
        app.main(["metadata", DEMO, "--output", str(tmp_path)])
        assert "## metadata" in (tmp_path / "report.md").read_text()

    def test_samplesheet_requires_an_assay(self, tmp_path):
        import geo_fetch as app

        with pytest.raises(SystemExit):
            app.main(["--command", "samplesheet", "--accession", DEMO,
                      "--output", str(tmp_path)])


class _ChunkedResponse:
    """Serves a body in pieces and refuses a whole-body read()."""

    def __init__(self, body, chunk=4, headers=None):
        self._body, self._chunk, self._pos = body, chunk, 0
        self.headers = headers or {}

    def read(self, n=-1):
        if n is None or n < 0:
            raise AssertionError("download() must stream, not read() the whole body")
        piece = self._body[self._pos:self._pos + n]
        self._pos += len(piece)
        return piece

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


class TestDownloadStreams:
    """Supplementary files run to many GB; the other four skills already stream."""

    def test_download_streams_to_disk(self, tmp_path):
        import geo_fetch_api as api

        body = b"x" * 10_000
        with patch.object(api.urllib.request, "urlopen", return_value=_ChunkedResponse(body)):
            dest = api.download("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE1nnn/GSE1/suppl/a.tar",
                                str(tmp_path))
        assert Path(dest).read_bytes() == body
        assert not list(tmp_path.glob("*.part"))

    def test_download_retries_a_failed_attempt(self, tmp_path):
        import geo_fetch_api as api

        responses = [OSError("reset"), _ChunkedResponse(b"ok")]

        def flaky(*a, **k):
            r = responses.pop(0)
            if isinstance(r, Exception):
                raise r
            return r

        with patch.object(api.urllib.request, "urlopen", side_effect=flaky), \
                patch.object(api.time, "sleep"):
            dest = api.download("https://ftp.ncbi.nlm.nih.gov/x/b.txt", str(tmp_path))
        assert Path(dest).read_bytes() == b"ok"

    def test_download_refuses_a_dot_dot_name(self, tmp_path):
        import geo_fetch_api as api

        with patch.object(api.urllib.request, "urlopen",
                          side_effect=AssertionError("must refuse before any request")):
            with pytest.raises(SystemExit, match="outside"):
                api.download("https://ftp.ncbi.nlm.nih.gov/x/..", str(tmp_path / "out"))


class TestDownloadSizeCheck:
    """GEO publishes no checksums (no md5 files, no ETag/Content-MD5), so size is
    the only verifiable property: Content-Length per file, and filelist.txt for
    the _RAW.tar. It catches truncation, not same-length corruption."""

    URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE1nnn/GSE1/suppl/GSE1_RAW.tar"

    def _serve(self, *responses):
        responses = list(responses)
        calls = []

        def urlopen(*a, **k):
            calls.append(a)
            return responses.pop(0)

        return urlopen, calls

    def test_a_short_body_is_retried_then_refused(self, tmp_path):
        import geo_fetch_api as api

        short = [_ChunkedResponse(b"12345", headers={"Content-Length": "10"}) for _ in range(3)]
        urlopen, calls = self._serve(*short)
        with patch.object(api.urllib.request, "urlopen", side_effect=urlopen), \
                patch.object(api.time, "sleep"):
            with pytest.raises(SystemExit, match="size"):
                api.download(self.URL, str(tmp_path))
        assert len(calls) == 3
        assert list(tmp_path.iterdir()) == []

    def test_a_short_first_attempt_is_recovered_by_the_retry(self, tmp_path):
        import geo_fetch_api as api

        urlopen, _ = self._serve(
            _ChunkedResponse(b"12345", headers={"Content-Length": "10"}),
            _ChunkedResponse(b"1234567890", headers={"Content-Length": "10"}))
        with patch.object(api.urllib.request, "urlopen", side_effect=urlopen), \
                patch.object(api.time, "sleep"):
            dest = api.download(self.URL, str(tmp_path))
        assert Path(dest).read_bytes() == b"1234567890"

    def test_a_matching_length_is_reported_as_verified(self, tmp_path, capsys):
        import geo_fetch_api as api

        urlopen, _ = self._serve(_ChunkedResponse(b"abc", headers={"Content-Length": "3"}))
        with patch.object(api.urllib.request, "urlopen", side_effect=urlopen):
            api.download(self.URL, str(tmp_path))
        assert "size verified" in capsys.readouterr().err

    def test_no_length_header_still_saves_but_says_so(self, tmp_path, capsys):
        import geo_fetch_api as api

        urlopen, _ = self._serve(_ChunkedResponse(b"abc"))
        with patch.object(api.urllib.request, "urlopen", side_effect=urlopen):
            dest = api.download(self.URL, str(tmp_path))
        assert Path(dest).read_bytes() == b"abc"
        err = capsys.readouterr().err
        assert "could not be verified" in err and "size verified" not in err

    def test_a_filelist_mismatch_fails_at_once(self, tmp_path):
        """Content-Length agreed with the body, so the archive itself disagrees
        with its filelist.txt; fetching a multi-GB tar again would not help."""
        import geo_fetch_api as api

        urlopen, calls = self._serve(_ChunkedResponse(b"abc", headers={"Content-Length": "3"}))
        with patch.object(api.urllib.request, "urlopen", side_effect=urlopen), \
                patch.object(api.time, "sleep"):
            with pytest.raises(SystemExit, match="filelist.txt"):
                api.download(self.URL, str(tmp_path), expected_size=5)
        assert len(calls) == 1
        assert list(tmp_path.iterdir()) == []

    def test_cmd_download_passes_the_filelist_size_for_the_raw_tar(self, tmp_path):
        import argparse

        import geo_fetch_api as api

        filelist = ("#Archive/File\tName\tTime\tSize\tType\n"
                    "Archive\tGSE30720_RAW.tar\t01/17/2013 15:25:56\t11786240\tTAR\n"
                    "File\tGSM762070_x.txt.gz\t07/15/2011 15:46:13\t123964\tTXT\n")
        seen = {}
        args = argparse.Namespace(accession="GSE30720", matrix=False, soft=False,
                                  miniml=False, suppl=True, out=str(tmp_path))
        with patch.object(api, "list_ftp_files",
                          return_value=["GSE30720_RAW.tar", "filelist.txt", "GSE30720_counts.txt.gz"]), \
                patch.object(api, "http_get", return_value=filelist.encode()), \
                patch.object(api, "download",
                             side_effect=lambda url, out, expected_size=None:
                             seen.__setitem__(url.rsplit("/", 1)[-1], expected_size)):
            api.cmd_download(args)
        assert seen == {"GSE30720_RAW.tar": 11786240, "filelist.txt": None,
                        "GSE30720_counts.txt.gz": None}

    def test_a_malformed_filelist_only_warns(self, tmp_path, capsys):
        import argparse

        import geo_fetch_api as api

        seen = {}
        args = argparse.Namespace(accession="GSE30720", matrix=False, soft=False,
                                  miniml=False, suppl=True, out=str(tmp_path))
        with patch.object(api, "list_ftp_files", return_value=["GSE30720_RAW.tar", "filelist.txt"]), \
                patch.object(api, "http_get", return_value=b"Archive\tGSE30720_RAW.tar\tnot-a-size\n"), \
                patch.object(api, "download",
                             side_effect=lambda url, out, expected_size=None:
                             seen.__setitem__(url.rsplit("/", 1)[-1], expected_size)):
            api.cmd_download(args)
        assert seen["GSE30720_RAW.tar"] is None
        assert "filelist.txt" in capsys.readouterr().err


class TestSafety:
    def test_warns_before_overwriting(self, tmp_path, capsys):
        import geo_fetch as app

        app.main(["--demo", "--output", str(tmp_path)])
        capsys.readouterr()
        app.main(["--demo", "--output", str(tmp_path)])
        assert "overwritten" in capsys.readouterr().err

    def test_demo_writes_nothing_outside_the_output_dir(self, tmp_path, monkeypatch):
        """Upstream's `runtable --out .` defaulted to the working directory."""
        import geo_fetch as app

        cwd = tmp_path / "cwd"
        cwd.mkdir()
        monkeypatch.chdir(cwd)
        app.main(["--demo", "--output", str(tmp_path / "out")])
        assert list(cwd.iterdir()) == []

    def test_report_holds_no_absolute_output_path(self, tmp_path):
        import geo_fetch as app

        app.main(["--demo", "--output", str(tmp_path)])
        assert str(tmp_path) not in (tmp_path / "report.md").read_text()


def _parse_output_contract(skill_md):
    """Extract files promised in the SKILL.md '## Output Structure' tree."""
    if not skill_md.exists():
        return []
    text = skill_md.read_text()
    m = re.search(r"##\s*Output Structure\s*\n+```[^\n]*\n(.*?)\n```", text, re.S)
    if not m:
        return []
    files, parents = [], {}
    for raw in m.group(1).splitlines():
        if not raw.strip():
            continue
        parts = re.split(r"\s+#", raw, maxsplit=1)
        entry, comment = parts[0], (parts[1] if len(parts) > 1 else "")
        mm = re.match(r"^([\s│├└─]*)(.*)$", entry)
        prefix, name = mm.group(1), mm.group(2).strip()
        if not name:
            continue
        depth = len(prefix) // 4
        if depth == 0:
            continue
        if name.endswith("/"):
            parents[depth] = name.rstrip("/")
            for d in [k for k in parents if k > depth]:
                del parents[d]
            continue
        if "optional" in comment.lower():
            continue
        rel = "/".join(parents[d] for d in sorted(parents) if d < depth)
        files.append(rel + "/" + name if rel else name)
    return files


class TestOutputContract:
    """Every artifact promised in SKILL.md '## Output Structure' must be produced."""

    def test_documented_outputs_are_produced(self, tmp_path):
        promised = _parse_output_contract(SKILL_DIR / "SKILL.md")
        if not promised:
            pytest.skip("No parseable '## Output Structure' section in SKILL.md")
        result = subprocess.run(
            [sys.executable, str(SCRIPT), "--demo", "--output", str(tmp_path)],
            capture_output=True, text=True)
        assert result.returncode == 0, f"demo run failed: {result.stderr}"
        missing = [p for p in promised if not (tmp_path / p).exists()]
        assert not missing, (
            "SKILL.md Output Structure promises artifacts the skill did not "
            "produce: " + ", ".join(missing))


# ── Demo output snapshot ──────────────────────────────────────────────────────
# Fails when a code change alters any demo number. If the change is intended,
# run `pytest skills/geo-fetch/tests --inline-snapshot=fix`; it rewrites the
# snapshot() values below, and the PR diff shows which numbers moved.

@pytest.mark.demo_snapshot
def test_demo_output_matches_snapshot(tmp_path):
    assert demo_runner.demo_outputs("geo-fetch", tmp_path) == snapshot({'download_geo.sh': """\
#!/bin/bash
#
#SBATCH --job-name=fastq_download
# #SBATCH --partition=<your_partition>   # set for your cluster
# #SBATCH --account=<your_account>        # set if required
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=24:00:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
# #SBATCH --mail-type=END,FAIL
# #SBATCH --mail-user=<you@example.org>

set -euo pipefail

OUTDIR=fastq
mkdir -p "$OUTDIR"

# curl >= 7.71 retries HTTP errors (403/429/503) too, not just network
# failures. Older curl rejects the flag outright, so probe before using it.
RETRY_ALL=""
if curl --help all 2>/dev/null | grep -q -- --retry-all-errors; then
  RETRY_ALL="--retry-all-errors"
fi

# sample: SAMN00668211
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/SRR309145.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309145/SRR309145.fastq.gz

# sample: SAMN00668212
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/SRR309146.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309146/SRR309146.fastq.gz

# sample: SAMN00668213
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/SRR309147.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309147/SRR309147.fastq.gz

# sample: SAMN00668214
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/SRR309148.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309148/SRR309148.fastq.gz

# sample: SAMN00668215
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/SRR309149.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309149/SRR309149.fastq.gz

# sample: SAMN00668216
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/SRR309150.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309150/SRR309150.fastq.gz

# sample: SAMN00668217
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/SRR309151.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309151/SRR309151.fastq.gz

# sample: SAMN00668218
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/SRR309152.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309152/SRR309152.fastq.gz

# sample: SAMN00668219
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/SRR309153.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309153/SRR309153.fastq.gz

# sample: SAMN00668220
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/SRR309154.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309154/SRR309154.fastq.gz

# sample: SAMN00668221
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/SRR309155.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309155/SRR309155.fastq.gz

# sample: SAMN00668222
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/SRR309156.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309156/SRR309156.fastq.gz

# sample: SAMN00668223
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/SRR309157.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309157/SRR309157.fastq.gz

# sample: SAMN00668224
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/SRR309158.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309158/SRR309158.fastq.gz

# sample: SAMN00668225
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/SRR309159.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309159/SRR309159.fastq.gz

# sample: SAMN00668226
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/SRR309160.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309160/SRR309160.fastq.gz

# sample: SAMN00668227
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/SRR309161.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309161/SRR309161.fastq.gz

# sample: SAMN00668228
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/SRR309162.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309162/SRR309162.fastq.gz

# sample: SAMN00668229
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/SRR309163.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309163/SRR309163.fastq.gz

# sample: SAMN00668230
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/SRR309164.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309164/SRR309164.fastq.gz

# sample: SAMN00668231
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/SRR309165.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309165/SRR309165.fastq.gz

# sample: SAMN00668232
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/SRR309166.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309166/SRR309166.fastq.gz

# sample: SAMN00668233
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/SRR309167.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309167/SRR309167.fastq.gz

# sample: SAMN00668234
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/SRR309168.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309168/SRR309168.fastq.gz

# sample: SAMN00668235
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/SRR309169.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309169/SRR309169.fastq.gz

# sample: SAMN00668236
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/SRR309170.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309170/SRR309170.fastq.gz

# sample: SAMN00668237
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/SRR309171.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309171/SRR309171.fastq.gz

# sample: SAMN00668238
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/SRR309172.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309172/SRR309172.fastq.gz

# sample: SAMN00668239
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/SRR309173.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309173/SRR309173.fastq.gz

# sample: SAMN00668240
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/SRR309174.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309174/SRR309174.fastq.gz

# sample: SAMN00668241
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/SRR309175.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309175/SRR309175.fastq.gz

# sample: SAMN00668242
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/SRR309176.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309176/SRR309176.fastq.gz

# sample: SAMN00668243
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/SRR309177.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309177/SRR309177.fastq.gz

# sample: SAMN00668244
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/SRR309178.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309178/SRR309178.fastq.gz

# sample: SAMN00668245
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/SRR309179.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309179/SRR309179.fastq.gz

# sample: SAMN00668246
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/SRR309180.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309180/SRR309180.fastq.gz

# sample: SAMN00668247
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/SRR309181.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309181/SRR309181.fastq.gz

# sample: SAMN00668248
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/SRR309182.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309182/SRR309182.fastq.gz

# sample: SAMN00668249
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/SRR309183.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309183/SRR309183.fastq.gz

# sample: SAMN00668250
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/SRR309184.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309184/SRR309184.fastq.gz

# sample: SAMN00668251
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/SRR309185.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309185/SRR309185.fastq.gz

# sample: SAMN00668252
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/SRR309186.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309186/SRR309186.fastq.gz

echo "Downloaded 42 file(s) to $OUTDIR"
""", 'report.md': """\
# GEO report — GSE30720

Source: [Gene Expression Omnibus](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE30720)

## metadata

```
accession   : GSE30720
title       : Seedling transcriptome sequencing of the Arabidopsis thaliana MAGIC founder accessions
taxon       : Arabidopsis thaliana
gdstype     : Expression profiling by high throughput sequencing
gpl         : 11221
n_samples   : 42
pdat        : 2011/08/29
summary     : We have sequenced messenger RNA isolated from seedling tissue for 19 accessions of Arabidopsis thaliana (with biological replication).  The 19 accessions for which RNA-Seq reads were collected have served as the founders for the MAGIC lines, a high-resolution recombinant inbred line mapping resource.  RNA sequencing data was used to examine differential gene expression among the accessions.
```

## samples

```
GSM762080	ct_1_seedling_br1
GSM762100	sf_2_seedling_br1
GSM762097	po_0_seedling_br2
GSM762103	tsu_0_seedling_br2
GSM762077	can_0_seedling_br2
GSM762083	edi_0_seedling_br2
GSM762109	wu_0_seedling_br2
GSM762086	kn_0_seedling_br1
GSM762089	ler_0_seedling_br2
GSM762108	wu_0_seedling_br1
GSM762106	ws_0_seedling_br1
GSM762082	edi_0_seedling_br1
GSM762102	tsu_0_seedling_br1
GSM762105	wil_2_seedling_br2
GSM762071	col_0_seedling_tr2
GSM762085	hi_0_seedling_br2
GSM762091	mt_0_seedling_br2
GSM762094	oy_0_seedling_br1
GSM762111	zu_0_seedling_br2
GSM762110	zu_0_seedling_br1
GSM762074	bur_0_seedling_br1
GSM762107	ws_0_seedling_br2
GSM762088	ler_0_seedling_br1
GSM762090	mt_0_seedling_br1
GSM762070	col_0_seedling_tr1
GSM762099	rsch_4_seedling_br2
GSM762073	col_0_seedling_tr4
GSM762076	can_0_seedling_br1
GSM762087	kn_0_seedling_br2
GSM762081	ct_1_seedling_br2
GSM762079	col_0_seedling_br2
GSM762084	hi_0_seedling_br1
GSM762092	no_0_seedling_br1
GSM762093	no_0_seedling_br2
GSM762072	col_0_seedling_tr3
GSM762078	col_0_seedling_br1
GSM762075	bur_0_seedling_br2
GSM762096	po_0_seedling_br1
GSM762095	oy_0_seedling_br2
GSM762101	sf_2_seedling_br2
GSM762104	wil_2_seedling_br1
GSM762098	rsch_4_seedling_br1
```

## files

```
[matrix]
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE30nnn/GSE30720/matrix/GSE30720_series_matrix.txt.gz
[soft]
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE30nnn/GSE30720/soft/GSE30720_family.soft.gz
[miniml]
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE30nnn/GSE30720/miniml/GSE30720_family.xml.tgz
[suppl]
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE30nnn/GSE30720/suppl/GSE30720_RAW.tar
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE30nnn/GSE30720/suppl/filelist.txt
```

## metadata-table

```
Wrote 42 sample row(s) x 10 column(s) for GSE30720 to <output>/tables/metadata.tsv
```

## samplesheet

```
GSE30720 -> SRA project PRJNA154687
Wrote 42 row(s) for 42 sample(s) to <output>/samplesheet.csv (nf-core/rnaseq)
```

## download-script

```
Wrote download_geo.sh: 42 download command(s) using curl.
Nothing has been downloaded. Run it with --run, submit it with --submit, or execute it yourself.
```

---

*ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions.*
""", 'result.json': """\
{
  "status": "ok",
  "ok": true,
  "skill": "geo-fetch",
  "version": "0.1.0",
  "completed_at": "<TIMESTAMP>",
  "input_checksum": "",
  "datasets": {},
  "summary": {
    "accession": "GSE30720",
    "commands": [
      "metadata",
      "samples",
      "files",
      "metadata-table",
      "samplesheet",
      "download-script"
    ],
    "demo": true,
    "ncbi_credentials_sent": false
  },
  "data": {
    "sections": {
      "metadata": "accession   : GSE30720\\ntitle       : Seedling transcriptome sequencing of the Arabidopsis thaliana MAGIC founder accessions\\ntaxon       : Arabidopsis thaliana\\ngdstype     : Expression profiling by high throughput sequencing\\ngpl         : 11221\\nn_samples   : 42\\npdat        : 2011/08/29\\nsummary     : We have sequenced messenger RNA isolated from seedling tissue for 19 accessions of Arabidopsis thaliana (with biological replication).  The 19 accessions for which RNA-Seq reads were collected have served as the founders for the MAGIC lines, a high-resolution recombinant inbred line mapping resource.  RNA sequencing data was used to examine differential gene expression among the accessions.\\n",
      "samples": "GSM762080\\tct_1_seedling_br1\\nGSM762100\\tsf_2_seedling_br1\\nGSM762097\\tpo_0_seedling_br2\\nGSM762103\\ttsu_0_seedling_br2\\nGSM762077\\tcan_0_seedling_br2\\nGSM762083\\tedi_0_seedling_br2\\nGSM762109\\twu_0_seedling_br2\\nGSM762086\\tkn_0_seedling_br1\\nGSM762089\\tler_0_seedling_br2\\nGSM762108\\twu_0_seedling_br1\\nGSM762106\\tws_0_seedling_br1\\nGSM762082\\tedi_0_seedling_br1\\nGSM762102\\ttsu_0_seedling_br1\\nGSM762105\\twil_2_seedling_br2\\nGSM762071\\tcol_0_seedling_tr2\\nGSM762085\\thi_0_seedling_br2\\nGSM762091\\tmt_0_seedling_br2\\nGSM762094\\toy_0_seedling_br1\\nGSM762111\\tzu_0_seedling_br2\\nGSM762110\\tzu_0_seedling_br1\\nGSM762074\\tbur_0_seedling_br1\\nGSM762107\\tws_0_seedling_br2\\nGSM762088\\tler_0_seedling_br1\\nGSM762090\\tmt_0_seedling_br1\\nGSM762070\\tcol_0_seedling_tr1\\nGSM762099\\trsch_4_seedling_br2\\nGSM762073\\tcol_0_seedling_tr4\\nGSM762076\\tcan_0_seedling_br1\\nGSM762087\\tkn_0_seedling_br2\\nGSM762081\\tct_1_seedling_br2\\nGSM762079\\tcol_0_seedling_br2\\nGSM762084\\thi_0_seedling_br1\\nGSM762092\\tno_0_seedling_br1\\nGSM762093\\tno_0_seedling_br2\\nGSM762072\\tcol_0_seedling_tr3\\nGSM762078\\tcol_0_seedling_br1\\nGSM762075\\tbur_0_seedling_br2\\nGSM762096\\tpo_0_seedling_br1\\nGSM762095\\toy_0_seedling_br2\\nGSM762101\\tsf_2_seedling_br2\\nGSM762104\\twil_2_seedling_br1\\nGSM762098\\trsch_4_seedling_br1\\n",
      "files": "[matrix]\\n  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE30nnn/GSE30720/matrix/GSE30720_series_matrix.txt.gz\\n[soft]\\n  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE30nnn/GSE30720/soft/GSE30720_family.soft.gz\\n[miniml]\\n  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE30nnn/GSE30720/miniml/GSE30720_family.xml.tgz\\n[suppl]\\n  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE30nnn/GSE30720/suppl/GSE30720_RAW.tar\\n  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE30nnn/GSE30720/suppl/filelist.txt\\n",
      "metadata-table": "Wrote 42 sample row(s) x 10 column(s) for GSE30720 to <output>/tables/metadata.tsv\\n",
      "samplesheet": "GSE30720 -> SRA project PRJNA154687\\nWrote 42 row(s) for 42 sample(s) to <output>/samplesheet.csv (nf-core/rnaseq)\\n",
      "download-script": "Wrote download_geo.sh: 42 download command(s) using curl.\\nNothing has been downloaded. Run it with --run, submit it with --submit, or execute it yourself."
    }
  }
}\
""", 'samplesheet.csv': """\
sample,fastq_1,fastq_2,strandedness
SAMN00668211,https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309145/SRR309145.fastq.gz,,auto
SAMN00668212,https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309146/SRR309146.fastq.gz,,auto
SAMN00668213,https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309147/SRR309147.fastq.gz,,auto
SAMN00668214,https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309148/SRR309148.fastq.gz,,auto
SAMN00668215,https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309149/SRR309149.fastq.gz,,auto
SAMN00668216,https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309150/SRR309150.fastq.gz,,auto
SAMN00668217,https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309151/SRR309151.fastq.gz,,auto
SAMN00668218,https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309152/SRR309152.fastq.gz,,auto
SAMN00668219,https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309153/SRR309153.fastq.gz,,auto
SAMN00668220,https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309154/SRR309154.fastq.gz,,auto
SAMN00668221,https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309155/SRR309155.fastq.gz,,auto
SAMN00668222,https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309156/SRR309156.fastq.gz,,auto
SAMN00668223,https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309157/SRR309157.fastq.gz,,auto
SAMN00668224,https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309158/SRR309158.fastq.gz,,auto
SAMN00668225,https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309159/SRR309159.fastq.gz,,auto
SAMN00668226,https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309160/SRR309160.fastq.gz,,auto
SAMN00668227,https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309161/SRR309161.fastq.gz,,auto
SAMN00668228,https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309162/SRR309162.fastq.gz,,auto
SAMN00668229,https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309163/SRR309163.fastq.gz,,auto
SAMN00668230,https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309164/SRR309164.fastq.gz,,auto
SAMN00668231,https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309165/SRR309165.fastq.gz,,auto
SAMN00668232,https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309166/SRR309166.fastq.gz,,auto
SAMN00668233,https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309167/SRR309167.fastq.gz,,auto
SAMN00668234,https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309168/SRR309168.fastq.gz,,auto
SAMN00668235,https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309169/SRR309169.fastq.gz,,auto
SAMN00668236,https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309170/SRR309170.fastq.gz,,auto
SAMN00668237,https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309171/SRR309171.fastq.gz,,auto
SAMN00668238,https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309172/SRR309172.fastq.gz,,auto
SAMN00668239,https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309173/SRR309173.fastq.gz,,auto
SAMN00668240,https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309174/SRR309174.fastq.gz,,auto
SAMN00668241,https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309175/SRR309175.fastq.gz,,auto
SAMN00668242,https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309176/SRR309176.fastq.gz,,auto
SAMN00668243,https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309177/SRR309177.fastq.gz,,auto
SAMN00668244,https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309178/SRR309178.fastq.gz,,auto
SAMN00668245,https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309179/SRR309179.fastq.gz,,auto
SAMN00668246,https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309180/SRR309180.fastq.gz,,auto
SAMN00668247,https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309181/SRR309181.fastq.gz,,auto
SAMN00668248,https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309182/SRR309182.fastq.gz,,auto
SAMN00668249,https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309183/SRR309183.fastq.gz,,auto
SAMN00668250,https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309184/SRR309184.fastq.gz,,auto
SAMN00668251,https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309185/SRR309185.fastq.gz,,auto
SAMN00668252,https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR309/SRR309186/SRR309186.fastq.gz,,auto
""", 'tables/metadata.tsv': """\
sample	replicate	species	sex	age	condition	genotype	treatment	tissue	accession
GSM762080	biological replicate 1	Arabidopsis thaliana	NA	11 days	NA	NA	NA	seedling	Ct-1
GSM762100	biological replicate 1	Arabidopsis thaliana	NA	11 days	NA	NA	NA	seedling	Sf-2
GSM762097	biological replicate 2	Arabidopsis thaliana	NA	11 days	NA	NA	NA	seedling	Po-0
GSM762103	biological replicate 2	Arabidopsis thaliana	NA	11 days	NA	NA	NA	seedling	Tsu-0
GSM762077	biological replicate 2	Arabidopsis thaliana	NA	11 days	NA	NA	NA	seedling	Can-0
GSM762083	biological replicate 2	Arabidopsis thaliana	NA	11 days	NA	NA	NA	seedling	Edi-0
GSM762109	biological replicate 2	Arabidopsis thaliana	NA	11 days	NA	NA	NA	seedling	Wu-0
GSM762086	biological replicate 1	Arabidopsis thaliana	NA	11 days	NA	NA	NA	seedling	Kn-0
GSM762089	biological replicate 2	Arabidopsis thaliana	NA	11 days	NA	NA	NA	seedling	Ler-0
GSM762108	biological replicate 1	Arabidopsis thaliana	NA	11 days	NA	NA	NA	seedling	Wu-0
GSM762106	biological replicate 1	Arabidopsis thaliana	NA	11 days	NA	NA	NA	seedling	Ws-0
GSM762082	biological replicate 1	Arabidopsis thaliana	NA	11 days	NA	NA	NA	seedling	Edi-0
GSM762102	biological replicate 1	Arabidopsis thaliana	NA	11 days	NA	NA	NA	seedling	Tsu-0
GSM762105	biological replicate 2	Arabidopsis thaliana	NA	11 days	NA	NA	NA	seedling	Wil-2
GSM762071	technical replicate 2	Arabidopsis thaliana	NA	11 days	NA	NA	NA	seedling	Col-0
GSM762085	biological replicate 2	Arabidopsis thaliana	NA	11 days	NA	NA	NA	seedling	Hi-0
GSM762091	biological replicate 2	Arabidopsis thaliana	NA	11 days	NA	NA	NA	seedling	Mt-0
GSM762094	biological replicate 1	Arabidopsis thaliana	NA	11 days	NA	NA	NA	seedling	Oy-0
GSM762111	biological replicate 2	Arabidopsis thaliana	NA	11 days	NA	NA	NA	seedling	Zu-0
GSM762110	biological replicate 1	Arabidopsis thaliana	NA	11 days	NA	NA	NA	seedling	Zu-0
GSM762074	biological replicate 1	Arabidopsis thaliana	NA	11 days	NA	NA	NA	seedling	Bur-0
GSM762107	biological replicate 2	Arabidopsis thaliana	NA	11 days	NA	NA	NA	seedling	Ws-0
GSM762088	biological replicate 1	Arabidopsis thaliana	NA	11 days	NA	NA	NA	seedling	Ler-0
GSM762090	biological replicate 1	Arabidopsis thaliana	NA	11 days	NA	NA	NA	seedling	Mt-0
GSM762070	technical replicate 1	Arabidopsis thaliana	NA	11 days	NA	NA	NA	seedling	Col-0
GSM762099	biological replicate 2	Arabidopsis thaliana	NA	11 days	NA	NA	NA	seedling	Rsch-4
GSM762073	technical replicate 4	Arabidopsis thaliana	NA	11 days	NA	NA	NA	seedling	Col-0
GSM762076	biological replicate 1	Arabidopsis thaliana	NA	11 days	NA	NA	NA	seedling	Can-0
GSM762087	biological replicate 2	Arabidopsis thaliana	NA	11 days	NA	NA	NA	seedling	Kn-0
GSM762081	biological replicate 2	Arabidopsis thaliana	NA	11 days	NA	NA	NA	seedling	Ct-1
GSM762079	biological replicate 2	Arabidopsis thaliana	NA	11 days	NA	NA	NA	seedling	Col-0
GSM762084	biological replicate 1	Arabidopsis thaliana	NA	11 days	NA	NA	NA	seedling	Hi-0
GSM762092	biological replicate 1	Arabidopsis thaliana	NA	11 days	NA	NA	NA	seedling	No-0
GSM762093	biological replicate 2	Arabidopsis thaliana	NA	11 days	NA	NA	NA	seedling	No-0
GSM762072	technical replicate 3	Arabidopsis thaliana	NA	11 days	NA	NA	NA	seedling	Col-0
GSM762078	biological replicate 1	Arabidopsis thaliana	NA	11 days	NA	NA	NA	seedling	Col-0
GSM762075	biological replicate 2	Arabidopsis thaliana	NA	11 days	NA	NA	NA	seedling	Bur-0
GSM762096	biological replicate 1	Arabidopsis thaliana	NA	11 days	NA	NA	NA	seedling	Po-0
GSM762095	biological replicate 2	Arabidopsis thaliana	NA	11 days	NA	NA	NA	seedling	Oy-0
GSM762101	biological replicate 2	Arabidopsis thaliana	NA	11 days	NA	NA	NA	seedling	Sf-2
GSM762104	biological replicate 1	Arabidopsis thaliana	NA	11 days	NA	NA	NA	seedling	Wil-2
GSM762098	biological replicate 1	Arabidopsis thaliana	NA	11 days	NA	NA	NA	seedling	Rsch-4
"""})
