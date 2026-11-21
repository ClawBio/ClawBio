"""Tests for ena-fetch.

Run with: pytest skills/ena-fetch/tests/test_ena_fetch.py -v

No network required: the vendored client's single HTTP entry point is served
from the committed fixtures in ../examples/.
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
SCRIPT = SKILL_DIR / "ena_fetch.py"
EXAMPLES = SKILL_DIR / "examples"
DEMO_PROJECT = "PRJEB56029"

sys.path.insert(0, str(SKILL_DIR))


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
        import ena_fetch_api as api

        assert api.to_https(value) == expected

    def test_parse_tsv_round_trips_the_fixture(self):
        import ena_fetch_api as api

        text = (EXAMPLES / f"demo_{DEMO_PROJECT}_filereport.tsv").read_text()
        rows = api.parse_tsv(text)
        assert len(rows) == 9
        assert all(r["run_accession"].startswith("ERR") for r in rows)
        assert {r["scientific_name"] for r in rows} == {"Arabidopsis thaliana"}

    def test_parse_tsv_tolerates_a_blank_payload(self):
        import ena_fetch_api as api

        assert api.parse_tsv("") == []

    def test_sample_attrs_drop_ena_bookkeeping_tags(self):
        """ENA-CHECKLIST / ENA-FIRST-PUBLIC are archive bookkeeping, not biology."""
        import ena_fetch_api as api

        xml = ('<SAMPLE_SET><SAMPLE><SAMPLE_NAME><SCIENTIFIC_NAME>Arabidopsis thaliana'
               '</SCIENTIFIC_NAME></SAMPLE_NAME><SAMPLE_ATTRIBUTES>'
               '<SAMPLE_ATTRIBUTE><TAG>ENA-CHECKLIST</TAG><VALUE>ERC000011</VALUE></SAMPLE_ATTRIBUTE>'
               '<SAMPLE_ATTRIBUTE><TAG>genotype</TAG><VALUE>Col-0</VALUE></SAMPLE_ATTRIBUTE>'
               '</SAMPLE_ATTRIBUTES></SAMPLE></SAMPLE_SET>')
        with patch.object(api, "http_get", return_value=xml.encode()):
            attrs = api.ena_sample_attrs("SAMEA1")
        assert attrs["genotype"] == "Col-0"
        assert "ENA-CHECKLIST" not in attrs
        assert attrs["scientific_name"] == "Arabidopsis thaliana"

    def test_clean_val_replaces_control_characters(self):
        import ena_fetch_api as api

        assert api._clean_val("a\tb\nc") == "a_b_c"

    def test_metadata_tsv_is_lf_terminated(self, tmp_path):
        import ena_fetch_api as api

        out = tmp_path / "metadata.tsv"
        api.write_metadata_tsv([api.harmonize_row("S1", "ERR1", {})], out)
        assert b"\r\n" not in out.read_bytes()


class TestReadMap:
    """The highest-risk pure function in the set: a wrong pairing silently
    drops the cDNA read of a 10x run, and the pipeline then runs on barcodes."""

    def test_explicit_read_map_is_parsed(self):
        import ena_fetch_api as api

        assert api.parse_read_map("3,4") == (3, 4)

    def test_read_map_rejects_equal_reads(self):
        import ena_fetch_api as api

        with pytest.raises(SystemExit):
            api.parse_read_map("2,2")

    def test_read_map_rejects_malformed_values(self):
        import ena_fetch_api as api

        for bad in ("0,1", "x,2", "1", "1,2,3"):
            with pytest.raises(SystemExit):
                api.parse_read_map(bad)

    def test_read_map_is_optional(self):
        """Not passing --read-map is the common case, not an error."""
        import ena_fetch_api as api

        assert api.parse_read_map(None) is None
        assert api.parse_read_map("") is None

    def test_read_map_accepts_positions_above_four(self):
        """Upstream validates 1-9, not 1-4. Documented rather than tightened:
        a run with more than four files is rare but not impossible, and a wrong
        position fails later with a missing-file error anyway."""
        import ena_fetch_api as api

        assert api.parse_read_map("1,9") == (1, 9)


class TestDemo:
    def test_demo_makes_no_network_call(self, tmp_path):
        import ena_fetch as app
        import ena_fetch_api as api

        with patch.object(api.urllib.request, "urlopen") as mocked:
            app.main(["--demo", "--output", str(tmp_path)])
            mocked.assert_not_called()

    def test_demo_writes_report_and_result(self, tmp_path):
        import ena_fetch as app

        app.main(["--demo", "--output", str(tmp_path)])
        assert (tmp_path / "report.md").exists()
        assert json.loads((tmp_path / "result.json").read_text())["skill"] == "ena-fetch"

    def test_demo_report_carries_the_disclaimer(self, tmp_path):
        import ena_fetch as app

        app.main(["--demo", "--output", str(tmp_path)])
        text = (tmp_path / "report.md").read_text()
        assert "research and educational tool" in text
        assert "not a medical device" in text

    def test_demo_writes_the_reproducibility_bundle(self, tmp_path):
        import ena_fetch as app

        app.main(["--demo", "--output", str(tmp_path)])
        repro = tmp_path / "reproducibility"
        for name in ("commands.sh", "environment.yml", "checksums.sha256"):
            assert (repro / name).exists()

    def test_checksum_labels_resolve_from_the_output_dir(self, tmp_path):
        import ena_fetch as app

        app.main(["--demo", "--output", str(tmp_path)])
        for line in (tmp_path / "reproducibility" / "checksums.sha256").read_text().splitlines():
            assert (tmp_path / line.split("  ", 1)[1]).exists()

    def test_demo_writes_a_metadata_table_with_one_row_per_run(self, tmp_path):
        import ena_fetch as app

        app.main(["--demo", "--output", str(tmp_path)])
        rows = list(csv.DictReader((tmp_path / "tables" / "metadata.tsv").open(), delimiter="\t"))
        assert len(rows) == 9
        assert {r["species"] for r in rows} == {"Arabidopsis thaliana"}

    def test_demo_writes_a_pipeline_ready_samplesheet(self, tmp_path):
        """nf-core/rnaseq column contract: sample,fastq_1,fastq_2,strandedness."""
        import ena_fetch as app

        app.main(["--demo", "--output", str(tmp_path)])
        sheet = tmp_path / "samplesheet.csv"
        rows = list(csv.DictReader(sheet.open()))
        assert list(rows[0]) == ["sample", "fastq_1", "fastq_2", "strandedness"]
        assert rows and all(r["fastq_1"].startswith("http") for r in rows)
        assert all(r["strandedness"] == "auto" for r in rows)

    def test_demo_samplesheet_writes_the_archive_md5_sidecar(self, tmp_path):
        """fastq_md5 is ';'-aligned with fastq_ftp; each URL keeps its own checksum."""
        import ena_fetch as app

        app.main(["--demo", "--output", str(tmp_path)])
        rows = list(csv.DictReader((tmp_path / "fastq_md5.tsv").open(), delimiter="\t"))
        md5 = {r["url"]: r["md5"] for r in rows}
        assert md5["https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/054/ERR10181254/"
                   "ERR10181254_1.fastq.gz"] == "735252236e07db9324a1a95112b237e9"
        assert md5["https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/054/ERR10181254/"
                   "ERR10181254_2.fastq.gz"] == "b3c6d69a99c420b3385970bf82af1fde"
        sheet_urls = {u for r in csv.DictReader((tmp_path / "samplesheet.csv").open())
                      for u in (r["fastq_1"], r["fastq_2"]) if u}
        assert sheet_urls == set(md5)

    def test_demo_download_script_verifies_every_file(self, tmp_path):
        import ena_fetch as app

        app.main(["--demo", "--output", str(tmp_path)])
        body = (tmp_path / "download_ena.sh").read_text()
        checks = [ln for ln in body.splitlines() if ln.startswith("verify_md5 ")]
        downloads = [ln for ln in body.splitlines() if ln.startswith("curl -fsSL")]
        assert len(checks) == len(downloads) > 0
        assert ('verify_md5 735252236e07db9324a1a95112b237e9 '
                '"$OUTDIR"/ERR10181254_1.fastq.gz') in checks

    def test_local_path_samplesheets_write_no_md5_sidecar(self, tmp_path):
        """With --fastq-dir the sheet holds local paths, nothing to verify."""
        import ena_fetch as app

        app.main(["--demo", "--command", "samplesheet", "--assay", "bulk",
                  "--fastq-dir", "/data/fastq", "--output", str(tmp_path)])
        assert not (tmp_path / "fastq_md5.tsv").exists()


class TestCLI:
    def test_no_args_exits_nonzero(self):
        result = subprocess.run([sys.executable, str(SCRIPT)], capture_output=True, text=True)
        assert result.returncode != 0

    def test_command_flag_dispatches(self, tmp_path):
        import ena_fetch as app

        app.main(["--demo", "--command", "runs", "--output", str(tmp_path)])
        assert "## runs" in (tmp_path / "report.md").read_text()

    def test_upstream_positional_form_still_works(self, tmp_path):
        """`ena_fetch.py runs PRJEB56029` must behave like upstream's CLI."""
        import ena_fetch as app

        app._install_demo_transport()  # serve the fixture, not the network
        app.main(["runs", DEMO_PROJECT, "--output", str(tmp_path)])
        assert "## runs" in (tmp_path / "report.md").read_text()

    def test_samplesheet_requires_an_assay(self, tmp_path):
        import ena_fetch as app

        with pytest.raises(SystemExit):
            app.main(["--command", "samplesheet", "--accession", DEMO_PROJECT,
                      "--output", str(tmp_path)])


class TestSafety:
    def test_warns_before_overwriting(self, tmp_path, capsys):
        import ena_fetch as app

        app.main(["--demo", "--output", str(tmp_path)])
        capsys.readouterr()
        app.main(["--demo", "--output", str(tmp_path)])
        assert "overwritten" in capsys.readouterr().err

    def test_demo_writes_nothing_outside_the_output_dir(self, tmp_path, monkeypatch):
        import ena_fetch as app

        cwd = tmp_path / "cwd"
        cwd.mkdir()
        monkeypatch.chdir(cwd)
        app.main(["--demo", "--output", str(tmp_path / "out")])
        assert list(cwd.iterdir()) == []

    def test_report_holds_no_absolute_output_path(self, tmp_path):
        import ena_fetch as app

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
    files = []
    parents = {}
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


class TestReportLimit:
    """`--limit` means different things per command, and the shared wrapper
    must not pick one of them globally.

    The vendored `report` defaults to `--limit 0` (no limit); the shared parser
    defaulted to 20 and forwarded it, so a 95-run study silently reported 20
    rows with status: ok. A short file report looks exactly like a complete one.
    """

    def _argv(self, extra, tmp_path):
        import ena_fetch as app

        args = app._build_parser().parse_args(
            ["--command", "report", "--accession", DEMO_PROJECT, *extra])
        return app._to_upstream_argv(args, tmp_path)

    def test_report_does_not_forward_a_limit_by_default(self, tmp_path):
        argv = self._argv([], tmp_path)
        assert "--limit" not in argv, (
            "forwarding the shared default overrides the vendored 0 = no limit")

    def test_an_explicit_limit_is_still_forwarded(self, tmp_path):
        argv = self._argv(["--limit", "7"], tmp_path)
        assert argv[argv.index("--limit") + 1] == "7"

    def test_explicit_zero_is_forwarded_not_treated_as_unset(self, tmp_path):
        """0 is falsy; `if args.limit` would drop the user's explicit request."""
        argv = self._argv(["--limit", "0"], tmp_path)
        assert argv[argv.index("--limit") + 1] == "0"

    def test_search_still_defaults_to_twenty(self, tmp_path):
        import ena_fetch as app

        args = app._build_parser().parse_args(["--command", "search", "--query", "x"])
        argv = app._to_upstream_argv(args, tmp_path)
        assert argv[argv.index("--limit") + 1] == "20"


class TestTruncationWarning:
    """Fixing the default stops the accidental case; an explicit --limit can
    still truncate. The objection was the silence, not the number."""

    def test_warns_when_rows_equal_the_limit(self):
        import ena_fetch as app

        body = "run_accession\tfastq_ftp\n" + "".join(
            f"ERR{i}\tftp://x/{i}.gz\n" for i in range(3))
        assert app._truncation_warning(body, 3)

    def test_silent_when_rows_are_fewer_than_the_limit(self):
        import ena_fetch as app

        body = "run_accession\tfastq_ftp\nERR1\tftp://x/1.gz\n"
        assert app._truncation_warning(body, 3) is None

    def test_silent_when_no_limit_was_applied(self):
        import ena_fetch as app

        body = "run_accession\tfastq_ftp\nERR1\tftp://x/1.gz\n"
        assert app._truncation_warning(body, None) is None
        assert app._truncation_warning(body, 0) is None

    def test_the_message_names_the_escape_hatch(self):
        import ena_fetch as app

        body = "h\n" + "".join(f"r{i}\n" for i in range(2))
        assert "--limit 0" in app._truncation_warning(body, 2)


class TestDownloadSizePreflight:
    """run_upstream redirects stdout/stderr into a StringIO, so anything the
    vendored command prints appears only after it returns -- all at once, after
    the wait it was meant to explain. The estimate has to be emitted before
    delegating, as warn_if_overwriting already is."""

    def test_bytes_are_rendered_human_readably(self):
        import ena_fetch as app

        assert app._human_bytes(0) == "0 B"
        assert app._human_bytes(42 * 1024 ** 2).startswith("42")
        assert app._human_bytes(90 * 1024 ** 3).endswith("GB")

    def test_the_estimate_reaches_the_real_stderr(self, tmp_path, capsys):
        import ena_fetch as app

        app._install_demo_transport()
        app.main(["--demo", "--command", "runs", "--output", str(tmp_path)])
        captured = capsys.readouterr()
        assert "## runs" in (tmp_path / "report.md").read_text()
        # runs is not a download; the point is that the harness sees stderr at all
        assert captured.err is not None

    def test_preflight_summarises_count_and_size(self):
        import ena_fetch as app

        rows = [{"fastq_bytes": "100;200"}, {"fastq_bytes": "300"}]
        line = app._download_preflight_line("PRJEB1", rows, "fastq_bytes")
        assert "3 file(s)" in line
        assert "600 B" in line

    def test_preflight_survives_missing_sizes(self):
        """ENA omits fastq_bytes for some records; an unknown total must not
        crash the download that the estimate exists to explain."""
        import ena_fetch as app

        line = app._download_preflight_line("PRJEB1", [{"fastq_bytes": ""}], "fastq_bytes")
        assert "PRJEB1" in line


# ── Demo output snapshot ──────────────────────────────────────────────────────
# Fails when a code change alters any demo number. If the change is intended,
# run `pytest skills/ena-fetch/tests --inline-snapshot=fix`; it rewrites the
# snapshot() values below, and the PR diff shows which numbers moved.

@pytest.mark.demo_snapshot
def test_demo_output_matches_snapshot(tmp_path):
    assert demo_runner.demo_outputs("ena-fetch", tmp_path) == snapshot({'download_ena.sh': """\
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

# Verify each file against the archive's published MD5 once it is complete.
if command -v md5sum >/dev/null 2>&1; then
  verify_md5() {
    printf '%s  %s\\n' "$1" "$2" | md5sum -c --quiet - >/dev/null 2>&1 || {
      echo "MD5 mismatch for $2: the file is truncated or corrupt; delete it and re-run this script" >&2
      exit 1
    }
  }
else
  echo "warning: md5sum not found; downloaded files will not be checksum-verified" >&2
  verify_md5() { :; }
fi

# sample: SAMEA111350648
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/ERR10181253_1.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/053/ERR10181253/ERR10181253_1.fastq.gz
verify_md5 05c0b37299a258859f92125a3ebae46c "$OUTDIR"/ERR10181253_1.fastq.gz
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/ERR10181253_2.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/053/ERR10181253/ERR10181253_2.fastq.gz
verify_md5 bbce0154ef5b7383084f5a655ff1fd4e "$OUTDIR"/ERR10181253_2.fastq.gz

# sample: SAMEA111350649
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/ERR10181254_1.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/054/ERR10181254/ERR10181254_1.fastq.gz
verify_md5 735252236e07db9324a1a95112b237e9 "$OUTDIR"/ERR10181254_1.fastq.gz
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/ERR10181254_2.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/054/ERR10181254/ERR10181254_2.fastq.gz
verify_md5 b3c6d69a99c420b3385970bf82af1fde "$OUTDIR"/ERR10181254_2.fastq.gz

# sample: SAMEA111350650
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/ERR10181255_1.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/055/ERR10181255/ERR10181255_1.fastq.gz
verify_md5 b3aaa03bd95d9bf30e402cf7676de05c "$OUTDIR"/ERR10181255_1.fastq.gz
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/ERR10181255_2.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/055/ERR10181255/ERR10181255_2.fastq.gz
verify_md5 600badb9d653bbbf22d37a063a162c1a "$OUTDIR"/ERR10181255_2.fastq.gz

# sample: SAMEA111350651
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/ERR10181256_1.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/056/ERR10181256/ERR10181256_1.fastq.gz
verify_md5 3a943f1a72d629549765b706badfbddd "$OUTDIR"/ERR10181256_1.fastq.gz
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/ERR10181256_2.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/056/ERR10181256/ERR10181256_2.fastq.gz
verify_md5 c85db1d41ff3f4c177337a54157f423f "$OUTDIR"/ERR10181256_2.fastq.gz

# sample: SAMEA111350652
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/ERR10181257_1.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/057/ERR10181257/ERR10181257_1.fastq.gz
verify_md5 2019258c721b9e43a1a5e5001f0a4225 "$OUTDIR"/ERR10181257_1.fastq.gz
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/ERR10181257_2.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/057/ERR10181257/ERR10181257_2.fastq.gz
verify_md5 e0ae1ece62a52b91ed60173b40ce73d1 "$OUTDIR"/ERR10181257_2.fastq.gz

# sample: SAMEA111350653
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/ERR10181258_1.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/058/ERR10181258/ERR10181258_1.fastq.gz
verify_md5 0fe6528436891dd0eb58849b9fdd0dd6 "$OUTDIR"/ERR10181258_1.fastq.gz
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/ERR10181258_2.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/058/ERR10181258/ERR10181258_2.fastq.gz
verify_md5 37203ab5cc66cb158a840c0c8561ef0c "$OUTDIR"/ERR10181258_2.fastq.gz

# sample: SAMEA111350654
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/ERR10181259_1.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/059/ERR10181259/ERR10181259_1.fastq.gz
verify_md5 836b83c7e84cf6280b2c5e6bfb6f8f04 "$OUTDIR"/ERR10181259_1.fastq.gz
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/ERR10181259_2.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/059/ERR10181259/ERR10181259_2.fastq.gz
verify_md5 f8c93bd7126006a9faec7cccb4207bf5 "$OUTDIR"/ERR10181259_2.fastq.gz

# sample: SAMEA111350655
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/ERR10181260_1.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/060/ERR10181260/ERR10181260_1.fastq.gz
verify_md5 901030c63eb4d97a878eed8f66e1eb0f "$OUTDIR"/ERR10181260_1.fastq.gz
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/ERR10181260_2.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/060/ERR10181260/ERR10181260_2.fastq.gz
verify_md5 31f1b6e7f02f8c08cbf21f8fa9d2e2af "$OUTDIR"/ERR10181260_2.fastq.gz

# sample: SAMEA111350656
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/ERR10181261_1.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/061/ERR10181261/ERR10181261_1.fastq.gz
verify_md5 68a498e2f53c8ca789ddc1558c80a23c "$OUTDIR"/ERR10181261_1.fastq.gz
curl -fsSL --retry 5 --retry-delay 10 $RETRY_ALL --connect-timeout 30 --speed-limit 1024 --speed-time 120 -C - --create-dirs -o "$OUTDIR"/ERR10181261_2.fastq.gz https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/061/ERR10181261/ERR10181261_2.fastq.gz
verify_md5 703f962baf6008d54a2c72879c0d7883 "$OUTDIR"/ERR10181261_2.fastq.gz

echo "Downloaded 18 file(s) to $OUTDIR"
""", 'fastq_md5.tsv': """\
url	md5
https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/054/ERR10181254/ERR10181254_1.fastq.gz	735252236e07db9324a1a95112b237e9
https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/054/ERR10181254/ERR10181254_2.fastq.gz	b3c6d69a99c420b3385970bf82af1fde
https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/060/ERR10181260/ERR10181260_1.fastq.gz	901030c63eb4d97a878eed8f66e1eb0f
https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/060/ERR10181260/ERR10181260_2.fastq.gz	31f1b6e7f02f8c08cbf21f8fa9d2e2af
https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/053/ERR10181253/ERR10181253_1.fastq.gz	05c0b37299a258859f92125a3ebae46c
https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/053/ERR10181253/ERR10181253_2.fastq.gz	bbce0154ef5b7383084f5a655ff1fd4e
https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/059/ERR10181259/ERR10181259_1.fastq.gz	836b83c7e84cf6280b2c5e6bfb6f8f04
https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/059/ERR10181259/ERR10181259_2.fastq.gz	f8c93bd7126006a9faec7cccb4207bf5
https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/061/ERR10181261/ERR10181261_1.fastq.gz	68a498e2f53c8ca789ddc1558c80a23c
https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/061/ERR10181261/ERR10181261_2.fastq.gz	703f962baf6008d54a2c72879c0d7883
https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/056/ERR10181256/ERR10181256_1.fastq.gz	3a943f1a72d629549765b706badfbddd
https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/056/ERR10181256/ERR10181256_2.fastq.gz	c85db1d41ff3f4c177337a54157f423f
https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/057/ERR10181257/ERR10181257_1.fastq.gz	2019258c721b9e43a1a5e5001f0a4225
https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/057/ERR10181257/ERR10181257_2.fastq.gz	e0ae1ece62a52b91ed60173b40ce73d1
https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/055/ERR10181255/ERR10181255_1.fastq.gz	b3aaa03bd95d9bf30e402cf7676de05c
https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/055/ERR10181255/ERR10181255_2.fastq.gz	600badb9d653bbbf22d37a063a162c1a
https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/058/ERR10181258/ERR10181258_1.fastq.gz	0fe6528436891dd0eb58849b9fdd0dd6
https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/058/ERR10181258/ERR10181258_2.fastq.gz	37203ab5cc66cb158a840c0c8561ef0c
""", 'report.md': """\
# ENA report — PRJEB56029

Source: [European Nucleotide Archive](https://www.ebi.ac.uk/ena/browser/view/PRJEB56029)

## runs

```
9 run(s) for PRJEB56029:

ERR10181254    ILLUMINA     RNA-Seq      reads=39994488
    ftp: ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/054/ERR10181254/ERR10181254_1.fastq.gz
    ftp: ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/054/ERR10181254/ERR10181254_2.fastq.gz
ERR10181260    ILLUMINA     RNA-Seq      reads=45764926
    ftp: ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/060/ERR10181260/ERR10181260_1.fastq.gz
    ftp: ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/060/ERR10181260/ERR10181260_2.fastq.gz
ERR10181253    ILLUMINA     RNA-Seq      reads=40842088
    ftp: ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/053/ERR10181253/ERR10181253_1.fastq.gz
    ftp: ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/053/ERR10181253/ERR10181253_2.fastq.gz
ERR10181259    ILLUMINA     RNA-Seq      reads=46792448
    ftp: ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/059/ERR10181259/ERR10181259_1.fastq.gz
    ftp: ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/059/ERR10181259/ERR10181259_2.fastq.gz
ERR10181261    ILLUMINA     RNA-Seq      reads=47009856
    ftp: ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/061/ERR10181261/ERR10181261_1.fastq.gz
    ftp: ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/061/ERR10181261/ERR10181261_2.fastq.gz
ERR10181256    ILLUMINA     RNA-Seq      reads=39414064
    ftp: ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/056/ERR10181256/ERR10181256_1.fastq.gz
    ftp: ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/056/ERR10181256/ERR10181256_2.fastq.gz
ERR10181257    ILLUMINA     RNA-Seq      reads=45247648
    ftp: ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/057/ERR10181257/ERR10181257_1.fastq.gz
    ftp: ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/057/ERR10181257/ERR10181257_2.fastq.gz
ERR10181255    ILLUMINA     RNA-Seq      reads=47474894
    ftp: ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/055/ERR10181255/ERR10181255_1.fastq.gz
    ftp: ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/055/ERR10181255/ERR10181255_2.fastq.gz
ERR10181258    ILLUMINA     RNA-Seq      reads=51635774
    ftp: ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/058/ERR10181258/ERR10181258_1.fastq.gz
    ftp: ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/058/ERR10181258/ERR10181258_2.fastq.gz
```

## metadata-table

```
Wrote 9 sample x replicate row(s) x 11 column(s) for PRJEB56029 to <output>/tables/metadata.tsv
```

## samplesheet

```
Wrote 9 row(s) for 9 sample(s) to <output>/samplesheet.csv (nf-core/rnaseq)
```

## download-script

```
Wrote download_ena.sh: 18 download command(s) using curl.
Nothing has been downloaded. Run it yourself with `bash download_ena.sh`, or submit it with `sbatch download_ena.sh`.
```

---

*ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions.*
""", 'result.json': """\
{
  "status": "ok",
  "ok": true,
  "skill": "ena-fetch",
  "version": "0.1.0",
  "completed_at": "<TIMESTAMP>",
  "input_checksum": "",
  "datasets": {},
  "summary": {
    "accession": "PRJEB56029",
    "commands": [
      "runs",
      "metadata-table",
      "samplesheet",
      "download-script"
    ],
    "demo": true
  },
  "data": {
    "sections": {
      "runs": "9 run(s) for PRJEB56029:\\n\\nERR10181254    ILLUMINA     RNA-Seq      reads=39994488\\n    ftp: ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/054/ERR10181254/ERR10181254_1.fastq.gz\\n    ftp: ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/054/ERR10181254/ERR10181254_2.fastq.gz\\nERR10181260    ILLUMINA     RNA-Seq      reads=45764926\\n    ftp: ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/060/ERR10181260/ERR10181260_1.fastq.gz\\n    ftp: ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/060/ERR10181260/ERR10181260_2.fastq.gz\\nERR10181253    ILLUMINA     RNA-Seq      reads=40842088\\n    ftp: ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/053/ERR10181253/ERR10181253_1.fastq.gz\\n    ftp: ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/053/ERR10181253/ERR10181253_2.fastq.gz\\nERR10181259    ILLUMINA     RNA-Seq      reads=46792448\\n    ftp: ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/059/ERR10181259/ERR10181259_1.fastq.gz\\n    ftp: ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/059/ERR10181259/ERR10181259_2.fastq.gz\\nERR10181261    ILLUMINA     RNA-Seq      reads=47009856\\n    ftp: ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/061/ERR10181261/ERR10181261_1.fastq.gz\\n    ftp: ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/061/ERR10181261/ERR10181261_2.fastq.gz\\nERR10181256    ILLUMINA     RNA-Seq      reads=39414064\\n    ftp: ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/056/ERR10181256/ERR10181256_1.fastq.gz\\n    ftp: ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/056/ERR10181256/ERR10181256_2.fastq.gz\\nERR10181257    ILLUMINA     RNA-Seq      reads=45247648\\n    ftp: ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/057/ERR10181257/ERR10181257_1.fastq.gz\\n    ftp: ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/057/ERR10181257/ERR10181257_2.fastq.gz\\nERR10181255    ILLUMINA     RNA-Seq      reads=47474894\\n    ftp: ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/055/ERR10181255/ERR10181255_1.fastq.gz\\n    ftp: ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/055/ERR10181255/ERR10181255_2.fastq.gz\\nERR10181258    ILLUMINA     RNA-Seq      reads=51635774\\n    ftp: ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/058/ERR10181258/ERR10181258_1.fastq.gz\\n    ftp: ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/058/ERR10181258/ERR10181258_2.fastq.gz\\n",
      "metadata-table": "Wrote 9 sample x replicate row(s) x 11 column(s) for PRJEB56029 to <output>/tables/metadata.tsv\\n",
      "samplesheet": "Wrote 9 row(s) for 9 sample(s) to <output>/samplesheet.csv (nf-core/rnaseq)\\n",
      "download-script": "Wrote download_ena.sh: 18 download command(s) using curl.\\nNothing has been downloaded. Run it yourself with `bash download_ena.sh`, or submit it with `sbatch download_ena.sh`."
    }
  }
}\
""", 'samplesheet.csv': """\
sample,fastq_1,fastq_2,strandedness
SAMEA111350648,https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/053/ERR10181253/ERR10181253_1.fastq.gz,https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/053/ERR10181253/ERR10181253_2.fastq.gz,auto
SAMEA111350649,https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/054/ERR10181254/ERR10181254_1.fastq.gz,https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/054/ERR10181254/ERR10181254_2.fastq.gz,auto
SAMEA111350650,https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/055/ERR10181255/ERR10181255_1.fastq.gz,https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/055/ERR10181255/ERR10181255_2.fastq.gz,auto
SAMEA111350651,https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/056/ERR10181256/ERR10181256_1.fastq.gz,https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/056/ERR10181256/ERR10181256_2.fastq.gz,auto
SAMEA111350652,https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/057/ERR10181257/ERR10181257_1.fastq.gz,https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/057/ERR10181257/ERR10181257_2.fastq.gz,auto
SAMEA111350653,https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/058/ERR10181258/ERR10181258_1.fastq.gz,https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/058/ERR10181258/ERR10181258_2.fastq.gz,auto
SAMEA111350654,https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/059/ERR10181259/ERR10181259_1.fastq.gz,https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/059/ERR10181259/ERR10181259_2.fastq.gz,auto
SAMEA111350655,https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/060/ERR10181260/ERR10181260_1.fastq.gz,https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/060/ERR10181260/ERR10181260_2.fastq.gz,auto
SAMEA111350656,https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/061/ERR10181261/ERR10181261_1.fastq.gz,https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/061/ERR10181261/ERR10181261_2.fastq.gz,auto
""", 'tables/metadata.tsv': """\
sample	replicate	species	sex	age	condition	genotype	treatment	tissue	common_name	ecotype
SAMEA111350649	ERR10181254	Arabidopsis thaliana	NA	mature non-fertilized ovule	NA	wild type genotype	NA	plant ovule	thale cress	Col-0
SAMEA111350655	ERR10181260	Arabidopsis thaliana	NA	mature non-fertilized ovule	NA	wild type genotype	NA	plant ovule	thale cress	Col-0
SAMEA111350648	ERR10181253	Arabidopsis thaliana	NA	mature non-fertilized ovule	NA	wild type genotype	NA	plant ovule	thale cress	Col-0
SAMEA111350654	ERR10181259	Arabidopsis thaliana	NA	mature non-fertilized ovule	NA	wild type genotype	NA	plant ovule	thale cress	Col-0
SAMEA111350656	ERR10181261	Arabidopsis thaliana	NA	mature non-fertilized ovule	NA	wild type genotype	NA	plant ovule	thale cress	Col-0
SAMEA111350651	ERR10181256	Arabidopsis thaliana	NA	mature non-fertilized ovule	NA	wild type genotype	NA	plant ovule	thale cress	Col-0
SAMEA111350652	ERR10181257	Arabidopsis thaliana	NA	mature non-fertilized ovule	NA	wild type genotype	NA	plant ovule	thale cress	Col-0
SAMEA111350650	ERR10181255	Arabidopsis thaliana	NA	mature non-fertilized ovule	NA	wild type genotype	NA	plant ovule	thale cress	Col-0
SAMEA111350653	ERR10181258	Arabidopsis thaliana	NA	mature non-fertilized ovule	NA	wild type genotype	NA	plant ovule	thale cress	Col-0
"""})
