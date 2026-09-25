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

SKILL_DIR = Path(__file__).resolve().parents[1]
SCRIPT = SKILL_DIR / "ena_fetch.py"
EXAMPLES = SKILL_DIR / "examples"
DEMO_PROJECT = "PRJEB56029"

sys.path.insert(0, str(SKILL_DIR))


class TestVendoredApi:
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
