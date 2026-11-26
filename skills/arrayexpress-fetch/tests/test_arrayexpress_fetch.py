"""Tests for arrayexpress-fetch.

Run with: pytest skills/arrayexpress-fetch/tests/test_arrayexpress_fetch.py -v

No network required: every HTTP call is mocked against the committed demo
fixtures in ../examples/.
"""

import csv
import json
import re
import subprocess
import sys
import urllib.parse
from pathlib import Path
from unittest.mock import patch

import pytest
from inline_snapshot import snapshot
from clawbio.common import demo_runner

SKILL_DIR = Path(__file__).resolve().parents[1]
SCRIPT = SKILL_DIR / "arrayexpress_fetch.py"
EXAMPLES = SKILL_DIR / "examples"
DEMO_ACCESSION = "E-MTAB-10030"
DEMO_STUDY = EXAMPLES / f"demo_{DEMO_ACCESSION}.json"
DEMO_SDRF = EXAMPLES / f"demo_{DEMO_ACCESSION}.sdrf.txt"

sys.path.insert(0, str(SKILL_DIR))


def _load(name):
    return json.loads((EXAMPLES / name).read_text())


class _FakeResponse:
    """Minimal stand-in for the object urlopen returns as a context manager."""

    def __init__(self, payload: bytes):
        self._payload = payload

    def read(self, n=None):
        if n is None:
            payload, self._payload = self._payload, b""
            return payload
        chunk, self._payload = self._payload[:n], self._payload[n:]
        return chunk

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


def _fake_urlopen(payload: dict | bytes):
    """Return a urlopen replacement that records the URLs it was asked for."""
    body = payload if isinstance(payload, bytes) else json.dumps(payload).encode()
    calls = []

    def _open(req, *args, **kwargs):
        calls.append(getattr(req, "full_url", req))
        return _FakeResponse(body)

    _open.calls = calls
    return _open


@pytest.fixture(autouse=True)
def _clear_info_cache():
    """file_url() memoises /info per accession; keep that out of other tests."""
    import arrayexpress_fetch_api as api

    getattr(api, "_INFO_CACHE", {}).clear()
    yield
    getattr(api, "_INFO_CACHE", {}).clear()


# --------------------------------------------------------------------------
# Vendored API module: upstream logic, behaviour preserved
# --------------------------------------------------------------------------


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
        import arrayexpress_fetch_api as api

        assert api.to_https(value) == expected

    def test_walk_files_finds_every_file_node(self):
        import arrayexpress_fetch_api as api

        files = list(api.walk_files(_load(f"demo_{DEMO_ACCESSION}.json")["section"]))
        assert files, "the fixture must expose at least one file node"
        assert all("path" in f for f in files)

    def test_classify_labels_magetab_and_raw(self):
        import arrayexpress_fetch_api as api

        assert api.classify({"path": "E-MTAB-10030.sdrf.txt"}) == "sdrf"
        assert api.classify({"path": "E-MTAB-10030.idf.txt"}) == "idf"
        assert api.classify({"path": "B1_S1_R1.fastq.gz"}) == "raw"

    def test_clean_sample_replaces_control_characters(self):
        import arrayexpress_fetch_api as api

        # "Clean output fields": control chars, including tab/CR/LF, become _
        assert "\t" not in api.clean_sample("Sample\t2")
        assert "\n" not in api.clean_sample("Sample\n2")

    def test_clean_sample_makes_ids_pipeline_safe(self):
        import arrayexpress_fetch_api as api

        assert api.clean_sample("Sample 2") == "Sample_2"

    def test_parse_read_map_accepts_a_pair(self):
        import arrayexpress_fetch_api as api

        assert api.parse_read_map("1,2") == (1, 2)

    def test_parse_read_map_rejects_nonsense(self):
        import arrayexpress_fetch_api as api

        with pytest.raises(SystemExit):
            api.parse_read_map("banana")

    def test_to_https_rewrites_ftp_urls(self):
        import arrayexpress_fetch_api as api

        # The skill never speaks FTP; every ftp:// URL is rewritten to https://
        assert api.to_https("ftp://ftp.ebi.ac.uk/x.gz") == "https://ftp.ebi.ac.uk/x.gz"

    def test_harmonize_row_maps_organism_onto_species(self):
        import arrayexpress_fetch_api as api

        row = api.harmonize_row("S1", "1", {"organism": "Rattus norvegicus"})
        assert row["species"] == "Rattus norvegicus"

    def test_harmonize_row_uses_na_for_absent_fields(self):
        import arrayexpress_fetch_api as api

        row = api.harmonize_row("S1", "1", {})
        assert row["sex"] == api.NA


class TestDownloadBase:
    """The download base is resolved from /studies/{acc}/info, not hardcoded.

    Same divergence as biostudies-fetch, for the same reason: `httpLink` names
    a different tree per collection (fire/ for E-MTAB, pub/databases/ for
    S-BSST), so no single constant is correct. Verified live 2026-09-22.
    """

    def test_file_url_is_resolved_from_the_info_endpoint(self):
        import arrayexpress_fetch_api as api

        info = {"httpLink": "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/030/E-MTAB-10030"}
        with patch.object(api, "get_json", return_value=info):
            url = api.file_url(DEMO_ACCESSION, "E-MTAB-10030.sdrf.txt")
        assert url == (
            "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/030/E-MTAB-10030"
            "/Files/E-MTAB-10030.sdrf.txt"
        )

    def test_file_url_percent_encodes_under_the_resolved_base(self):
        import arrayexpress_fetch_api as api

        info = {"httpLink": "https://ftp.ebi.ac.uk/biostudies/fire/E-M/001/E-MTAB-1"}
        with patch.object(api, "get_json", return_value=info):
            assert api.file_url("E-MTAB-1", "a dir/b.txt").endswith("/Files/a%20dir/b.txt")

    def test_file_url_falls_back_when_info_omits_httplink(self):
        import arrayexpress_fetch_api as api

        with patch.object(api, "get_json", return_value={}):
            assert api.file_url("E-MTAB-1", "b.txt") == f"{api.FILES}/E-MTAB-1/b.txt"

    @pytest.mark.parametrize("link", [
        "https://evil.example/x", "https://ebi.ac.uk.evil.com/x",
        "http://ftp.ebi.ac.uk/x", "file:///etc",
    ])
    def test_untrusted_httplink_falls_back_to_the_legacy_path(self, link):
        """The info payload is server data: a link off ebi.ac.uk, or not https,
        must not decide where bytes are fetched from."""
        import arrayexpress_fetch_api as api

        with patch.object(api, "get_json", return_value={"httpLink": link}):
            assert api.file_url("E-MTAB-1", "b.txt") == f"{api.FILES}/E-MTAB-1/b.txt"


# --------------------------------------------------------------------------
# CLI contract
# --------------------------------------------------------------------------


class TestCLI:
    def test_demo_requires_no_subcommand(self, tmp_path):
        import arrayexpress_fetch as app

        assert app.main(["--demo", "--output", str(tmp_path)]) == 0

    def test_command_form_is_runner_reachable(self, tmp_path):
        """clawbio/cli.py drops bare positionals, so --command must exist."""
        import arrayexpress_fetch as app

        parser = app._build_parser()
        args = parser.parse_args(["--command", "metadata", "--accession", DEMO_ACCESSION])
        assert args.command == "metadata"

    def test_upstream_positional_form_still_works(self):
        import arrayexpress_fetch as app
        from clawbio.common import archive_fetch as af

        argv = af.expand_positional(["metadata", DEMO_ACCESSION], app.COMMANDS)
        assert argv[:2] == ["--command", "metadata"]
        assert DEMO_ACCESSION in argv

    def test_samplesheet_requires_an_assay(self):
        import arrayexpress_fetch as app

        parser = app._build_parser()
        args = parser.parse_args(["--command", "samplesheet", "--accession", DEMO_ACCESSION])
        with pytest.raises(SystemExit):
            app._to_upstream_argv(args, Path("/tmp"))

    def test_search_without_query_is_rejected(self):
        import arrayexpress_fetch as app

        parser = app._build_parser()
        args = parser.parse_args(["--command", "search"])
        with pytest.raises(SystemExit):
            app._to_upstream_argv(args, Path("/tmp"))

    def test_metadata_without_accession_is_rejected(self):
        import arrayexpress_fetch as app

        parser = app._build_parser()
        args = parser.parse_args(["--command", "metadata"])
        with pytest.raises(SystemExit):
            app._to_upstream_argv(args, Path("/tmp"))

    def test_missing_output_is_rejected(self):
        import arrayexpress_fetch as app

        with pytest.raises(SystemExit):
            app.main(["--command", "metadata", "--accession", DEMO_ACCESSION])

    def test_relative_out_is_anchored_under_output(self, tmp_path):
        """Upstream defaulted to cwd; nothing may escape --output."""
        from clawbio.common import archive_fetch as af

        resolved = af.resolve_out("sheet.csv", tmp_path, "samplesheet.csv")
        assert resolved == tmp_path / "sheet.csv"


# --------------------------------------------------------------------------
# Demo mode
# --------------------------------------------------------------------------


class TestDemo:
    def test_demo_writes_the_documented_tree(self, tmp_path):
        import arrayexpress_fetch as app

        app.main(["--demo", "--output", str(tmp_path)])
        for rel in ("report.md", "result.json", "tables/metadata.tsv",
                    "reproducibility/commands.sh",
                    "reproducibility/environment.yml",
                    "reproducibility/checksums.sha256"):
            assert (tmp_path / rel).exists(), f"missing {rel}"

    def test_demo_makes_no_network_call(self, tmp_path):
        import arrayexpress_fetch as app

        with patch("urllib.request.urlopen",
                   side_effect=AssertionError("demo mode must not hit the network")):
            assert app.main(["--demo", "--output", str(tmp_path)]) == 0

    def test_demo_samplesheet_is_written(self, tmp_path):
        import arrayexpress_fetch as app

        app.main(["--demo", "--output", str(tmp_path)])
        sheet = tmp_path / "samplesheet.csv"
        assert sheet.exists(), "the demo must exercise the samplesheet path"
        header = sheet.read_text().splitlines()[0]
        assert header.startswith("sample,fastq_1,fastq_2")

    def test_demo_result_json_has_the_standard_envelope(self, tmp_path):
        import arrayexpress_fetch as app

        app.main(["--demo", "--output", str(tmp_path)])
        payload = json.loads((tmp_path / "result.json").read_text())
        for key in ("skill", "version", "summary", "data"):
            assert key in payload
        assert payload["skill"] == "arrayexpress-fetch"

    def test_demo_report_names_the_accession(self, tmp_path):
        import arrayexpress_fetch as app

        app.main(["--demo", "--output", str(tmp_path)])
        assert DEMO_ACCESSION in (tmp_path / "report.md").read_text()

    def test_demo_metadata_table_carries_the_species(self, tmp_path):
        import arrayexpress_fetch as app

        app.main(["--demo", "--output", str(tmp_path)])
        text = (tmp_path / "tables" / "metadata.tsv").read_text()
        assert "Rattus norvegicus" in text


class TestSdrfAndSamplesheet:
    """The SDRF-driven paths — the skill's main value, and the ones §0b blocked."""

    def test_get_sdrf_text_reads_the_attached_file(self):
        import arrayexpress_fetch_api as api

        sdrf = DEMO_SDRF.read_bytes()
        # file_url() resolves /info first, so stub it rather than http_get alone.
        with patch.object(api, "http_get", return_value=sdrf), \
             patch.object(api, "file_url", return_value="https://example/x.sdrf.txt"), \
             patch.object(api, "list_files",
                          return_value=[{"path": f"{DEMO_ACCESSION}.sdrf.txt"}]):
            text = api.get_sdrf_text(DEMO_ACCESSION)
        assert text.startswith("Source Name\t")

    def test_get_sdrf_text_errors_when_no_sdrf_is_attached(self):
        import arrayexpress_fetch_api as api

        with patch.object(api, "list_files", return_value=[{"path": "readme.txt"}]):
            with pytest.raises(SystemExit):
                api.get_sdrf_text(DEMO_ACCESSION)

    def test_finalize_sample_ids_cleans_in_place(self):
        """Cleaning is finalize_sample_ids' job, not write_samplesheet's."""
        import arrayexpress_fetch_api as api

        rows = [{"sample": "Sample 1", "fastq_1": "a_R1.gz", "fastq_2": "a_R2.gz"}]
        api.finalize_sample_ids(rows)
        assert rows[0]["sample"] == "Sample_1"

    def test_finalize_sample_ids_refuses_a_collision(self):
        """Two distinct samples normalising to one id would be silently pooled."""
        import arrayexpress_fetch_api as api

        rows = [{"sample": "Sample 1"}, {"sample": "Sample/1"}]
        with pytest.raises(SystemExit):
            api.finalize_sample_ids(rows)

    def test_bulk_samplesheet_adds_strandedness(self, tmp_path):
        import arrayexpress_fetch_api as api

        rows = [{"sample": "S1", "fastq_1": "a_R1.gz", "fastq_2": "a_R2.gz"}]
        out = tmp_path / "s.csv"
        api.write_samplesheet(rows, out, "bulk", strandedness="reverse")
        text = out.read_text()
        assert "strandedness" in text.splitlines()[0]
        assert "reverse" in text

    def test_fastq_pair_prefers_r1_r2(self):
        import arrayexpress_fetch_api as api

        urls = ["ftp://h/x_R2.fastq.gz", "ftp://h/x_R1.fastq.gz"]
        r1, r2 = api.fastq_pair(urls)
        assert r1.endswith("_R1.fastq.gz") and r2.endswith("_R2.fastq.gz")


# --------------------------------------------------------------------------
# Safety
# --------------------------------------------------------------------------


class TestFastqUriValues:
    def test_a_bare_filename_uri_falls_back_to_ena(self, tmp_path, capsys):
        """Comment[FASTQ_URI] is submitter free text. A bare filename is not a URL:
        it is reported and ignored, so the run's ENA links are used instead."""
        import argparse

        import arrayexpress_fetch_api as api

        sdrf = ("Source Name\tComment[ENA_RUN]\tComment[FASTQ_URI]\n"
                "S1\tERR1\tS1_1.fastq.gz\n"
                "S1\tERR1\tS1_2.fastq.gz\n")
        ena = ["https://ftp.sra.ebi.ac.uk/vol1/ERR1_1.fastq.gz",
               "https://ftp.sra.ebi.ac.uk/vol1/ERR1_2.fastq.gz"]
        out = tmp_path / "samplesheet.csv"
        args = argparse.Namespace(accession="E-MTAB-1", assay="bulk", strandedness="auto",
                                  out=str(out), local_dir=None, fastq_dir=None,
                                  fastq_naming=None, read_map=None)
        with patch.object(api, "get_sdrf_text", return_value=sdrf), \
                patch.object(api, "ena_run_fastq", return_value=ena) as fallback:
            api.cmd_samplesheet(args)
        fallback.assert_called_once_with("ERR1")
        row = next(csv.DictReader(out.open()))
        assert (row["fastq_1"], row["fastq_2"]) == tuple(ena)
        assert "S1_1.fastq.gz" in capsys.readouterr().err


class TestSafety:
    @pytest.mark.parametrize("hostile", ["../escaped.txt", "a/../../escaped.txt", "/tmp/escaped.txt"])
    def test_download_refuses_a_path_that_escapes_the_output_dir(self, tmp_path, hostile):
        """The file path comes from the archive API response, not from the user."""
        import arrayexpress_fetch_api as api

        out = tmp_path / "out"
        with patch.object(api.urllib.request, "urlopen",
                          side_effect=AssertionError("must refuse before any request")), \
                patch.object(api, "study_info", return_value={}):
            with pytest.raises(SystemExit, match="outside"):
                api.download_file("S-TEST1", hostile, str(out))
        assert not (tmp_path / "escaped.txt").exists()

    def test_report_carries_the_disclaimer(self, tmp_path):
        import arrayexpress_fetch as app
        from clawbio.common.report import DISCLAIMER

        app.main(["--demo", "--output", str(tmp_path)])
        assert DISCLAIMER.strip()[:40] in (tmp_path / "report.md").read_text()

    def test_warns_before_overwriting(self, tmp_path, capsys):
        import arrayexpress_fetch as app

        app.main(["--demo", "--output", str(tmp_path)])
        capsys.readouterr()
        app.main(["--demo", "--output", str(tmp_path)])
        assert "overwrit" in capsys.readouterr().err.lower()

    def test_skill_sends_no_user_data(self):
        """Only public accessions and query terms go out. Nothing local."""
        import arrayexpress_fetch_api as api

        url = urllib.parse.urlsplit(api.API)
        assert (url.scheme, url.hostname) == ("https", "www.ebi.ac.uk")
        assert "ArrayExpress" == api.COLLECTION


# --------------------------------------------------------------------------
# Output contract (copied from scaffold_skill.py — no shared base class exists)
# --------------------------------------------------------------------------


def _parse_output_contract(skill_md):
    """Extract files promised in the SKILL.md '## Output Structure' tree.

    Returns output-relative file paths. Directory lines, and any entry whose
    inline comment contains 'optional', are skipped. Returns [] when there is no
    parseable section, so skills without the section are simply not gated.
    """
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
            continue  # the root output_directory/ line
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
    """Every artifact promised in SKILL.md '## Output Structure' must be produced.

    Guards against doc/code drift: a skill documenting an output it never writes.
    Mark conditional artifacts with '(optional)' in the tree comment to exempt them.
    """

    def test_documented_outputs_are_produced(self, tmp_path):
        promised = _parse_output_contract(SKILL_DIR / "SKILL.md")
        if not promised:
            pytest.skip("No parseable '## Output Structure' section in SKILL.md")
        result = subprocess.run(
            [sys.executable, str(SCRIPT), "--demo", "--output", str(tmp_path)],
            capture_output=True, text=True
        )
        assert result.returncode == 0, f"demo run failed: {result.stderr}"
        missing = [p for p in promised if not (tmp_path / p).exists()]
        assert not missing, (
            "SKILL.md Output Structure promises artifacts the skill did not "
            "produce: " + ", ".join(missing) + ". Write them, mark them "
            "'(optional)' in the SKILL.md tree, or remove them from the "
            "documented Output Structure."
        )


class TestDeterminism:
    """The report must describe what was produced, not where it was put."""

    def test_report_holds_no_absolute_output_path(self, tmp_path):
        import arrayexpress_fetch as app

        app.main(["--demo", "--output", str(tmp_path)])
        assert str(tmp_path) not in (tmp_path / "report.md").read_text()

    def test_two_output_dirs_give_identical_reports(self, tmp_path):
        import arrayexpress_fetch as app

        a, b = tmp_path / "a", tmp_path / "b"
        app.main(["--demo", "--output", str(a)])
        app.main(["--demo", "--output", str(b)])
        assert (a / "report.md").read_text() == (b / "report.md").read_text()


class TestMetadataTableRowUnit:
    """One row per sample x replicate, not one per SDRF line.

    The demo fixture cannot catch a regression here: E-MTAB-10030 is a 10x
    study, so its reads are packed into read1/read2/index1 column pairs on a
    single SDRF row per sample. A conventional bulk paired-end SDRF puts each
    FASTQ on its own row, which is the shape that exposed the double-counting
    (E-MTAB-5688: 24 SDRF rows, 12 samples). Hence a synthetic SDRF here,
    following the pattern in pride-fetch's tests.
    """

    SDRF = (
        "Source Name\tCharacteristics[organism]\tComment[ENA_RUN]\t"
        "Comment[LIBRARY_LAYOUT]\tComment[FASTQ_URI]\n"
        "SampleA\tHomo sapiens\tERR1\tPAIRED\tftp://ftp.sra.ebi.ac.uk/a_1.fastq.gz\n"
        "SampleA\tHomo sapiens\tERR1\tPAIRED\tftp://ftp.sra.ebi.ac.uk/a_2.fastq.gz\n"
        "SampleB\tHomo sapiens\tERR2\tPAIRED\tftp://ftp.sra.ebi.ac.uk/b_1.fastq.gz\n"
        "SampleB\tHomo sapiens\tERR2\tPAIRED\tftp://ftp.sra.ebi.ac.uk/b_2.fastq.gz\n")

    def _run(self, tmp_path, sdrf):
        import csv as _csv
        from types import SimpleNamespace

        import arrayexpress_fetch_api as api

        out = tmp_path / "metadata.tsv"
        args = SimpleNamespace(accession="E-MTAB-0000", out=str(out))
        with patch.object(api, "get_sdrf_text", return_value=sdrf), \
                patch.object(api, "merge_biosample", side_effect=lambda s, a: a):
            api.cmd_metadata_table(args)
        return list(_csv.DictReader(out.open(), delimiter="\t"))

    def test_paired_end_rows_collapse_to_one_row_per_run(self, tmp_path):
        rows = self._run(tmp_path, self.SDRF)
        assert [r["sample"] for r in rows] == ["SampleA", "SampleB"], (
            "each FASTQ became its own row; group on Comment[ENA_RUN] first")

    def test_a_sample_with_two_runs_keeps_both(self, tmp_path):
        """'sample x replicate' means multi-run samples stay distinguishable."""
        sdrf = self.SDRF.replace(
            "SampleB\tHomo sapiens\tERR2", "SampleA\tHomo sapiens\tERR2")
        rows = self._run(tmp_path, sdrf)
        assert [r["sample"] for r in rows] == ["SampleA", "SampleA"]
        assert len({r["replicate"] for r in rows}) == 2

    def test_biosample_is_looked_up_once_per_group_not_per_fastq(self, tmp_path):
        """Halves the API calls a paired-end study makes; arrayexpress has no
        cache, unlike ena-fetch."""
        from types import SimpleNamespace

        import arrayexpress_fetch_api as api

        sdrf = self.SDRF.replace("Comment[ENA_RUN]", "Comment[BioSD_SAMPLE]") \
                        .replace("ERR1", "SAMEA1").replace("ERR2", "SAMEA2")
        args = SimpleNamespace(accession="E-MTAB-0000", out=str(tmp_path / "m.tsv"))
        with patch.object(api, "get_sdrf_text", return_value=sdrf), \
                patch.object(api, "merge_biosample",
                             side_effect=lambda s, a: a) as merge:
            api.cmd_metadata_table(args)
        assert merge.call_count == 2, f"{merge.call_count} lookups for 2 samples"


class TestDownloadScriptRoutesToEna:
    """ArrayExpress brokers sequencing reads to ENA, so this skill does not
    emit a FASTQ download script -- but it must say so, loudly. It previously
    died inside argparse with no message at all."""

    def test_it_names_ena_fetch_and_exits_non_zero(self, tmp_path):
        proc = subprocess.run(
            [sys.executable, str(SCRIPT), "--command", "download-script",
             "--accession", DEMO_ACCESSION, "--output", str(tmp_path)],
            capture_output=True, text=True)
        assert proc.returncode != 0
        assert "ena-fetch" in proc.stderr, proc.stderr
        assert "sra-tools" in proc.stderr

    def test_it_does_not_advertise_slurm_or_exec_flags(self):
        """--no-slurm/--run/--submit existed only for the emitter that is not
        here. --run and --submit were never read by anything."""
        proc = subprocess.run([sys.executable, str(SCRIPT), "--help"],
                              capture_output=True, text=True)
        for flag in ("--no-slurm", "--run", "--submit"):
            assert flag not in proc.stdout, f"{flag} is still advertised"


# ── Demo output snapshot ──────────────────────────────────────────────────────
# Fails when a code change alters any demo number. If the change is intended,
# run `pytest skills/arrayexpress-fetch/tests --inline-snapshot=fix`; it rewrites the
# snapshot() values below, and the PR diff shows which numbers moved.

@pytest.mark.demo_snapshot
def test_demo_output_matches_snapshot(tmp_path):
    assert demo_runner.demo_outputs("arrayexpress-fetch", tmp_path) == snapshot({'report.md': """\
# ArrayExpress report — E-MTAB-10030

Source: [EMBL-EBI ArrayExpress](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-10030)

## metadata

```
accession   : E-MTAB-10030
title       : Single-cell RNA-seq of rat microglia in monoculture and in coculture with neurons and astrocytes
release     : 2021-03-23
organism    : Rattus norvegicus

description : Microglia, brain resident macrophages, require instruction from the central nervous system microenvironment to maintain their identity, morphology, and to regulate inflammatory responses. We investigated the heterogeneity of response of microglia to the presence of neurons and astrocytes by performing single-cell sequencing of microglia in both monoculture, and in coculture with neurons and astrocytes.

files       : 20 total  idf=1  raw=18  sdrf=1
```

## files

```
20 file(s) for E-MTAB-10030:

raw            2608854551  B1_S2_R1.fastq.gz
raw            2823794346  B1_S2_R2.fastq.gz
raw             427160570  B1_S2_I1.fastq.gz
raw            4991718428  B1_S4_R1.fastq.gz
raw            5216968590  B1_S4_R2.fastq.gz
raw             777240263  B1_S4_I1.fastq.gz
raw            1483407903  B2_S2_R1.fastq.gz
raw            1719034285  B2_S2_R2.fastq.gz
raw             247703776  B2_S2_I1.fastq.gz
raw            1649147648  B1_S1_R1.fastq.gz
raw            1806913320  B1_S1_R2.fastq.gz
raw             277961552  B1_S1_I1.fastq.gz
raw            8960636020  B1_S3_R1.fastq.gz
raw            9299329622  B1_S3_R2.fastq.gz
raw            1378440463  B1_S3_I1.fastq.gz
raw            1915465169  B2_S1_R1.fastq.gz
raw            2252667303  B2_S1_R2.fastq.gz
raw             315655053  B2_S1_I1.fastq.gz
idf                  5544  E-MTAB-10030.idf.txt
sdrf                 7242  E-MTAB-10030.sdrf.txt
```

## sdrf

```
Source Name	Comment[ENA_SAMPLE]	Comment[BioSD_SAMPLE]	Characteristics[organism]	Characteristics[strain]	Characteristics[age]	Unit[time unit]	Term Source REF	Term Accession Number	Characteristics[developmental stage]	Characteristics[sex]	Characteristics[genotype]	Characteristics[organism part]	Characteristics[individual]	Characteristics[cell type]	Characteristics[growth condition]	Characteristics[biological sample]	Characteristics[culture batch]	Characteristics[sequencing batch]	Material Type	Protocol REF	Protocol REF	Protocol REF	Extract Name	Comment[LIBRARY_LAYOUT]	Comment[LIBRARY_SELECTION]	Comment[LIBRARY_SOURCE]	Comment[LIBRARY_STRAND]	Comment[LIBRARY_STRATEGY]	Comment[cdna read]	Comment[cdna read offset]	Comment[cdna read size]	Comment[cell barcode offset]	Comment[cell barcode read]	Comment[cell barcode size]	Comment[end bias]	Comment[input molecule]	Comment[library construction]	Comment[primer]	Comment[sample barcode offset]	Comment[sample barcode read]	Comment[sample barcode size]	Comment[single cell isolation]	Comment[spike in]	Comment[umi barcode offset]	Comment[umi barcode read]	Comment[umi barcode size]	Protocol REF	Performer	Assay Name	Technology Type	Comment[ENA_EXPERIMENT]	Scan Name	Comment[SUBMITTED_FILE_NAME]	Comment[ENA_RUN]	Comment[BAM_URI]	Comment[read1 file]	Comment[FASTQ_URI]	Comment[read2 file]	Comment[FASTQ_URI]	Comment[index1 file]	Comment[FASTQ_URI]	Factor Value[growth condition]
Sample 2	ERS5592419	SAMEA7845252	Rattus norvegicus	Sprague Dawley	0 to 2	day	EFO	UO_0000033	postnatal	mixed	wild type genotype	neocortex	pooled	microglial cell	co-culture with neurons and astrocytes	2	1	1	cell	P-MTAB-105596	P-MTAB-105597	P-MTAB-105598	Sample 2	PAIRED	Oligo-dT	TRANSCRIPTOMIC SINGLE CELL	second strand	ssRNA-seq	read 2	0	91	0	read 1	16	3 prime tag	polyA RNA	10xV3	oligo-dT	0	index 1	8	10x	none	16	read 1	12	P-MTAB-105599	Edinburgh Clinical Research Facility	Sample 2	sequencing assay	ERX4970960	Sample_2.bam	Sample_2.bam	ERR5166018	ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR516/ERR5166018/Sample_2.bam	B1_S2_R1.fastq.gz	ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10030/B1_S2_R1.fastq.gz	B1_S2_R2.fastq.gz	ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10030/B1_S2_R2.fastq.gz	B1_S2_I1.fastq.gz	ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10030/B1_S2_I1.fastq.gz	co-culture with neurons and astrocytes
Sample 4	ERS5592420	SAMEA7845253	Rattus norvegicus	Sprague Dawley	0 to 2	day	EFO	UO_0000033	postnatal	mixed	wild type genotype	neocortex	pooled	microglial cell	co-culture with neurons and astrocytes	4	2	1	cell	P-MTAB-105596	P-MTAB-105597	P-MTAB-105598	Sample 4	PAIRED	Oligo-dT	TRANSCRIPTOMIC SINGLE CELL	second strand	ssRNA-seq	read 2	0	91	0	read 1	16	3 prime tag	polyA RNA	10xV3	oligo-dT	0	index 1	8	10x	none	16	read 1	12	P-MTAB-105599	Edinburgh Clinical Research Facility	Sample 4	sequencing assay	ERX4970961	Sample_4.bam	Sample_4.bam	ERR5166019	ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR516/ERR5166019/Sample_4.bam	B1_S4_R1.fastq.gz	ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10030/B1_S4_R1.fastq.gz	B1_S4_R2.fastq.gz	ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10030/B1_S4_R2.fastq.gz	B1_S4_I1.fastq.gz	ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10030/B1_S4_I1.fastq.gz	co-culture with neurons and astrocytes
Sample 6	ERS5592421	SAMEA7845254	Rattus norvegicus	Sprague Dawley	0 to 2	day	EFO	UO_0000033	postnatal	mixed	wild type genotype	neocortex	pooled	microglial cell	co-culture with neurons and astrocytes	2	1	2	cell	P-MTAB-105596	P-MTAB-105597	P-MTAB-105598	Sample 6	PAIRED	Oligo-dT	TRANSCRIPTOMIC SINGLE CELL	second strand	ssRNA-seq	read 2	0	91	0	read 1	16	3 prime tag	polyA RNA	10xV3	oligo-dT	0	index 1	8	10x	none	16	read 1	12	P-MTAB-105599	Edinburgh Clinical Research Facility	Sample 6	sequencing assay	ERX4970962	Sample_6.bam	Sample_6.bam	ERR5166020	ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR516/ERR5166020/Sample_6.bam	B2_S2_R1.fastq.gz	ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10030/B2_S2_R1.fastq.gz	B2_S2_R2.fastq.gz	ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10030/B2_S2_R2.fastq.gz	B2_S2_I1.fastq.gz	ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10030/B2_S2_I1.fastq.gz	co-culture with neurons and astrocytes
Sample 1	ERS5592422	SAMEA7845255	Rattus norvegicus	Sprague Dawley	0 to 2	day	EFO	UO_0000033	postnatal	mixed	wild type genotype	neocortex	pooled	microglial cell	monoculture	1	1	1	cell	P-MTAB-105596	P-MTAB-105597	P-MTAB-105598	Sample 1	PAIRED	Oligo-dT	TRANSCRIPTOMIC SINGLE CELL	second strand	ssRNA-seq	read 2	0	91	0	read 1	16	3 prime tag	polyA RNA	10xV3	oligo-dT	0	index 1	8	10x	none	16	read 1	12	P-MTAB-105599	Edinburgh Clinical Research Facility	Sample 1	sequencing assay	ERX4970963	Sample_1.bam	Sample_1.bam	ERR5166021	ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR516/ERR5166021/Sample_1.bam	B1_S1_R1.fastq.gz	ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10030/B1_S1_R1.fastq.gz	B1_S1_R2.fastq.gz	ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10030/B1_S1_R2.fastq.gz	B1_S1_I1.fastq.gz	ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10030/B1_S1_I1.fastq.gz	monoculture
Sample 3	ERS5592423	SAMEA7845256	Rattus norvegicus	Sprague Dawley	0 to 2	day	EFO	UO_0000033	postnatal	mixed	wild type genotype	neocortex	pooled	microglial cell	monoculture	3	2	1	cell	P-MTAB-105596	P-MTAB-105597	P-MTAB-105598	Sample 3	PAIRED	Oligo-dT	TRANSCRIPTOMIC SINGLE CELL	second strand	ssRNA-seq	read 2	0	91	0	read 1	16	3 prime tag	polyA RNA	10xV3	oligo-dT	0	index 1	8	10x	none	16	read 1	12	P-MTAB-105599	Edinburgh Clinical Research Facility	Sample 3	sequencing assay	ERX4970964	Sample_3.bam	Sample_3.bam	ERR5166022	ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR516/ERR5166022/Sample_3.bam	B1_S3_R1.fastq.gz	ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10030/B1_S3_R1.fastq.gz	B1_S3_R2.fastq.gz	ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10030/B1_S3_R2.fastq.gz	B1_S3_I1.fastq.gz	ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10030/B1_S3_I1.fastq.gz	monoculture
Sample 5	ERS5592424	SAMEA7845257	Rattus norvegicus	Sprague Dawley	0 to 2	day	EFO	UO_0000033	postnatal	mixed	wild type genotype	neocortex	pooled	microglial cell	monoculture	1	1	2	cell	P-MTAB-105596	P-MTAB-105597	P-MTAB-105598	Sample 5	PAIRED	Oligo-dT	TRANSCRIPTOMIC SINGLE CELL	second strand	ssRNA-seq	read 2	0	91	0	read 1	16	3 prime tag	polyA RNA	10xV3	oligo-dT	0	index 1	8	10x	none	16	read 1	12	P-MTAB-105599	Edinburgh Clinical Research Facility	Sample 5	sequencing assay	ERX4970965	Sample_5.bam	Sample_5.bam	ERR5166023	ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR516/ERR5166023/Sample_5.bam	B2_S1_R1.fastq.gz	ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10030/B2_S1_R1.fastq.gz	B2_S1_R2.fastq.gz	ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10030/B2_S1_R2.fastq.gz	B2_S1_I1.fastq.gz	ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10030/B2_S1_I1.fastq.gz	monoculture
```

## metadata-table

```
Wrote 6 sample x replicate row(s) x 15 column(s) for E-MTAB-10030 to <output>/tables/metadata.tsv
```

## samplesheet

```
warning: sample 'Sample 2' -> 'Sample_2'
warning: sample 'Sample 4' -> 'Sample_4'
warning: sample 'Sample 6' -> 'Sample_6'
warning: sample 'Sample 1' -> 'Sample_1'
warning: sample 'Sample 3' -> 'Sample_3'
warning: sample 'Sample 5' -> 'Sample_5'
Wrote 6 row(s) for 6 sample(s) to <output>/samplesheet.csv (nf-core/scrnaseq)
```

## search

```
223 total hits; showing 5:

E-GEOD-43808     Transcriptome signature of mouse embryonic stem cell-derived microglia is closely related to primary cultured microglia
E-GEOD-45384     Effects of Th1- and Th17-derived factors on microglia.
E-GEOD-18602     Microglia in ischemic brain injury
E-GEOD-57181     In vivo response of microglia to Aβ and prostaglandin-E2 EP2 receptor signaling
E-GEOD-68376     Interleukin-1 signaling is involved in local self-renewal and maintenance of microglia
```

---

*ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions.*
""", 'result.json': """\
{
  "status": "ok",
  "ok": true,
  "skill": "arrayexpress-fetch",
  "version": "0.1.0",
  "completed_at": "<TIMESTAMP>",
  "input_checksum": "",
  "datasets": {},
  "summary": {
    "accession": "E-MTAB-10030",
    "commands": [
      "metadata",
      "files",
      "sdrf",
      "metadata-table",
      "samplesheet",
      "search"
    ],
    "demo": true
  },
  "data": {
    "sections": {
      "metadata": "accession   : E-MTAB-10030\\ntitle       : Single-cell RNA-seq of rat microglia in monoculture and in coculture with neurons and astrocytes\\nrelease     : 2021-03-23\\norganism    : Rattus norvegicus\\n\\ndescription : Microglia, brain resident macrophages, require instruction from the central nervous system microenvironment to maintain their identity, morphology, and to regulate inflammatory responses. We investigated the heterogeneity of response of microglia to the presence of neurons and astrocytes by performing single-cell sequencing of microglia in both monoculture, and in coculture with neurons and astrocytes.\\n\\nfiles       : 20 total  idf=1  raw=18  sdrf=1\\n",
      "files": "20 file(s) for E-MTAB-10030:\\n\\nraw            2608854551  B1_S2_R1.fastq.gz\\nraw            2823794346  B1_S2_R2.fastq.gz\\nraw             427160570  B1_S2_I1.fastq.gz\\nraw            4991718428  B1_S4_R1.fastq.gz\\nraw            5216968590  B1_S4_R2.fastq.gz\\nraw             777240263  B1_S4_I1.fastq.gz\\nraw            1483407903  B2_S2_R1.fastq.gz\\nraw            1719034285  B2_S2_R2.fastq.gz\\nraw             247703776  B2_S2_I1.fastq.gz\\nraw            1649147648  B1_S1_R1.fastq.gz\\nraw            1806913320  B1_S1_R2.fastq.gz\\nraw             277961552  B1_S1_I1.fastq.gz\\nraw            8960636020  B1_S3_R1.fastq.gz\\nraw            9299329622  B1_S3_R2.fastq.gz\\nraw            1378440463  B1_S3_I1.fastq.gz\\nraw            1915465169  B2_S1_R1.fastq.gz\\nraw            2252667303  B2_S1_R2.fastq.gz\\nraw             315655053  B2_S1_I1.fastq.gz\\nidf                  5544  E-MTAB-10030.idf.txt\\nsdrf                 7242  E-MTAB-10030.sdrf.txt\\n",
      "sdrf": "Source Name\\tComment[ENA_SAMPLE]\\tComment[BioSD_SAMPLE]\\tCharacteristics[organism]\\tCharacteristics[strain]\\tCharacteristics[age]\\tUnit[time unit]\\tTerm Source REF\\tTerm Accession Number\\tCharacteristics[developmental stage]\\tCharacteristics[sex]\\tCharacteristics[genotype]\\tCharacteristics[organism part]\\tCharacteristics[individual]\\tCharacteristics[cell type]\\tCharacteristics[growth condition]\\tCharacteristics[biological sample]\\tCharacteristics[culture batch]\\tCharacteristics[sequencing batch]\\tMaterial Type\\tProtocol REF\\tProtocol REF\\tProtocol REF\\tExtract Name\\tComment[LIBRARY_LAYOUT]\\tComment[LIBRARY_SELECTION]\\tComment[LIBRARY_SOURCE]\\tComment[LIBRARY_STRAND]\\tComment[LIBRARY_STRATEGY]\\tComment[cdna read]\\tComment[cdna read offset]\\tComment[cdna read size]\\tComment[cell barcode offset]\\tComment[cell barcode read]\\tComment[cell barcode size]\\tComment[end bias]\\tComment[input molecule]\\tComment[library construction]\\tComment[primer]\\tComment[sample barcode offset]\\tComment[sample barcode read]\\tComment[sample barcode size]\\tComment[single cell isolation]\\tComment[spike in]\\tComment[umi barcode offset]\\tComment[umi barcode read]\\tComment[umi barcode size]\\tProtocol REF\\tPerformer\\tAssay Name\\tTechnology Type\\tComment[ENA_EXPERIMENT]\\tScan Name\\tComment[SUBMITTED_FILE_NAME]\\tComment[ENA_RUN]\\tComment[BAM_URI]\\tComment[read1 file]\\tComment[FASTQ_URI]\\tComment[read2 file]\\tComment[FASTQ_URI]\\tComment[index1 file]\\tComment[FASTQ_URI]\\tFactor Value[growth condition]\\nSample 2\\tERS5592419\\tSAMEA7845252\\tRattus norvegicus\\tSprague Dawley\\t0 to 2\\tday\\tEFO\\tUO_0000033\\tpostnatal\\tmixed\\twild type genotype\\tneocortex\\tpooled\\tmicroglial cell\\tco-culture with neurons and astrocytes\\t2\\t1\\t1\\tcell\\tP-MTAB-105596\\tP-MTAB-105597\\tP-MTAB-105598\\tSample 2\\tPAIRED\\tOligo-dT\\tTRANSCRIPTOMIC SINGLE CELL\\tsecond strand\\tssRNA-seq\\tread 2\\t0\\t91\\t0\\tread 1\\t16\\t3 prime tag\\tpolyA RNA\\t10xV3\\toligo-dT\\t0\\tindex 1\\t8\\t10x\\tnone\\t16\\tread 1\\t12\\tP-MTAB-105599\\tEdinburgh Clinical Research Facility\\tSample 2\\tsequencing assay\\tERX4970960\\tSample_2.bam\\tSample_2.bam\\tERR5166018\\tftp://ftp.sra.ebi.ac.uk/vol1/run/ERR516/ERR5166018/Sample_2.bam\\tB1_S2_R1.fastq.gz\\tftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10030/B1_S2_R1.fastq.gz\\tB1_S2_R2.fastq.gz\\tftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10030/B1_S2_R2.fastq.gz\\tB1_S2_I1.fastq.gz\\tftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10030/B1_S2_I1.fastq.gz\\tco-culture with neurons and astrocytes\\nSample 4\\tERS5592420\\tSAMEA7845253\\tRattus norvegicus\\tSprague Dawley\\t0 to 2\\tday\\tEFO\\tUO_0000033\\tpostnatal\\tmixed\\twild type genotype\\tneocortex\\tpooled\\tmicroglial cell\\tco-culture with neurons and astrocytes\\t4\\t2\\t1\\tcell\\tP-MTAB-105596\\tP-MTAB-105597\\tP-MTAB-105598\\tSample 4\\tPAIRED\\tOligo-dT\\tTRANSCRIPTOMIC SINGLE CELL\\tsecond strand\\tssRNA-seq\\tread 2\\t0\\t91\\t0\\tread 1\\t16\\t3 prime tag\\tpolyA RNA\\t10xV3\\toligo-dT\\t0\\tindex 1\\t8\\t10x\\tnone\\t16\\tread 1\\t12\\tP-MTAB-105599\\tEdinburgh Clinical Research Facility\\tSample 4\\tsequencing assay\\tERX4970961\\tSample_4.bam\\tSample_4.bam\\tERR5166019\\tftp://ftp.sra.ebi.ac.uk/vol1/run/ERR516/ERR5166019/Sample_4.bam\\tB1_S4_R1.fastq.gz\\tftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10030/B1_S4_R1.fastq.gz\\tB1_S4_R2.fastq.gz\\tftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10030/B1_S4_R2.fastq.gz\\tB1_S4_I1.fastq.gz\\tftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10030/B1_S4_I1.fastq.gz\\tco-culture with neurons and astrocytes\\nSample 6\\tERS5592421\\tSAMEA7845254\\tRattus norvegicus\\tSprague Dawley\\t0 to 2\\tday\\tEFO\\tUO_0000033\\tpostnatal\\tmixed\\twild type genotype\\tneocortex\\tpooled\\tmicroglial cell\\tco-culture with neurons and astrocytes\\t2\\t1\\t2\\tcell\\tP-MTAB-105596\\tP-MTAB-105597\\tP-MTAB-105598\\tSample 6\\tPAIRED\\tOligo-dT\\tTRANSCRIPTOMIC SINGLE CELL\\tsecond strand\\tssRNA-seq\\tread 2\\t0\\t91\\t0\\tread 1\\t16\\t3 prime tag\\tpolyA RNA\\t10xV3\\toligo-dT\\t0\\tindex 1\\t8\\t10x\\tnone\\t16\\tread 1\\t12\\tP-MTAB-105599\\tEdinburgh Clinical Research Facility\\tSample 6\\tsequencing assay\\tERX4970962\\tSample_6.bam\\tSample_6.bam\\tERR5166020\\tftp://ftp.sra.ebi.ac.uk/vol1/run/ERR516/ERR5166020/Sample_6.bam\\tB2_S2_R1.fastq.gz\\tftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10030/B2_S2_R1.fastq.gz\\tB2_S2_R2.fastq.gz\\tftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10030/B2_S2_R2.fastq.gz\\tB2_S2_I1.fastq.gz\\tftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10030/B2_S2_I1.fastq.gz\\tco-culture with neurons and astrocytes\\nSample 1\\tERS5592422\\tSAMEA7845255\\tRattus norvegicus\\tSprague Dawley\\t0 to 2\\tday\\tEFO\\tUO_0000033\\tpostnatal\\tmixed\\twild type genotype\\tneocortex\\tpooled\\tmicroglial cell\\tmonoculture\\t1\\t1\\t1\\tcell\\tP-MTAB-105596\\tP-MTAB-105597\\tP-MTAB-105598\\tSample 1\\tPAIRED\\tOligo-dT\\tTRANSCRIPTOMIC SINGLE CELL\\tsecond strand\\tssRNA-seq\\tread 2\\t0\\t91\\t0\\tread 1\\t16\\t3 prime tag\\tpolyA RNA\\t10xV3\\toligo-dT\\t0\\tindex 1\\t8\\t10x\\tnone\\t16\\tread 1\\t12\\tP-MTAB-105599\\tEdinburgh Clinical Research Facility\\tSample 1\\tsequencing assay\\tERX4970963\\tSample_1.bam\\tSample_1.bam\\tERR5166021\\tftp://ftp.sra.ebi.ac.uk/vol1/run/ERR516/ERR5166021/Sample_1.bam\\tB1_S1_R1.fastq.gz\\tftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10030/B1_S1_R1.fastq.gz\\tB1_S1_R2.fastq.gz\\tftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10030/B1_S1_R2.fastq.gz\\tB1_S1_I1.fastq.gz\\tftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10030/B1_S1_I1.fastq.gz\\tmonoculture\\nSample 3\\tERS5592423\\tSAMEA7845256\\tRattus norvegicus\\tSprague Dawley\\t0 to 2\\tday\\tEFO\\tUO_0000033\\tpostnatal\\tmixed\\twild type genotype\\tneocortex\\tpooled\\tmicroglial cell\\tmonoculture\\t3\\t2\\t1\\tcell\\tP-MTAB-105596\\tP-MTAB-105597\\tP-MTAB-105598\\tSample 3\\tPAIRED\\tOligo-dT\\tTRANSCRIPTOMIC SINGLE CELL\\tsecond strand\\tssRNA-seq\\tread 2\\t0\\t91\\t0\\tread 1\\t16\\t3 prime tag\\tpolyA RNA\\t10xV3\\toligo-dT\\t0\\tindex 1\\t8\\t10x\\tnone\\t16\\tread 1\\t12\\tP-MTAB-105599\\tEdinburgh Clinical Research Facility\\tSample 3\\tsequencing assay\\tERX4970964\\tSample_3.bam\\tSample_3.bam\\tERR5166022\\tftp://ftp.sra.ebi.ac.uk/vol1/run/ERR516/ERR5166022/Sample_3.bam\\tB1_S3_R1.fastq.gz\\tftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10030/B1_S3_R1.fastq.gz\\tB1_S3_R2.fastq.gz\\tftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10030/B1_S3_R2.fastq.gz\\tB1_S3_I1.fastq.gz\\tftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10030/B1_S3_I1.fastq.gz\\tmonoculture\\nSample 5\\tERS5592424\\tSAMEA7845257\\tRattus norvegicus\\tSprague Dawley\\t0 to 2\\tday\\tEFO\\tUO_0000033\\tpostnatal\\tmixed\\twild type genotype\\tneocortex\\tpooled\\tmicroglial cell\\tmonoculture\\t1\\t1\\t2\\tcell\\tP-MTAB-105596\\tP-MTAB-105597\\tP-MTAB-105598\\tSample 5\\tPAIRED\\tOligo-dT\\tTRANSCRIPTOMIC SINGLE CELL\\tsecond strand\\tssRNA-seq\\tread 2\\t0\\t91\\t0\\tread 1\\t16\\t3 prime tag\\tpolyA RNA\\t10xV3\\toligo-dT\\t0\\tindex 1\\t8\\t10x\\tnone\\t16\\tread 1\\t12\\tP-MTAB-105599\\tEdinburgh Clinical Research Facility\\tSample 5\\tsequencing assay\\tERX4970965\\tSample_5.bam\\tSample_5.bam\\tERR5166023\\tftp://ftp.sra.ebi.ac.uk/vol1/run/ERR516/ERR5166023/Sample_5.bam\\tB2_S1_R1.fastq.gz\\tftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10030/B2_S1_R1.fastq.gz\\tB2_S1_R2.fastq.gz\\tftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10030/B2_S1_R2.fastq.gz\\tB2_S1_I1.fastq.gz\\tftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10030/B2_S1_I1.fastq.gz\\tmonoculture\\n",
      "metadata-table": "Wrote 6 sample x replicate row(s) x 15 column(s) for E-MTAB-10030 to <output>/tables/metadata.tsv\\n",
      "samplesheet": "warning: sample 'Sample 2' -> 'Sample_2'\\nwarning: sample 'Sample 4' -> 'Sample_4'\\nwarning: sample 'Sample 6' -> 'Sample_6'\\nwarning: sample 'Sample 1' -> 'Sample_1'\\nwarning: sample 'Sample 3' -> 'Sample_3'\\nwarning: sample 'Sample 5' -> 'Sample_5'\\nWrote 6 row(s) for 6 sample(s) to <output>/samplesheet.csv (nf-core/scrnaseq)\\n",
      "search": "223 total hits; showing 5:\\n\\nE-GEOD-43808     Transcriptome signature of mouse embryonic stem cell-derived microglia is closely related to primary cultured microglia\\nE-GEOD-45384     Effects of Th1- and Th17-derived factors on microglia.\\nE-GEOD-18602     Microglia in ischemic brain injury\\nE-GEOD-57181     In vivo response of microglia to A\\u03b2 and prostaglandin-E2 EP2 receptor signaling\\nE-GEOD-68376     Interleukin-1 signaling is involved in local self-renewal and maintenance of microglia\\n"
    }
  }
}\
""", 'samplesheet.csv': """\
sample,fastq_1,fastq_2
Sample_1,https://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10030/B1_S1_R1.fastq.gz,https://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10030/B1_S1_R2.fastq.gz
Sample_2,https://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10030/B1_S2_R1.fastq.gz,https://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10030/B1_S2_R2.fastq.gz
Sample_3,https://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10030/B1_S3_R1.fastq.gz,https://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10030/B1_S3_R2.fastq.gz
Sample_4,https://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10030/B1_S4_R1.fastq.gz,https://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10030/B1_S4_R2.fastq.gz
Sample_5,https://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10030/B2_S1_R1.fastq.gz,https://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10030/B2_S1_R2.fastq.gz
Sample_6,https://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10030/B2_S2_R1.fastq.gz,https://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10030/B2_S2_R2.fastq.gz
""", 'tables/metadata.tsv': """\
sample	replicate	species	sex	age	condition	genotype	treatment	tissue	individual	cell_type	growth_condition	biological_sample	culture_batch	sequencing_batch
Sample 2	ERR5166018	Rattus norvegicus	mixed	0 to 2	NA	wild type genotype	NA	neocortex	pooled	microglial cell	co-culture with neurons and astrocytes	2	1	1
Sample 4	ERR5166019	Rattus norvegicus	mixed	0 to 2	NA	wild type genotype	NA	neocortex	pooled	microglial cell	co-culture with neurons and astrocytes	4	2	1
Sample 6	ERR5166020	Rattus norvegicus	mixed	0 to 2	NA	wild type genotype	NA	neocortex	pooled	microglial cell	co-culture with neurons and astrocytes	2	1	2
Sample 1	ERR5166021	Rattus norvegicus	mixed	0 to 2	NA	wild type genotype	NA	neocortex	pooled	microglial cell	monoculture	1	1	1
Sample 3	ERR5166022	Rattus norvegicus	mixed	0 to 2	NA	wild type genotype	NA	neocortex	pooled	microglial cell	monoculture	3	2	1
Sample 5	ERR5166023	Rattus norvegicus	mixed	0 to 2	NA	wild type genotype	NA	neocortex	pooled	microglial cell	monoculture	1	1	2
"""})
