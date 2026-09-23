"""Tests for arrayexpress-fetch.

Run with: pytest skills/arrayexpress-fetch/tests/test_arrayexpress_fetch.py -v

No network required: every HTTP call is mocked against the committed demo
fixtures in ../examples/.
"""

import json
import re
import subprocess
import sys
from pathlib import Path
from unittest.mock import patch

import pytest

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


class TestVendoredApi:
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


class TestSafety:
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

        assert api.API.startswith("https://www.ebi.ac.uk")
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
