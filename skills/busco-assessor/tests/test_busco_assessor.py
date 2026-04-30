"""
TDD test suite for busco_assessor.py — all tests must pass without the BUSCO binary.
Run: pytest skills/busco-assessor/tests/ -v
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from unittest.mock import patch

import pytest

SKILL_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SKILL_DIR))

import gzip
from unittest.mock import MagicMock, patch

from busco_assessor import (  # noqa: E402
    DEMO_COMPLETENESS,
    DEMO_FULL_TABLE_ROWS,
    DEMO_LINEAGE,
    DEMO_LIVE_COMPLETENESS,
    DEMO_LIVE_ORGANISM,
    DEMO_LIVE_URL,
    DISCLAIMER,
    MAX_DOWNLOAD_BYTES,
    NCBI_USER_AGENT,
    build_busco_command,
    download_and_decompress,
    infer_lineage,
    make_demo_busco_outputs,
    make_demo_fasta,
    ncbi_lineage_to_busco,
    ncbi_taxonomy_lookup,
    parse_full_table,
    parse_short_summary,
    run_demo,
    run_demo_live,
    validate_args,
    write_reproducibility_bundle,
    write_report,
    write_result_json,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_args(**kwargs) -> argparse.Namespace:
    defaults = dict(
        input=None,
        output="/tmp/busco_test",
        demo=False,
        demo_live=False,
        mode="genome",
        lineage=None,
        auto_lineage=False,
        auto_lineage_euk=False,
        auto_lineage_prok=False,
        augustus=False,
        cpu=4,
        download_path=None,
        organism=None,
    )
    defaults.update(kwargs)
    return argparse.Namespace(**defaults)


# ---------------------------------------------------------------------------
# TestInferLineage
# ---------------------------------------------------------------------------

class TestInferLineage:
    def test_bacteria_hint_gives_auto_lineage_prok(self):
        flag, value = infer_lineage("my bacteria genome")
        assert flag == "--auto-lineage-prok"
        assert value is None

    def test_ecoli_hint_gives_auto_lineage_prok(self):
        flag, value = infer_lineage("E. coli K-12 assembly")
        assert flag == "--auto-lineage-prok"
        assert value is None

    def test_human_hint_gives_primates_odb10(self):
        flag, value = infer_lineage("human assembly hg38")
        assert flag == "--lineage"
        assert value == "primates_odb10"

    def test_drosophila_hint_gives_diptera(self):
        flag, value = infer_lineage("fruit fly transcriptome")
        assert flag == "--lineage"
        assert value == "diptera_odb10"

    def test_unknown_hint_gives_auto_lineage(self):
        flag, value = infer_lineage("mystery organism XJ-9")
        assert flag == "--auto-lineage"
        assert value is None

    def test_case_insensitive(self):
        flag, value = infer_lineage("BACTERIA genome")
        assert flag == "--auto-lineage-prok"
        assert value is None

    def test_archaea_gives_archaea_odb12(self):
        flag, value = infer_lineage("archaea genome assembly")
        assert flag == "--lineage"
        assert value == "archaea_odb12"

    def test_generic_eukaryote_gives_auto_lineage_euk(self):
        flag, value = infer_lineage("eukaryote genome")
        assert flag == "--auto-lineage-euk"
        assert value is None

    def test_mouse_gives_mammalia(self):
        flag, value = infer_lineage("mouse genome Mus musculus")
        assert flag == "--lineage"
        assert value == "mammalia_odb10"

    def test_plant_gives_embryophyta(self):
        flag, value = infer_lineage("Arabidopsis thaliana genome")
        assert flag == "--lineage"
        assert value == "embryophyta_odb10"

    def test_yeast_gives_fungi(self):
        flag, value = infer_lineage("Saccharomyces cerevisiae assembly")
        assert flag == "--lineage"
        assert value == "fungi_odb10"

    def test_bacteriophage_not_routed_to_prok(self):
        # "bacteria" is a substring of "bacteriophage" but NOT a whole word — must not route to prok
        flag, value = infer_lineage("bacteriophage assembly")
        assert flag == "--auto-lineage"
        assert value is None

    def test_rational_not_routed_to_mammalia(self):
        # "rat" is a substring of "rational" — must not trigger mammalia routing
        flag, value = infer_lineage("rational design organism")
        assert flag == "--auto-lineage"
        assert value is None

    def test_transplant_not_routed_to_embryophyta(self):
        # "plant" is a substring of "transplant"
        flag, value = infer_lineage("transplant organ genome")
        assert flag == "--auto-lineage"
        assert value is None


# ---------------------------------------------------------------------------
# TestParseShortSummary
# ---------------------------------------------------------------------------

DEMO_SHORT_SUMMARY_TEXT = """\
# BUSCO version is: 6.0.0
# The lineage dataset is: bacteria_odb12 (Creation date: 2021-09-02, number of species: 3857, number of BUSCOs: 124)
# Summarized benchmarking in ALL lineages: /tmp/busco_run
# BUSCO was run in mode: genome
# Gene predictor used: prodigal

        --------------------------------------------------
        |Results from dataset bacteria_odb12              |
        --------------------------------------------------
        |C:95.2%[S:93.1%,D:2.1%],F:2.3%,M:2.5%,n:124    |
        |118    Complete BUSCOs (C)                       |
        |115    Complete and single-copy BUSCOs (S)       |
        |3      Complete and duplicated BUSCOs (D)        |
        |3      Fragmented BUSCOs (F)                     |
        |3      Missing BUSCOs (M)                        |
        |124    Total BUSCO groups searched               |
        --------------------------------------------------
"""


class TestParseShortSummary:
    def test_parses_completeness_values(self, tmp_path):
        summary = tmp_path / "short_summary.txt"
        summary.write_text(DEMO_SHORT_SUMMARY_TEXT)
        scores = parse_short_summary(summary)
        assert scores["C"] == pytest.approx(95.2)
        assert scores["S"] == pytest.approx(93.1)
        assert scores["D"] == pytest.approx(2.1)
        assert scores["F"] == pytest.approx(2.3)
        assert scores["M"] == pytest.approx(2.5)
        assert scores["n"] == 124

    def test_score_string_roundtrip(self, tmp_path):
        summary = tmp_path / "short_summary.txt"
        summary.write_text(DEMO_SHORT_SUMMARY_TEXT)
        scores = parse_short_summary(summary)
        assert scores["score_string"] == "C:95.2%[S:93.1%,D:2.1%],F:2.3%,M:2.5%,n:124"

    def test_missing_file_raises(self):
        with pytest.raises(FileNotFoundError):
            parse_short_summary(Path("/nonexistent/short_summary.txt"))


# ---------------------------------------------------------------------------
# TestParseFullTable
# ---------------------------------------------------------------------------

DEMO_FULL_TABLE_TSV = """\
# Busco id\tStatus\tSequence\tScore\tLength
1098at2\tComplete\tseq1\t742.3\t312
1099at2\tComplete\tseq1\t698.1\t287
1103at2\tFragmented\tseq2\t341.2\t98
"""


class TestParseFullTable:
    def test_parses_tsv_rows(self, tmp_path):
        tsv = tmp_path / "full_table.tsv"
        tsv.write_text(DEMO_FULL_TABLE_TSV)
        rows = parse_full_table(tsv)
        assert len(rows) == 3
        assert rows[0]["status"] == "Complete"
        assert rows[2]["status"] == "Fragmented"

    def test_skips_comment_lines(self, tmp_path):
        tsv = tmp_path / "full_table.tsv"
        tsv.write_text(DEMO_FULL_TABLE_TSV)
        rows = parse_full_table(tsv)
        for row in rows:
            assert not row["busco_id"].startswith("#")


# ---------------------------------------------------------------------------
# TestValidateArgs
# ---------------------------------------------------------------------------

class TestValidateArgs:
    def test_demo_skips_input_requirement(self):
        args = _make_args(demo=True)
        validate_args(args)  # should not raise

    def test_input_required_without_demo(self):
        args = _make_args(demo=False, input=None)
        with pytest.raises(ValueError, match="--input"):
            validate_args(args)

    def test_conflicting_lineage_flags(self):
        args = _make_args(demo=True, lineage="bacteria_odb12", auto_lineage=True)
        with pytest.raises(ValueError, match="lineage"):
            validate_args(args)

    def test_proteins_mode_with_fna_warns(self, tmp_path, capsys):
        fna = tmp_path / "assembly.fna"
        fna.write_text(">seq1\nACGT\n")
        args = _make_args(demo=False, input=str(fna), mode="proteins")
        validate_args(args)
        captured = capsys.readouterr()
        assert "WARNING" in captured.out or "WARNING" in captured.err


# ---------------------------------------------------------------------------
# TestMakeDemoOutputs
# ---------------------------------------------------------------------------

class TestMakeDemoOutputs:
    def test_demo_fasta_has_five_sequences(self, tmp_path):
        fasta = make_demo_fasta(tmp_path)
        content = fasta.read_text()
        assert content.count(">") == 5

    def test_demo_busco_output_files_exist(self, tmp_path):
        make_demo_busco_outputs(tmp_path, DEMO_COMPLETENESS, DEMO_FULL_TABLE_ROWS)
        assert (tmp_path / "short_summary.txt").exists()
        assert (tmp_path / "short_summary.json").exists()
        assert (tmp_path / "full_table.tsv").exists()

    def test_demo_short_summary_parseable(self, tmp_path):
        make_demo_busco_outputs(tmp_path, DEMO_COMPLETENESS, DEMO_FULL_TABLE_ROWS)
        scores = parse_short_summary(tmp_path / "short_summary.txt")
        assert scores["C"] == pytest.approx(DEMO_COMPLETENESS["C"])
        assert scores["n"] == DEMO_COMPLETENESS["n"]


# ---------------------------------------------------------------------------
# TestRunDemo
# ---------------------------------------------------------------------------

class TestRunDemo:
    def test_run_demo_produces_all_outputs(self, tmp_path):
        args = _make_args(demo=True, output=str(tmp_path))
        run_demo(args)
        assert (tmp_path / "report.md").exists()
        assert (tmp_path / "result.json").exists()
        assert (tmp_path / "busco_run" / "short_summary.txt").exists()
        assert (tmp_path / "busco_run" / "full_table.tsv").exists()
        assert (tmp_path / "reproducibility" / "commands.sh").exists()
        assert (tmp_path / "reproducibility" / "environment.yml").exists()
        assert (tmp_path / "reproducibility" / "checksums.sha256").exists()

    def test_run_demo_report_contains_disclaimer(self, tmp_path):
        args = _make_args(demo=True, output=str(tmp_path))
        run_demo(args)
        report = (tmp_path / "report.md").read_text()
        assert DISCLAIMER in report

    def test_run_demo_result_json_has_scores(self, tmp_path):
        args = _make_args(demo=True, output=str(tmp_path))
        run_demo(args)
        data = json.loads((tmp_path / "result.json").read_text())
        assert data["C"] == pytest.approx(DEMO_COMPLETENESS["C"])
        assert data["n"] == DEMO_COMPLETENESS["n"]

    def test_run_demo_no_busco_binary_needed(self, tmp_path):
        args = _make_args(demo=True, output=str(tmp_path))
        with patch("shutil.which", return_value=None):
            run_demo(args)  # must not raise SystemExit


# ---------------------------------------------------------------------------
# TestWriteReproducibility
# ---------------------------------------------------------------------------

class TestWriteReproducibility:
    def test_commands_sh_written(self, tmp_path):
        args = _make_args(demo=True, output=str(tmp_path))
        write_reproducibility_bundle(
            tmp_path,
            cmd=["busco", "-i", "demo.fna", "-m", "genome"],
            args=args,
            input_path=None,
            is_demo=True,
        )
        commands_sh = tmp_path / "reproducibility" / "commands.sh"
        assert commands_sh.exists()
        assert "busco_assessor.py" in commands_sh.read_text()

    def test_environment_yml_has_busco(self, tmp_path):
        args = _make_args(demo=True, output=str(tmp_path))
        write_reproducibility_bundle(
            tmp_path,
            cmd=["busco", "-i", "demo.fna", "-m", "genome"],
            args=args,
            input_path=None,
            is_demo=True,
        )
        env_yml = tmp_path / "reproducibility" / "environment.yml"
        assert env_yml.exists()
        assert "busco" in env_yml.read_text()

    def test_checksums_sha256_written(self, tmp_path):
        args = _make_args(demo=True, output=str(tmp_path))
        write_reproducibility_bundle(
            tmp_path,
            cmd=["busco", "-i", "demo.fna", "-m", "genome"],
            args=args,
            input_path=None,
            is_demo=True,
        )
        assert (tmp_path / "reproducibility" / "checksums.sha256").exists()


# ---------------------------------------------------------------------------
# TestBuildBuscoCommand
# ---------------------------------------------------------------------------

class TestBuildBuscoCommand:
    def test_genome_mode_lineage_command(self, tmp_path):
        args = _make_args(input="assembly.fna", mode="genome", cpu=4)
        cmd = build_busco_command(args, "--lineage", "bacteria_odb12", tmp_path)
        assert "-m" in cmd
        assert "genome" in cmd
        assert "-l" in cmd
        assert "bacteria_odb12" in cmd

    def test_auto_lineage_prok_command(self, tmp_path):
        args = _make_args(input="assembly.fna", mode="genome", cpu=4)
        cmd = build_busco_command(args, "--auto-lineage-prok", None, tmp_path)
        assert "--auto-lineage-prok" in cmd

    def test_cpu_flag_passed(self, tmp_path):
        args = _make_args(input="assembly.fna", mode="genome", cpu=8)
        cmd = build_busco_command(args, "--lineage", "bacteria_odb12", tmp_path)
        assert "-c" in cmd
        cpu_index = cmd.index("-c")
        assert cmd[cpu_index + 1] == "8"


# ---------------------------------------------------------------------------
# NCBI taxonomy XML fixture (mirrors real efetch response for taxid 4932)
# ---------------------------------------------------------------------------

SACCHAROMYCES_XML = """\
<?xml version="1.0" ?>
<!DOCTYPE TaxaSet PUBLIC "-//NLM//DTD Taxon, 14th January 2002//EN"
    "https://www.ncbi.nlm.nih.gov/entrez/query/DTD/taxon.dtd">
<TaxaSet>
  <Taxon>
    <TaxId>4932</TaxId>
    <ScientificName>Saccharomyces cerevisiae</ScientificName>
    <LineageEx>
      <Taxon>
        <TaxId>131567</TaxId>
        <ScientificName>cellular organisms</ScientificName>
        <Rank>no rank</Rank>
      </Taxon>
      <Taxon>
        <TaxId>2759</TaxId>
        <ScientificName>Eukaryota</ScientificName>
        <Rank>superkingdom</Rank>
      </Taxon>
      <Taxon>
        <TaxId>4751</TaxId>
        <ScientificName>Fungi</ScientificName>
        <Rank>kingdom</Rank>
      </Taxon>
      <Taxon>
        <TaxId>4890</TaxId>
        <ScientificName>Ascomycota</ScientificName>
        <Rank>phylum</Rank>
      </Taxon>
      <Taxon>
        <TaxId>147537</TaxId>
        <ScientificName>Saccharomycotina</ScientificName>
        <Rank>subphylum</Rank>
      </Taxon>
      <Taxon>
        <TaxId>4891</TaxId>
        <ScientificName>Saccharomycetes</ScientificName>
        <Rank>class</Rank>
      </Taxon>
      <Taxon>
        <TaxId>4892</TaxId>
        <ScientificName>Saccharomycetales</ScientificName>
        <Rank>order</Rank>
      </Taxon>
    </LineageEx>
  </Taxon>
</TaxaSet>
"""

BACTERIA_XML = """\
<?xml version="1.0" ?>
<TaxaSet>
  <Taxon>
    <TaxId>562</TaxId>
    <ScientificName>Escherichia coli</ScientificName>
    <LineageEx>
      <Taxon>
        <TaxId>131567</TaxId>
        <ScientificName>cellular organisms</ScientificName>
        <Rank>no rank</Rank>
      </Taxon>
      <Taxon>
        <TaxId>2</TaxId>
        <ScientificName>Bacteria</ScientificName>
        <Rank>superkingdom</Rank>
      </Taxon>
      <Taxon>
        <TaxId>1236</TaxId>
        <ScientificName>Gammaproteobacteria</ScientificName>
        <Rank>class</Rank>
      </Taxon>
    </LineageEx>
  </Taxon>
</TaxaSet>
"""

SEARCH_JSON = b'{"esearchresult":{"idlist":["4932"]}}'


# ---------------------------------------------------------------------------
# TestNCBITaxonomy
# ---------------------------------------------------------------------------

class TestNCBITaxonomy:
    def _mock_urlopen(self, search_response: bytes, fetch_response: bytes):
        """Return a context manager mock that returns search then fetch responses."""
        search_cm = MagicMock()
        search_cm.__enter__ = MagicMock(return_value=MagicMock(read=MagicMock(return_value=search_response)))
        search_cm.__exit__ = MagicMock(return_value=False)

        fetch_cm = MagicMock()
        fetch_cm.__enter__ = MagicMock(return_value=MagicMock(read=MagicMock(return_value=fetch_response)))
        fetch_cm.__exit__ = MagicMock(return_value=False)

        mock = MagicMock(side_effect=[search_cm, fetch_cm])
        return mock

    def test_saccharomyces_lookup_returns_saccharomycetes(self):
        mock_urlopen = self._mock_urlopen(SEARCH_JSON, SACCHAROMYCES_XML.encode())
        with patch("urllib.request.urlopen", mock_urlopen):
            lineage = ncbi_taxonomy_lookup("Saccharomyces cerevisiae")
        names = [e["name"].lower() for e in lineage]
        assert "saccharomycetes" in names

    def test_lookup_returns_rank_name_pairs(self):
        mock_urlopen = self._mock_urlopen(SEARCH_JSON, SACCHAROMYCES_XML.encode())
        with patch("urllib.request.urlopen", mock_urlopen):
            lineage = ncbi_taxonomy_lookup("Saccharomyces cerevisiae")
        assert all("rank" in e and "name" in e for e in lineage)

    def test_network_error_returns_empty_list(self):
        import urllib.error
        with patch("urllib.request.urlopen", side_effect=urllib.error.URLError("timeout")):
            lineage = ncbi_taxonomy_lookup("Saccharomyces cerevisiae")
        assert lineage == []

    def test_ncbi_lineage_to_busco_saccharomyces(self):
        lineage = [
            {"rank": "superkingdom", "name": "Eukaryota"},
            {"rank": "kingdom",      "name": "Fungi"},
            {"rank": "phylum",       "name": "Ascomycota"},
            {"rank": "class",        "name": "Saccharomycetes"},
            {"rank": "order",        "name": "Saccharomycetales"},
        ]
        flag, value = ncbi_lineage_to_busco(lineage)
        assert flag == "--lineage"
        assert value == "saccharomycetes_odb10"

    def test_ncbi_lineage_to_busco_bacteria(self):
        lineage = [
            {"rank": "superkingdom", "name": "Bacteria"},
            {"rank": "class",        "name": "Gammaproteobacteria"},
        ]
        flag, value = ncbi_lineage_to_busco(lineage)
        assert flag == "--auto-lineage-prok"
        assert value is None

    def test_ncbi_lineage_to_busco_fungi_generic(self):
        lineage = [
            {"rank": "superkingdom", "name": "Eukaryota"},
            {"rank": "kingdom",      "name": "Fungi"},
        ]
        flag, value = ncbi_lineage_to_busco(lineage)
        assert flag == "--lineage"
        assert value == "fungi_odb10"

    def test_ncbi_lineage_to_busco_empty_falls_back_to_auto(self):
        flag, value = ncbi_lineage_to_busco([])
        assert flag == "--auto-lineage"
        assert value is None

    def test_user_agent_header_is_set(self):
        import urllib.error
        captured = []
        def fake_urlopen(req, timeout=None):
            captured.append(req)
            raise urllib.error.URLError("mock — short circuit")
        with patch("urllib.request.urlopen", side_effect=fake_urlopen):
            ncbi_taxonomy_lookup("anything")
        assert captured, "urlopen was never called"
        assert captured[0].get_header("User-agent") == NCBI_USER_AGENT

    def test_api_key_appended_when_set(self):
        import urllib.error
        captured = []
        def fake_urlopen(req, timeout=None):
            captured.append(req)
            raise urllib.error.URLError("mock — short circuit")
        with patch("urllib.request.urlopen", side_effect=fake_urlopen), \
             patch("busco_assessor.NCBI_API_KEY", "testkey123"):
            ncbi_taxonomy_lookup("Saccharomyces cerevisiae")
        assert captured, "urlopen was never called"
        assert "api_key=testkey123" in captured[0].get_full_url()

    def test_no_api_key_when_empty(self):
        import urllib.error
        captured = []
        def fake_urlopen(req, timeout=None):
            captured.append(req)
            raise urllib.error.URLError("mock — short circuit")
        with patch("urllib.request.urlopen", side_effect=fake_urlopen), \
             patch("busco_assessor.NCBI_API_KEY", ""):
            ncbi_taxonomy_lookup("Saccharomyces cerevisiae")
        assert captured, "urlopen was never called"
        assert "api_key" not in captured[0].get_full_url()

    def test_ncbi_lineage_to_busco_mammalia(self):
        lineage = [
            {"rank": "superkingdom", "name": "Eukaryota"},
            {"rank": "phylum",       "name": "Chordata"},
            {"rank": "class",        "name": "Mammalia"},
            {"rank": "order",        "name": "Primates"},
        ]
        flag, value = ncbi_lineage_to_busco(lineage)
        assert flag == "--lineage"
        assert value == "primates_odb10"


# ---------------------------------------------------------------------------
# TestDownloadDecompress
# ---------------------------------------------------------------------------

class TestDownloadDecompress:
    def _make_gz_bytes(self, content: str = ">Mito_seq\nACGTACGT\n") -> bytes:
        import io
        buf = io.BytesIO()
        with gzip.GzipFile(fileobj=buf, mode="wb") as fh:
            fh.write(content.encode())
        return buf.getvalue()

    def _mock_urlopen(self, data: bytes):
        """Return a urlopen mock that streams *data* in 65536-byte chunks."""
        class FakeResp:
            def __init__(self):
                self._data = data
                self._pos = 0
            def read(self, n=-1):
                if n == -1:
                    chunk = self._data[self._pos:]
                    self._pos = len(self._data)
                else:
                    chunk = self._data[self._pos:self._pos + n]
                    self._pos += len(chunk)
                return chunk
            def __enter__(self): return self
            def __exit__(self, *a): pass
        return MagicMock(return_value=FakeResp())

    def test_creates_decompressed_fasta(self, tmp_path):
        gz_bytes = self._make_gz_bytes()
        with patch("urllib.request.urlopen", self._mock_urlopen(gz_bytes)):
            result = download_and_decompress(
                "http://fake.example.com/test.fa.gz", tmp_path / "out"
            )
        assert result.exists()
        assert result.name.endswith(".fa")
        assert ">Mito_seq" in result.read_text()

    def test_removes_gz_after_decompress(self, tmp_path):
        gz_bytes = self._make_gz_bytes()
        with patch("urllib.request.urlopen", self._mock_urlopen(gz_bytes)):
            result = download_and_decompress(
                "http://fake.example.com/test.fa.gz", tmp_path / "out"
            )
        gz_result = result.parent / (result.name + ".gz")
        assert not gz_result.exists()

    def test_raises_on_oversized_download(self, tmp_path):
        oversized = b"X" * (MAX_DOWNLOAD_BYTES + 1)
        with patch("urllib.request.urlopen", self._mock_urlopen(oversized)), \
             pytest.raises(RuntimeError, match="exceeded"):
            download_and_decompress(
                "http://fake.example.com/huge.fa.gz", tmp_path / "out"
            )

    def test_gz_cleaned_up_on_size_cap_exceeded(self, tmp_path):
        oversized = b"X" * (MAX_DOWNLOAD_BYTES + 1)
        out_dir = tmp_path / "out"
        with patch("urllib.request.urlopen", self._mock_urlopen(oversized)), \
             pytest.raises(RuntimeError):
            download_and_decompress("http://fake.example.com/huge.fa.gz", out_dir)
        leftover = list(out_dir.glob("*.gz")) if out_dir.exists() else []
        assert leftover == []


# ---------------------------------------------------------------------------
# TestRunDemoLive
# ---------------------------------------------------------------------------

class TestRunDemoLive:
    def _setup_mock_download(self, tmp_path):
        """Patch download_and_decompress to write a tiny 5-seq FASTA instead."""
        def fake_download(url, dest_dir):
            dest_dir.mkdir(parents=True, exist_ok=True)
            fasta = dest_dir / "Saccharomyces_cerevisiae.R64-1-1.dna.chromosome.Mito.fa"
            fasta.write_text(">Mito_seq1\nACGTACGT\n>Mito_seq2\nGCATGCAT\n")
            return fasta
        return fake_download

    def test_run_demo_live_produces_all_outputs(self, tmp_path):
        args = _make_args(demo=False, output=str(tmp_path))
        fake_dl = self._setup_mock_download(tmp_path)
        with patch("busco_assessor.download_and_decompress", side_effect=fake_dl), \
             patch("busco_assessor.check_busco", return_value=False), \
             patch("busco_assessor.ncbi_taxonomy_lookup", return_value=[
                 {"rank": "kingdom", "name": "Fungi"},
                 {"rank": "class",   "name": "Saccharomycetes"},
             ]):
            run_demo_live(args)

        assert (tmp_path / "report.md").exists()
        assert (tmp_path / "result.json").exists()
        assert (tmp_path / "busco_run" / "short_summary.txt").exists()
        assert (tmp_path / "busco_run" / "full_table.tsv").exists()
        assert (tmp_path / "reproducibility" / "commands.sh").exists()
        assert (tmp_path / "reproducibility" / "environment.yml").exists()

    def test_run_demo_live_report_contains_mito_note(self, tmp_path):
        args = _make_args(demo=False, output=str(tmp_path))
        fake_dl = self._setup_mock_download(tmp_path)
        with patch("busco_assessor.download_and_decompress", side_effect=fake_dl), \
             patch("busco_assessor.check_busco", return_value=False), \
             patch("busco_assessor.ncbi_taxonomy_lookup", return_value=[]):
            run_demo_live(args)
        report = (tmp_path / "report.md").read_text()
        assert "mitochondri" in report.lower()

    def test_run_demo_live_low_completeness_scores(self, tmp_path):
        args = _make_args(demo=False, output=str(tmp_path))
        fake_dl = self._setup_mock_download(tmp_path)
        with patch("busco_assessor.download_and_decompress", side_effect=fake_dl), \
             patch("busco_assessor.check_busco", return_value=False), \
             patch("busco_assessor.ncbi_taxonomy_lookup", return_value=[]):
            run_demo_live(args)
        data = json.loads((tmp_path / "result.json").read_text())
        assert data["C"] < 10.0  # mito genome gives low completeness

    def test_run_demo_live_disclaimer_present(self, tmp_path):
        args = _make_args(demo=False, output=str(tmp_path))
        fake_dl = self._setup_mock_download(tmp_path)
        with patch("busco_assessor.download_and_decompress", side_effect=fake_dl), \
             patch("busco_assessor.check_busco", return_value=False), \
             patch("busco_assessor.ncbi_taxonomy_lookup", return_value=[]):
            run_demo_live(args)
        report = (tmp_path / "report.md").read_text()
        assert DISCLAIMER in report

    def test_run_demo_live_result_json_has_lineage(self, tmp_path):
        args = _make_args(demo=False, output=str(tmp_path))
        fake_dl = self._setup_mock_download(tmp_path)
        with patch("busco_assessor.download_and_decompress", side_effect=fake_dl), \
             patch("busco_assessor.check_busco", return_value=False), \
             patch("busco_assessor.ncbi_taxonomy_lookup", return_value=[
                 {"rank": "class", "name": "Saccharomycetes"}
             ]):
            run_demo_live(args)
        data = json.loads((tmp_path / "result.json").read_text())
        assert "lineage" in data
