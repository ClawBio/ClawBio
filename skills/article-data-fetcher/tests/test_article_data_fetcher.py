"""Tests for article-data-fetcher skill."""

import json
import re
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
import sys

sys.path.insert(0, str(Path(__file__).parent.parent))
from article_data_fetcher import (
    ACCESSION_PATTERNS,
    discover_accessions,
    file_matches_types,
    parse_identifier,
)


# ---------------------------------------------------------------------------
# parse_identifier
# ---------------------------------------------------------------------------

class TestParseIdentifier:
    def test_doi(self):
        kind, value = parse_identifier("10.1038/s41586-021-03819-2")
        assert kind == "doi"
        assert value == "10.1038/s41586-021-03819-2"

    def test_pmid_bare(self):
        kind, value = parse_identifier("34613072")
        assert kind == "pmid"
        assert value == "34613072"

    def test_pmid_prefixed(self):
        kind, value = parse_identifier("PMID:34613072")
        assert kind == "pmid"
        assert value == "34613072"

    def test_url(self):
        kind, value = parse_identifier("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123456")
        assert kind == "url"

    def test_invalid(self):
        with pytest.raises(ValueError):
            parse_identifier("not-an-id")


# ---------------------------------------------------------------------------
# discover_accessions
# ---------------------------------------------------------------------------

class TestDiscoverAccessions:
    def test_geo(self):
        text = "Data available at GSE145926 and GSE200001."
        found = discover_accessions(text)
        assert "geo" in found
        assert "GSE145926" in found["geo"]
        assert "GSE200001" in found["geo"]

    def test_zenodo(self):
        text = "Deposited at 10.5281/zenodo.7654321"
        found = discover_accessions(text)
        assert "zenodo" in found

    def test_sra(self):
        text = "Raw reads: PRJNA654321"
        found = discover_accessions(text)
        assert "sra" in found

    def test_no_accessions(self):
        text = "This paper has no data."
        found = discover_accessions(text)
        assert found == {}

    def test_multiple_repos(self):
        text = "GEO: GSE185224. Zenodo: 10.5281/zenodo.1234567."
        found = discover_accessions(text)
        assert "geo" in found
        assert "zenodo" in found


# ---------------------------------------------------------------------------
# file_matches_types
# ---------------------------------------------------------------------------

class TestFileMatchesTypes:
    def test_vcf_match(self):
        assert file_matches_types("variants.vcf", {"vcf"})

    def test_vcf_gz_match(self):
        assert file_matches_types("variants.vcf.gz", {"vcf.gz"})

    def test_all_matches_everything(self):
        assert file_matches_types("anything.bam", {"all"})

    def test_no_match(self):
        assert not file_matches_types("matrix.h5ad", {"vcf", "fasta"})

    def test_csv_match(self):
        assert file_matches_types("metadata.csv", {"csv", "vcf"})

    def test_case_insensitive(self):
        assert file_matches_types("Sample.VCF", {"vcf"})

    def test_h5ad_match(self):
        assert file_matches_types("counts.h5ad", {"h5ad"})


# ---------------------------------------------------------------------------
# _parse_size (Apache listing sizes)
# ---------------------------------------------------------------------------

from hypothesis import example, given  # noqa: E402
from hypothesis import strategies as st  # noqa: E402

import article_data_fetcher as adf  # noqa: E402

_UNITS = {"": 1, "K": 1024, "M": 1024**2, "G": 1024**3, "T": 1024**4}


@given(n=st.integers(0, 10_000), unit=st.sampled_from(sorted(_UNITS)))
@example(n=945, unit="K")  # sizes as GEO's listing prints them
@example(n=14, unit="M")
def test_parse_size_reads_whole_numbers_with_units(n: int, unit: str) -> None:
    assert adf._parse_size(f"{n}{unit}") == n * _UNITS[unit]
    assert adf._parse_size(f"{n}{unit.lower()}") == n * _UNITS[unit]


@given(a=st.integers(0, 9_999), unit=st.sampled_from(["K", "M", "G"]))
def test_parse_size_is_monotonic_in_tenths(a: int, unit: str) -> None:
    """'1.4M' < '1.5M': a larger listed size never parses to fewer bytes."""
    lo, hi = f"{a // 10}.{a % 10}{unit}", f"{(a + 1) // 10}.{(a + 1) % 10}{unit}"
    assert adf._parse_size(lo) <= adf._parse_size(hi)


@pytest.mark.parametrize("text", ["", "-", "  ", "K", "1.2.3M", "abc"])
def test_parse_size_unknown_is_none(text: str) -> None:
    assert adf._parse_size(text) is None


# ---------------------------------------------------------------------------
# Demo: offline, non-interactive, only under --output
# ---------------------------------------------------------------------------

import os  # noqa: E402
import shutil  # noqa: E402
import subprocess  # noqa: E402

from clawbio.common.checksums import sha256_file  # noqa: E402

SCRIPT = Path(__file__).resolve().parents[1] / "article_data_fetcher.py"
DEMO_FILE = "GSE100866_CD8_merged-ADT_umi.csv"


@pytest.fixture(autouse=True)
def _cwd_in_tmp(tmp_path: Path, monkeypatch) -> None:
    """Run each test, and any subprocess it starts, from its own tmp_path (AGENTS.md)."""
    monkeypatch.chdir(tmp_path)


def _offline_env() -> dict[str, str]:
    """The current environment with every proxy pointed at a closed port."""
    env = {k: v for k, v in os.environ.items() if "proxy" not in k.lower()}
    dead = "http://127.0.0.1:9"
    env.update(HTTP_PROXY=dead, HTTPS_PROXY=dead, http_proxy=dead, https_proxy=dead, NO_PROXY="", no_proxy="")
    return env


@pytest.fixture(scope="module")
def demo_run(tmp_path_factory) -> tuple[Path, Path]:
    cwd = tmp_path_factory.mktemp("cwd")
    out = cwd / "out"
    subprocess.run([sys.executable, str(SCRIPT), "--demo", "--output", str(out)],
                   check=True, capture_output=True, text=True, cwd=cwd,
                   stdin=subprocess.DEVNULL, env=_offline_env(), timeout=120)
    return cwd, out


def test_demo_runs_offline_without_prompting(demo_run) -> None:
    cwd, out = demo_run
    assert (out / "GSE100866" / DEMO_FILE).is_file()
    assert sorted(p.name for p in cwd.iterdir()) == ["out"], "demo wrote outside --output"
    report = (out / "report.md").read_text()
    assert "GSE100866" in report and DEMO_FILE in report
    assert adf.DISCLAIMER in report


def test_demo_manifest_records_no_absolute_paths(demo_run) -> None:
    _, out = demo_run
    manifest = json.loads((out / "manifest.json").read_text())
    assert [f["local_path"] for f in manifest["files"] if f["downloaded"]] == [f"GSE100866/{DEMO_FILE}"]
    assert str(out) not in (out / "manifest.json").read_text()


def test_demo_bundle_verifies_and_replays(demo_run, tmp_path: Path) -> None:
    _, out = demo_run
    lines = (out / "reproducibility" / "checksums.sha256").read_text().splitlines()
    labels = {label for _, label in (line.split("  ", 1) for line in lines)}
    assert {"report.md", "manifest.json", "result.json", f"GSE100866/{DEMO_FILE}",
            "reproducibility/commands.sh", "reproducibility/environment.yml"} <= labels
    for digest, label in (line.split("  ", 1) for line in lines):
        assert digest == sha256_file(out / label), label

    copy, elsewhere = tmp_path / "copy", tmp_path / "elsewhere"
    shutil.copytree(out, copy)
    elsewhere.mkdir()
    env = _offline_env() | {"PYTHON": sys.executable}
    subprocess.run(["bash", str(copy / "reproducibility" / "commands.sh")], check=True,
                   capture_output=True, cwd=elsewhere, stdin=subprocess.DEVNULL, env=env, timeout=120)
    # The replay rewrites the bundle; its checksums verify against its own outputs
    # (result.json carries a run timestamp), and the downloaded data is identical.
    for line in (copy / "reproducibility" / "checksums.sha256").read_text().splitlines():
        digest, label = line.split("  ", 1)
        assert digest == sha256_file(copy / label), label
    data = f"GSE100866/{DEMO_FILE}"
    assert sha256_file(copy / data) == sha256_file(out / data)
    assert not any(elsewhere.iterdir()), "replay wrote into the working directory"


# ---------------------------------------------------------------------------
# Demo output snapshot
# ---------------------------------------------------------------------------
# Fails when a code change alters any demo number. If the change is intended,
# run `uv run pytest skills/article-data-fetcher/tests --inline-snapshot=fix`; it
# rewrites the snapshot() values below, and the PR diff shows which numbers moved.

from inline_snapshot import snapshot  # noqa: E402

from clawbio.common import demo_runner  # noqa: E402


@pytest.mark.demo_snapshot
def test_demo_output_matches_snapshot(tmp_path):
    assert demo_runner.demo_outputs("article-data-fetcher", tmp_path) == snapshot({'manifest.json': """\
{
  "demo": true,
  "accession": "GSE100866",
  "files": [
    {
      "filename": "GSE100866_CD8_merged-ADT_umi.csv.gz",
      "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE100nnn/GSE100866/suppl/GSE100866_CD8_merged-ADT_umi.csv.gz",
      "repository": "GEO",
      "accession": "GSE100866",
      "size_bytes": 19456,
      "data_type": "count matrix",
      "local_path": "GSE100866/GSE100866_CD8_merged-ADT_umi.csv",
      "downloaded": true
    }
  ]
}
""", 'report.md': """\
# Article Data Fetcher — Demo Report

**Accession**: GSE100866 (GEO, recorded responses)
**Files listed**: 9
**Files downloaded**: 1

## Files Downloaded

- `GSE100866/GSE100866_CD8_merged-ADT_umi.csv` (19.0 KB, count matrix)

---
*ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions.*
""", 'result.json': """\
{
  "skill": "article-data-fetcher",
  "version": "0.1.0",
  "completed_at": "<TIMESTAMP>",
  "input_checksum": "",
  "datasets": {},
  "summary": {
    "accession": "GSE100866",
    "files_listed": 9,
    "files_downloaded": 1
  },
  "data": {
    "files": [
      {
        "filename": "GSE100866_CD8_merged-ADT_umi.csv.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE100nnn/GSE100866/suppl/GSE100866_CD8_merged-ADT_umi.csv.gz",
        "repository": "GEO",
        "accession": "GSE100866",
        "size_bytes": 19456,
        "data_type": "count matrix",
        "local_path": "GSE100866/GSE100866_CD8_merged-ADT_umi.csv",
        "downloaded": true
      }
    ]
  }
}\
"""})
