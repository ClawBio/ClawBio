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
