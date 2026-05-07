"""Unit tests for the eQTL Catalogue region fetcher (tabix-on-FTP).

Mocks pysam.TabixFile + the REST metadata endpoint so tests run offline. Live
smoke test (test_live.py) hits the real FTP server when RUN_LIVE_TESTS=1.
"""

from __future__ import annotations

import sys
from pathlib import Path
from unittest.mock import MagicMock

import pytest

# Inject the skill dir so the bare-name `eqtl_catalogue_region_fetch`
# resolves under ClawBio's kebab-cased skill directory (`skills/eqtl-
# catalogue-region-fetch/`), which is not a valid Python identifier.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from eqtl_catalogue_region_fetch import (  # noqa: E402
    EQTLCatalogueAPIError,
    EQTLCatalogueClient,
    FTP_COLUMNS,
    RegionVariant,
    ftp_url_for,
)


# ----------------------- URL construction -----------------------


def test_ftp_url_for_canonical():
    url = ftp_url_for("QTS000015", "QTD000266")
    assert url == (
        "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/"
        "QTS000015/QTD000266/QTD000266.all.tsv.gz"
    )


def test_ftp_url_for_strips_trailing_slash_on_base():
    url = ftp_url_for("QTS000015", "QTD000266",
                      ftp_base="https://example.org/sumstats/")
    assert url == "https://example.org/sumstats/QTS000015/QTD000266/QTD000266.all.tsv.gz"


# ----------------------- FTP row parsing helpers -----------------------


SORT1_ROW_FIELDS = [
    "ENSG00000134243",  # molecular_trait_id (SORT1, canonical 1p13.3 LDL/CHD locus)
    "1",                # chromosome
    "109274968",        # position (SORT1 lead variant, Musunuru 2010)
    "T",                # ref
    "G",                # alt
    "chr1_109274968_T_G",# variant
    "5",                # ma_samples
    "0.32",             # maf
    "2.5e-15",          # pvalue
    "0.444622",         # beta
    "0.0477403",        # se
    "SNP",              # type
    "127",              # ac
    "396",              # an
    "0.999",            # r2
    "ENSG00000134243",  # molecular_trait_object_id
    "ENSG00000134243",  # gene_id
    "27.42",            # median_tpm
    "rs12740374",       # rsid
]


def _row_str(*overrides):
    """Build a tab-joined FTP row string with optional field overrides
    (positional, in FTP_COLUMNS order)."""
    fields = list(SORT1_ROW_FIELDS)
    for i, val in enumerate(overrides):
        if val is not None:
            fields[i] = str(val)
    return "\t".join(fields)


def _mock_tabix(rows: list[str], *, raise_on_unprefixed: bool = False):
    tbx = MagicMock()
    tbx.header = []  # FTP files do not start with #-prefixed header lines
    tbx.contigs = ["1", "2", "3"]
    if raise_on_unprefixed:
        # Force a chr-prefix retry path.
        def fetch(chrom, *args, **kwargs):
            if not chrom.startswith("chr"):
                raise ValueError("contig not found")
            return iter(rows)
        tbx.fetch = fetch
    else:
        tbx.fetch = MagicMock(return_value=iter(rows))
    tbx.close = MagicMock()
    return tbx


@pytest.fixture
def patched_pysam(monkeypatch):
    """Patch pysam.TabixFile for the duration of one test, returning a holder
    where the test can install a MagicMock per call.
    """
    holder = {"tbx": None}

    def factory(url):
        return holder["tbx"]

    import pysam
    monkeypatch.setattr(pysam, "TabixFile", factory)
    return holder


@pytest.fixture
def client_with_metadata(monkeypatch):
    """EQTLCatalogueClient whose REST `_get` is mocked to return a stable
    GTEx liver ge-eQTL (QTD000266) metadata dict."""
    client = EQTLCatalogueClient()
    monkeypatch.setattr(client, "_get", lambda path, params=None: {
        "study_id": "QTS000015",
        "dataset_id": "QTD000266",
        "study_label": "GTEx",
        "tissue_label": "liver",
        "condition_label": "naive",
        "sample_group": "liver",
        "quant_method": "ge",
        "release": "7.0",
    })
    return client


# ----------------------- fetch_region happy path -----------------------


def test_fetch_region_returns_harmonised_variants(client_with_metadata, patched_pysam):
    rows = [
        _row_str(),                           # SORT1 row
        _row_str("ENSG00000134243", "2", "109280000", "A", "G",
                 "chr1_109280000_A_G", None, None, "1.2e-15",
                 "0.445", "0.046", "SNP"),
        _row_str("ENSG00000999999", "2", "109265000", "A", "G",
                 "chr1_109265000_A_G", None, None, "0.5",
                 "0.001", "0.05", "SNP"),
    ]
    patched_pysam["tbx"] = _mock_tabix(rows)
    result = client_with_metadata.fetch_region(
        dataset_id="QTD000266", chromosome="1",
        start_bp=108_774_968, end_bp=109_774_968,
        molecular_trait_id="ENSG00000134243",
    )
    # Filter by molecular_trait_id should drop the third row.
    assert result.n_variants == 2
    v0 = result.variants[0]
    assert isinstance(v0, RegionVariant)
    assert v0.variant_id == "1_109274968_T_G"  # chr prefix stripped
    assert v0.chromosome == "1"
    assert v0.position == 109_274_968
    assert v0.ref == "T" and v0.alt == "G"
    assert v0.beta == pytest.approx(0.444622)
    assert v0.se == pytest.approx(0.0477403)
    assert v0.p_value == pytest.approx(2.5e-15)
    assert v0.maf == pytest.approx(0.32)
    patched_pysam["tbx"].close.assert_called_once()


def test_fetch_region_records_release(client_with_metadata, patched_pysam):
    patched_pysam["tbx"] = _mock_tabix([_row_str()])
    result = client_with_metadata.fetch_region(
        dataset_id="QTD000266", chromosome="1",
        start_bp=108_774_968, end_bp=109_774_968,
        molecular_trait_id="ENSG00000134243",
    )
    assert result.release.dataset_release == "7.0"
    assert result.release.fetched_at_utc.endswith("Z")
    assert result.release.study_label == "GTEx"
    assert result.release.sample_group == "liver"


def test_fetch_region_no_molecular_trait_filter_returns_all(client_with_metadata, patched_pysam):
    rows = [_row_str(), _row_str("ENSG00000999999", "2", "109265000", "A", "G",
                                  "chr1_109265000_A_G")]
    patched_pysam["tbx"] = _mock_tabix(rows)
    result = client_with_metadata.fetch_region(
        dataset_id="QTD000266", chromosome="1",
        start_bp=108_774_968, end_bp=109_774_968,
    )
    assert result.n_variants == 2


def test_fetch_region_explicit_study_id_skips_metadata_call(monkeypatch, patched_pysam):
    """When study_id is passed explicitly the client should not make a REST call."""
    client = EQTLCatalogueClient()
    calls = {"n": 0}

    def fake_get(path, params=None):
        calls["n"] += 1
        return {"study_id": "QTS000015", "release": "7.0"}

    monkeypatch.setattr(client, "_get", fake_get)
    patched_pysam["tbx"] = _mock_tabix([_row_str()])
    result = client.fetch_region(
        dataset_id="QTD000266", chromosome="1",
        start_bp=108_774_968, end_bp=109_774_968,
        study_id="QTS000015",
    )
    # Expect 1 metadata call (for release info), not 2.
    assert calls["n"] == 1
    assert result.n_variants == 1


def test_fetch_region_raises_on_unopenable_index(monkeypatch):
    client = EQTLCatalogueClient()
    monkeypatch.setattr(client, "_get", lambda *a, **kw: {"study_id": "QTS000015"})

    def boom(url):
        raise OSError("tabix open failed")
    import pysam
    monkeypatch.setattr(pysam, "TabixFile", boom)
    with pytest.raises(EQTLCatalogueAPIError, match="could not open tabix"):
        client.fetch_region(
            dataset_id="QTD000266", chromosome="1",
            start_bp=1, end_bp=1000,
        )


def test_fetch_region_skips_rows_with_wrong_column_count(client_with_metadata, patched_pysam):
    rows = [
        _row_str(),                                    # 19 cols, OK
        "\t".join(["ENSG", "2", "109265000", "C", "T"]) # 5 cols, malformed
    ]
    patched_pysam["tbx"] = _mock_tabix(rows)
    result = client_with_metadata.fetch_region(
        dataset_id="QTD000266", chromosome="1",
        start_bp=108_774_968, end_bp=109_774_968,
    )
    assert result.n_variants == 1
    assert any("schema may have drifted" in n for n in result.notes)


def test_fetch_region_chr_prefix_retry(monkeypatch, client_with_metadata, patched_pysam):
    """Tabix indexes can use 'chr' prefix; client retries with the alt form."""
    patched_pysam["tbx"] = _mock_tabix([_row_str()], raise_on_unprefixed=True)
    result = client_with_metadata.fetch_region(
        dataset_id="QTD000266", chromosome="1",
        start_bp=108_774_968, end_bp=109_774_968,
    )
    assert result.n_variants == 1


def test_resolve_study_id_falls_back_to_metadata(monkeypatch, patched_pysam):
    client = EQTLCatalogueClient()
    monkeypatch.setattr(client, "_get",
                        lambda *a, **kw: {"study_id": "QTS000015", "release": "7.0"})
    patched_pysam["tbx"] = _mock_tabix([_row_str()])
    result = client.fetch_region(
        dataset_id="QTD000266", chromosome="1",
        start_bp=108_774_968, end_bp=109_774_968,
    )
    assert result.n_variants == 1


def test_resolve_study_id_raises_when_metadata_lacks_study_id(monkeypatch):
    client = EQTLCatalogueClient()
    monkeypatch.setattr(client, "_get", lambda *a, **kw: {"release": "7.0"})
    with pytest.raises(EQTLCatalogueAPIError, match="could not resolve study_id"):
        client.fetch_region(
            dataset_id="QTD000266", chromosome="1",
            start_bp=1, end_bp=1000,
        )


# ----------------------- FTP_COLUMNS schema check -----------------------


def test_ftp_columns_count_is_19():
    # If eQTL Catalogue ever changes the schema this test fails loudly.
    assert len(FTP_COLUMNS) == 19
    assert FTP_COLUMNS[0] == "molecular_trait_id"
    assert FTP_COLUMNS[5] == "variant"
    assert FTP_COLUMNS[8] == "pvalue"
    assert FTP_COLUMNS[9] == "beta"
