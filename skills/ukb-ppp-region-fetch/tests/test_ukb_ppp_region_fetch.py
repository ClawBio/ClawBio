"""Unit tests for the UKB-PPP region fetcher (Synapse-backed).

Mocks the Synapse REST listing endpoint + the per-protein tar archive so
tests run offline. Live smoke test (test_live_ukb_ppp_region_fetch.py)
hits the real Synapse backend when RUN_LIVE_TESTS=1 + SYNAPSE_AUTH_TOKEN
is set.
"""

from __future__ import annotations

import gzip
import io
import json
import sys
import tarfile
from pathlib import Path
from unittest.mock import MagicMock

import pytest
from hypothesis import given, strategies as st

# Inject the skill dir so the bare-name `ukb_ppp_region_fetch` resolves
# under ClawBio's kebab-cased skill directory (`skills/ukb-ppp-region-fetch/`),
# which is not a valid Python identifier.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from ukb_ppp_region_fetch import (  # noqa: E402
    ANCESTRY_FOLDER_IDS,
    ANCESTRY_LABELS,
    PAT_REQUIRED_MESSAGE,
    PROTEIN_FILE_RE,
    REGENIE_COLUMNS,
    RegionResult,
    RegionVariant,
    SynapseListFetcher,
    UKBPPPAccessError,
    UKBPPPClient,
    UKBPPPRelease,
    UKB_PPP_COHORT_N,
    UKB_PPP_COHORT_N_SOURCE,
    save_region_result_as_bundled_slice,
)
from ukb_ppp_region_fetch import (  # noqa: E402
    _bundled_slice_key,
    _bundled_slices_dir,
    _extract_region_from_tar,
    _load_bundled_slice,
    _maybe_int,
    _n_samples_from_rows,
    _normalise_row,
    _ProteinIndex,
    main,
)


# ----------------------- File-name regex -----------------------


def test_protein_file_re_matches_canonical_sort1():
    """SORT1's actual UKB-PPP file name (verified against Synapse listing).

    SORT1 is Sortilin / NTR3; UniProt Q99523; sits on the Cardiometabolic
    Olink panel.
    """
    name = "SORT1_Q99523_OID20684_v1_Cardiometabolic.tar"
    m = PROTEIN_FILE_RE.match(name)
    assert m is not None
    assert m.group("hgnc") == "SORT1"
    assert m.group("uniprot") == "Q99523"
    assert m.group("olink_id") == "OID20684"
    assert m.group("panel") == "Cardiometabolic"


def test_protein_file_re_matches_complex_panel_names():
    """Some panels carry roman-numeral suffixes (e.g. Inflammation_II)."""
    name = "A1BG_P04217_OID30771_v1_Inflammation_II.tar"
    m = PROTEIN_FILE_RE.match(name)
    assert m is not None
    assert m.group("panel") == "Inflammation_II"


def test_protein_file_re_rejects_non_protein_files():
    """README / metadata files should not be parsed as protein archives."""
    for bad in ("README.md", "metadata.tsv", "snp_map_chr1.tsv.gz",
                "SORT1_no_uniprot.tar"):
        assert PROTEIN_FILE_RE.match(bad) is None


# ----------------------- Static enums -----------------------


def test_ancestry_folder_ids_cover_all_published_strata():
    """Sun 2023 Nature reports 7 ancestry strata (incl. Combined)."""
    expected = {"EUR", "AFR", "EAS", "SAS", "CSA", "MID", "AMR", "ALL"}
    assert set(ANCESTRY_FOLDER_IDS) == expected


def test_ancestry_labels_user_friendly_for_all_codes():
    for code in ANCESTRY_FOLDER_IDS:
        assert code in ANCESTRY_LABELS
        assert len(ANCESTRY_LABELS[code]) > len(code)  # not just echoing the code


def test_regenie_columns_count_is_14():
    """Per Sun 2023 methods + REGENIE step-2 output. If the schema ever
    drifts upstream the parser's `len(fields) != len(header)` check
    flags drift loudly, not silently.
    """
    assert len(REGENIE_COLUMNS) == 14


def test_cohort_n_matches_sun_2023():
    """Published cohort sizes, quoted from Sun BB et al. Nature 2023;622:329-338
    (doi:10.1038/s41586-023-06592-6, PMID 37794186), Results:

    - "Discovery pQTL analyses were performed in participants of European
      ancestry from the randomly selected baseline cohort ( n = 34,557)"
    - "ancestry-specific pQTL analyses within five ancestry groups as defined
      by pan-UKBB: African (AFR, n = 931), Central/South Asian (CSA, n = 920),
      Middle Eastern (MID, n = 308), East Asian (EAS, n = 262) and admixed
      American (AMR, n = 97)"
    - "we also performed pQTL mapping for the full UKB-PPP cohort ( n = 52,363)"

    SAS is an alias of CSA. These are cohort sizes, not the N of any slice.
    """
    assert UKB_PPP_COHORT_N == {
        "EUR": 34_557, "AFR": 931, "CSA": 920, "SAS": 920,
        "EAS": 262, "MID": 308, "AMR": 97, "ALL": 52_363,
    }
    assert "37794186" in UKB_PPP_COHORT_N_SOURCE
    assert "10.1038/s41586-023-06592-6" in UKB_PPP_COHORT_N_SOURCE


def test_ancestry_label_all_is_not_a_meta_analysis():
    """The full-cohort (ALL) analysis is one GWAS of all participants
    (n = 52,363), not a meta-analysis of the strata."""
    assert "meta" not in ANCESTRY_LABELS["ALL"].lower()



# ----------------------- Protein index resolution -----------------------


def _stub_listing(name_uniprot_panel_id: list[tuple[str, str, str, str]]):
    """Build a Synapse listing-shape payload from (hgnc, uniprot, olink_id, panel)."""
    return [
        {
            "name": f"{hgnc}_{uni}_{oid}_v1_{panel}.tar",
            "id": f"syn{i + 60000000}",
        }
        for i, (hgnc, uni, oid, panel) in enumerate(name_uniprot_panel_id)
    ]


@pytest.fixture
def stub_lister():
    lister = SynapseListFetcher()
    lister.list_folder = MagicMock()  # type: ignore[assignment]
    return lister


def test_protein_index_resolves_by_hgnc(stub_lister):
    stub_lister.list_folder.return_value = _stub_listing([
        ("SORT1", "Q99523", "OID20684", "Cardiometabolic"),
        ("IL6R", "P08887", "OID20677", "Cardiometabolic"),
    ])
    idx = _ProteinIndex(stub_lister, "EUR")
    rec = idx.resolve("SORT1")
    assert rec["hgnc"] == "SORT1"
    assert rec["uniprot"] == "Q99523"
    assert rec["olink_id"] == "OID20684"


def test_protein_index_resolves_by_uniprot(stub_lister):
    stub_lister.list_folder.return_value = _stub_listing([
        ("SORT1", "Q99523", "OID20684", "Cardiometabolic"),
    ])
    idx = _ProteinIndex(stub_lister, "EUR")
    rec = idx.resolve("Q99523")
    assert rec["hgnc"] == "SORT1"


def test_protein_index_raises_for_unknown(stub_lister):
    stub_lister.list_folder.return_value = _stub_listing([
        ("SORT1", "Q99523", "OID20684", "Cardiometabolic"),
    ])
    idx = _ProteinIndex(stub_lister, "EUR")
    with pytest.raises(UKBPPPAccessError, match="no UKB-PPP protein file"):
        idx.resolve("NOTAPROTEIN")


def test_protein_index_lazy_loads_once(stub_lister):
    """Repeated lookups in the same ancestry hit Synapse once."""
    stub_lister.list_folder.return_value = _stub_listing([
        ("SORT1", "Q99523", "OID20684", "Cardiometabolic"),
        ("IL6R", "P08887", "OID20677", "Cardiometabolic"),
    ])
    idx = _ProteinIndex(stub_lister, "EUR")
    idx.resolve("SORT1")
    idx.resolve("IL6R")
    assert stub_lister.list_folder.call_count == 1


def test_protein_index_raises_on_unknown_ancestry(stub_lister):
    idx = _ProteinIndex(stub_lister, "ZZZ")
    with pytest.raises(UKBPPPAccessError, match="unknown UKB-PPP ancestry"):
        idx.resolve("SORT1")


def test_protein_index_skips_non_protein_entries(stub_lister):
    """README + metadata files in the folder don't break index build."""
    stub_lister.list_folder.return_value = [
        {"name": "SORT1_Q99523_OID20684_v1_Cardiometabolic.tar", "id": "syn1"},
        {"name": "README.md", "id": "syn2"},
        {"name": "ancestry_metadata.tsv", "id": "syn3"},
    ]
    idx = _ProteinIndex(stub_lister, "EUR")
    assert idx.resolve("SORT1")["synapse_id"] == "syn1"


# ----------------------- REGENIE row normalisation -----------------------


def test_normalise_row_harmonises_to_alt_effect():
    """OT convention: chr_pos_ref_alt, ALT = effect. REGENIE matches:
    ALLELE0 = REF, ALLELE1 = ALT (effect)."""
    row = {
        "CHROM": "1", "GENPOS": "109274968",
        "ALLELE0": "T", "ALLELE1": "G",
        "A1FREQ": "0.32", "INFO": "0.999", "N": "33628", "TEST": "ADD",
        "BETA": "0.44", "SE": "0.05", "CHISQ": "75.8",
        "LOG10P": "12.5", "EXTRA": "",
    }
    v = _normalise_row(row)
    assert v.variant_id == "1_109274968_T_G"
    assert v.chromosome == "1"
    assert v.position == 109274968
    assert v.ref == "T"
    assert v.alt == "G"
    assert v.beta == 0.44
    assert v.se == 0.05
    assert v.effect_allele_frequency == 0.32
    assert v.maf == 0.32  # already <= 0.5
    # 10^-12.5 ~= 3.16e-13
    assert v.p_value is not None
    assert abs(v.p_value - 10 ** -12.5) < 1e-25


def test_normalise_row_maf_fold_when_eaf_above_half():
    """A1FREQ > 0.5 should fold to MAF = 1 - A1FREQ."""
    row = {
        "CHROM": "1", "GENPOS": "109274968",
        "ALLELE0": "T", "ALLELE1": "G",
        "A1FREQ": "0.78", "INFO": "0.999", "N": "33628", "TEST": "ADD",
        "BETA": "0.44", "SE": "0.05", "CHISQ": "75.8",
        "LOG10P": "12.5", "EXTRA": "",
    }
    v = _normalise_row(row)
    assert v.effect_allele_frequency == 0.78
    assert v.maf == pytest.approx(0.22)


def test_normalise_row_handles_missing_logp():
    row = {
        "CHROM": "1", "GENPOS": "109274968",
        "ALLELE0": "T", "ALLELE1": "G",
        "A1FREQ": "NA", "INFO": "0.5", "N": "33628", "TEST": "ADD",
        "BETA": "NA", "SE": "NA", "CHISQ": "NA",
        "LOG10P": "NA", "EXTRA": "",
    }
    v = _normalise_row(row)
    assert v.beta is None
    assert v.se is None
    assert v.p_value is None
    assert v.effect_allele_frequency is None
    assert v.maf is None


def test_normalise_row_clamps_overflow_to_zero():
    """Very small p-values (LOG10P > ~300) overflow to 0.0 rather than err."""
    row = {
        "CHROM": "1", "GENPOS": "109274968",
        "ALLELE0": "T", "ALLELE1": "G",
        "A1FREQ": "0.32", "INFO": "0.999", "N": "33628", "TEST": "ADD",
        "BETA": "0.44", "SE": "0.05", "CHISQ": "9999",
        "LOG10P": "500", "EXTRA": "",
    }
    v = _normalise_row(row)
    # 10^-500 underflows Python float; we treat as 0.0 not error.
    assert v.p_value == 0.0


# ----------------------- Tar extraction + region filter -----------------------


def _build_test_tar(tmp_path: Path, *, chrom: str, rows: list[dict[str, str]]) -> Path:
    """Build an in-memory tar with one chr<N> gzipped REGENIE file."""
    out_path = tmp_path / "SORT1_Q99523_OID20684_v1_Cardiometabolic.tar"
    header_line = " ".join(REGENIE_COLUMNS)
    data_lines = [
        " ".join(r[c] for c in REGENIE_COLUMNS) for r in rows
    ]
    payload = "\n".join([header_line] + data_lines) + "\n"
    gz_blob = gzip.compress(payload.encode("utf-8"))
    inner_name = f"discovery_chr{chrom}_SORT1_Q99523_OID20684_v1_Cardiometabolic.regenie.gz"
    with tarfile.open(out_path, mode="w") as tf:
        info = tarfile.TarInfo(name=inner_name)
        info.size = len(gz_blob)
        tf.addfile(info, io.BytesIO(gz_blob))
    return out_path


@pytest.fixture
def synthetic_tar(tmp_path: Path) -> Path:
    # REGENIE step-2 fills EXTRA with a non-empty token (e.g. "0" or
    # "NA"); the parser splits on whitespace so an empty trailing field
    # would collapse and break column count. Use "0" here to mirror
    # actual REGENIE output.
    rows = [
        # In-window lead variant
        {
            "CHROM": "1", "GENPOS": "109274968", "ID": "rs12740374",
            "ALLELE0": "T", "ALLELE1": "G",
            "A1FREQ": "0.32", "INFO": "0.999", "N": "33628", "TEST": "ADD",
            "BETA": "0.44", "SE": "0.05", "CHISQ": "75.8",
            "LOG10P": "12.5", "EXTRA": "0",
        },
        # In-window, weak signal
        {
            "CHROM": "1", "GENPOS": "109274970", "ID": "rsX",
            "ALLELE0": "A", "ALLELE1": "C",
            "A1FREQ": "0.05", "INFO": "0.95", "N": "33628", "TEST": "ADD",
            "BETA": "0.02", "SE": "0.04", "CHISQ": "0.25",
            "LOG10P": "0.4", "EXTRA": "0",
        },
        # Out-of-window (below)
        {
            "CHROM": "1", "GENPOS": "108774967", "ID": "rsY",
            "ALLELE0": "G", "ALLELE1": "A",
            "A1FREQ": "0.1", "INFO": "0.9", "N": "33628", "TEST": "ADD",
            "BETA": "0.0", "SE": "0.05", "CHISQ": "0.01",
            "LOG10P": "0.01", "EXTRA": "0",
        },
        # Out-of-window (above)
        {
            "CHROM": "1", "GENPOS": "109774969", "ID": "rsZ",
            "ALLELE0": "T", "ALLELE1": "C",
            "A1FREQ": "0.1", "INFO": "0.9", "N": "33628", "TEST": "ADD",
            "BETA": "0.0", "SE": "0.05", "CHISQ": "0.01",
            "LOG10P": "0.01", "EXTRA": "0",
        },
        # Different chromosome (should not appear)
        {
            "CHROM": "2", "GENPOS": "109274968", "ID": "rsOtherChrom",
            "ALLELE0": "T", "ALLELE1": "G",
            "A1FREQ": "0.32", "INFO": "0.999", "N": "33628", "TEST": "ADD",
            "BETA": "0.0", "SE": "0.05", "CHISQ": "0.1",
            "LOG10P": "0.1", "EXTRA": "0",
        },
    ]
    return _build_test_tar(tmp_path, chrom="1", rows=rows)


def test_extract_region_filters_to_window(synthetic_tar):
    """Only in-window variants on the requested chromosome appear."""
    variants, notes = _extract_region_from_tar(
        tar_path=synthetic_tar,
        chromosome="1",
        start_bp=108_774_968,
        end_bp=109_774_968,
    )
    positions = sorted(v.position for v in variants)
    assert positions == [109_274_968, 109_274_970]
    assert all(v.chromosome == "1" for v in variants)
    assert notes == []


def test_extract_region_returns_empty_when_window_misses(synthetic_tar):
    """Window entirely outside the source range returns no rows + no error."""
    variants, _ = _extract_region_from_tar(
        tar_path=synthetic_tar,
        chromosome="1",
        start_bp=1, end_bp=100,
    )
    assert variants == []


def test_extract_region_chr_filter_is_strict(synthetic_tar):
    """chr1 query should not match the chr2 row in the same tar member."""
    variants, _ = _extract_region_from_tar(
        tar_path=synthetic_tar,
        chromosome="1",
        start_bp=109_000_000, end_bp=109_500_000,
    )
    assert all(v.chromosome == "1" for v in variants)


def test_extract_region_raises_when_chr_member_missing(tmp_path):
    """Tar without a chr<N> member surfaces a clear error."""
    bad_tar = tmp_path / "empty.tar"
    with tarfile.open(bad_tar, mode="w") as tf:
        info = tarfile.TarInfo(name="README.txt")
        info.size = 0
        tf.addfile(info, io.BytesIO(b""))
    with pytest.raises(UKBPPPAccessError, match="no member matching chr1"):
        _extract_region_from_tar(
            tar_path=bad_tar, chromosome="1",
            start_bp=1, end_bp=1000,
        )


def test_extract_region_word_boundary_avoids_chr1_chr10_collision(tmp_path):
    """A chr<N> query for chr1 should not accidentally match chr10/11/etc.

    Verifies the strict word-boundary regex that distinguishes the
    chr<N> token from longer-number neighbours inside the same tar.
    """
    chr10_only_tar = tmp_path / "chr10only.tar"
    payload = (
        " ".join(REGENIE_COLUMNS) + "\n"
        + "10 1000 rs0 A T 0.1 0.95 100 ADD 0.0 0.05 1 0.1 _" + "\n"
    )
    gz_blob = gzip.compress(payload.encode("utf-8"))
    with tarfile.open(chr10_only_tar, mode="w") as tf:
        info = tarfile.TarInfo(name="discovery_chr10_PROT_X_v1_Panel.regenie.gz")
        info.size = len(gz_blob)
        tf.addfile(info, io.BytesIO(gz_blob))
    with pytest.raises(UKBPPPAccessError, match="no member matching chr1"):
        _extract_region_from_tar(
            tar_path=chr10_only_tar, chromosome="1",
            start_bp=1, end_bp=2000,
        )


# ----------------------- Client end-to-end with mocks -----------------------


def test_live_fetch_n_samples_comes_from_rows(synthetic_tar, stub_lister, tmp_path):
    """End-to-end: client resolves protein, downloads tar, extracts region.

    `bundled_slices_dir` is pointed at an empty tmp dir so the client
    falls through to the live-fetch (mocked) path, bypassing any real
    bundled slices that might exist in the skill's bundled_slices/ dir.
    """
    stub_lister.list_folder.return_value = [
        {"name": "SORT1_Q99523_OID20684_v1_Cardiometabolic.tar",
         "id": "syn52362400"},
    ]
    downloader = MagicMock()
    downloader.download_tar = MagicMock(return_value=synthetic_tar)
    client = UKBPPPClient(
        lister=stub_lister, downloader=downloader,
        bundled_slices_dir=tmp_path / "empty_bundle",
    )

    result = client.fetch_region(
        protein_label="SORT1", ancestry="EUR",
        chromosome="1", start_bp=108_774_968, end_bp=109_774_968,
    )

    assert isinstance(result, RegionResult)
    assert result.n_variants == 2
    assert result.release.protein_hgnc == "SORT1"
    assert result.release.protein_uniprot == "Q99523"
    assert result.release.olink_reagent_id == "OID20684"
    assert result.release.ancestry == "EUR"
    assert result.release.ancestry_label == "European (discovery)"
    # N is the rows' own REGENIE N, not a per-ancestry constant; the
    # published cohort size rides beside it.
    assert result.release.n_samples == 33_628
    assert result.release.n_samples_min == 33_628
    assert result.release.n_samples_label == "33,628"
    assert result.release.cohort_n == 34_557
    assert result.release.synapse_id == "syn52362400"
    assert "syn52362400" in result.release.source_url
    # Verify downloader received the right synapse id
    downloader.download_tar.assert_called_once_with("syn52362400")


def _two_row_live_result(tmp_path, stub_lister, n_values: tuple[str, str]):
    """Mocked live fetch over a tar with two in-window rows carrying `n_values`."""
    base = {
        "CHROM": "1", "ID": "rs", "ALLELE0": "T", "ALLELE1": "G",
        "A1FREQ": "0.3", "INFO": "0.99", "TEST": "ADD", "BETA": "0.1",
        "SE": "0.05", "CHISQ": "4", "LOG10P": "1.3", "EXTRA": "0",
    }
    rows = [
        {**base, "GENPOS": "109274968", "N": n_values[0]},
        {**base, "GENPOS": "109274970", "N": n_values[1]},
    ]
    tar = _build_test_tar(tmp_path, chrom="1", rows=rows)
    stub_lister.list_folder.return_value = [
        {"name": "SORT1_Q99523_OID20684_v1_Cardiometabolic.tar", "id": "syn1"},
    ]
    downloader = MagicMock()
    downloader.download_tar = MagicMock(return_value=tar)
    client = UKBPPPClient(
        lister=stub_lister, downloader=downloader,
        bundled_slices_dir=tmp_path / "empty_bundle",
    )
    return client.fetch_region(
        protein_label="SORT1", ancestry="EUR",
        chromosome="1", start_bp=108_774_968, end_bp=109_774_968,
    )


def test_n_samples_is_a_range_when_rows_disagree(stub_lister, tmp_path):
    """REGENIE N is per variant; when the rows disagree the release carries
    the range (max as `n_samples`, min beside it) and the label says so."""
    result = _two_row_live_result(tmp_path, stub_lister, ("33000", "33628"))
    assert result.n_variants == 2
    assert result.release.n_samples == 33_628
    assert result.release.n_samples_min == 33_000
    assert result.release.n_samples_label == "33,000 to 33,628 (varies by variant)"


def test_n_samples_none_when_rows_state_no_n(stub_lister, tmp_path):
    """Rows with no usable N leave the N unknown (None), never a table value;
    the rows themselves are kept and the published cohort size still rides."""
    result = _two_row_live_result(tmp_path, stub_lister, ("NA", "NA"))
    assert result.n_variants == 2
    assert result.release.n_samples is None
    assert result.release.n_samples_min is None
    assert result.release.n_samples_label is None
    assert result.release.cohort_n == 34_557


_N_CELL = st.one_of(
    st.integers(min_value=0, max_value=60_000).map(str),
    st.sampled_from(["NA", "", "n/a", "3.5e4", "33628.0", "33.5", "-5", "nan"]),
)


@given(cells=st.lists(_N_CELL, max_size=20))
def test_n_samples_from_rows_is_the_range_of_the_parsable_ns(cells):
    """Property: the (min, max) is taken over exactly the rows whose N is a
    positive whole number (written as an integer or a float); no such N gives
    (None, None)."""
    variants = [
        RegionVariant(
            variant_id=f"1_{i}_A_C", chromosome="1", position=i, ref="A",
            alt="C", beta=None, se=None, p_value=None, maf=None,
            effect_allele_frequency=None, raw={"N": c},
        )
        for i, c in enumerate(cells)
    ]
    def whole_positive(c):
        try:
            f = float(c)
        except ValueError:
            return None
        return int(f) if f == f and f > 0 and f == int(f) else None

    parsed = [n for n in map(whole_positive, cells) if n is not None]
    lo, hi = _n_samples_from_rows(variants)
    if parsed:
        assert (lo, hi) == (min(parsed), max(parsed))
    else:
        assert (lo, hi) == (None, None)


# ----------------------- Bundled-slice loader -----------------------


def _sort1_eur_test_result(chrom="1", start=108_774_968, end=109_774_968):
    release = UKBPPPRelease(
        release_label="UKB-PPP r1 2023 (Sun 2023)",
        fetched_at_utc="2026-05-15T14:30:00Z",
        protein_hgnc="SORT1",
        protein_uniprot="Q99523",
        olink_reagent_id="OID20684",
        olink_panel="Cardiometabolic",
        ancestry="EUR",
        ancestry_label="European (discovery)",
        n_samples=33_628,
        synapse_id="syn52362400",
        source_url="https://www.synapse.org/Synapse:syn52362400",
        n_samples_min=33_628,
        cohort_n=34_557,
    )
    variants = [
        RegionVariant(
            variant_id="1_109274968_T_G", chromosome=chrom, position=109_274_968,
            ref="T", alt="G", beta=0.44, se=0.05, p_value=3.16e-13,
            maf=0.32, effect_allele_frequency=0.32, raw={"N": "33628"},
        ),
    ]
    return RegionResult(
        protein_label_short="SORT1", ancestry="EUR",
        chromosome=chrom, region_start_bp=start, region_end_bp=end,
        n_variants=len(variants), variants=variants, release=release,
    )


def test_bundled_slice_key_is_deterministic():
    """The on-disk filename for a (protein, ancestry, region) is stable
    across case + whitespace + `chr` prefix variants."""
    a = _bundled_slice_key("SORT1", "EUR", "1", 108_774_968, 109_774_968)
    b = _bundled_slice_key(" sort1 ", "eur", "chr1", 108_774_968, 109_774_968)
    assert a == b == "SORT1__EUR__chr1__108774968_109774968.json.gz"


def test_save_and_load_round_trip(tmp_path):
    """Writing a RegionResult to bundled_slices/ and reading it back via
    `fetch_region` returns an equivalent result (no network)."""
    original = _sort1_eur_test_result()
    save_region_result_as_bundled_slice(
        original, bundled_slices_dir=tmp_path,
    )

    client = UKBPPPClient(
        lister=MagicMock(),       # not used
        downloader=MagicMock(),   # not used
        bundled_slices_dir=tmp_path,
        allow_live_fetch=False,   # prove no network is touched
    )
    loaded = client.fetch_region(
        protein_label="SORT1", ancestry="EUR",
        chromosome="1", start_bp=108_774_968, end_bp=109_774_968,
    )
    assert loaded.n_variants == 1
    assert loaded.variants[0].variant_id == "1_109274968_T_G"
    assert loaded.release.protein_hgnc == "SORT1"
    assert loaded.release.synapse_id == "syn52362400"
    assert loaded.release.n_samples == 33_628
    assert loaded.release.n_samples_min == 33_628
    assert loaded.release.cohort_n == 34_557


def test_slice_saved_before_0_1_2_is_rederived_from_its_rows(tmp_path):
    """A slice written before 0.1.2 stores a per-ancestry table value as its N
    and has no `n_samples_min`. The loader derives the N from the rows instead,
    so a caller's own older slice never shows the retired number under the
    "per-variant N" label."""
    path = save_region_result_as_bundled_slice(
        _sort1_eur_test_result(), bundled_slices_dir=tmp_path,
    )
    with gzip.open(path, "rt", encoding="utf-8") as f:
        d = json.load(f)
    d["release"]["n_samples"] = 46_673
    for key in ("n_samples_min", "cohort_n"):
        d["release"].pop(key, None)
    with gzip.open(path, "wt", encoding="utf-8") as f:
        json.dump(d, f)
    r = _load_bundled_slice(path)
    assert r.release.n_samples == 33_628
    assert r.release.n_samples_min == 33_628
    assert r.release.cohort_n == 34_557


def test_maybe_int_reads_only_positive_whole_numbers():
    assert _maybe_int("33628") == 33_628
    assert _maybe_int("33628.0") == 33_628
    assert _maybe_int(33_628) == 33_628
    for bad in ("", "NA", None, "0", "-5", "33.5", "nan", "n/a"):
        assert _maybe_int(bad) is None, bad


def test_every_shipped_slice_stores_the_rows_n():
    """Every bundled slice's stored N equals the per-variant N of its own rows.
    The loader trusts the stored block, so this is the invariant that keeps it
    honest: a slice saved with a stale or invented N fails here."""
    slices = sorted(_bundled_slices_dir().glob("*.json.gz"))
    assert slices, "no bundled slices found"
    for path in slices:
        r = _load_bundled_slice(path)
        lo, hi = _n_samples_from_rows(r.variants)
        assert hi is not None, path.name
        assert r.release.n_samples == hi, path.name
        assert r.release.n_samples_min == lo, path.name
        assert r.release.cohort_n == UKB_PPP_COHORT_N[r.ancestry], path.name
    sort1 = _load_bundled_slice(
        _bundled_slices_dir() / "SORT1__EUR__chr1__108774968_109774968.json.gz"
    )
    assert _n_samples_from_rows(sort1.variants) == (33_628, 33_628)
    assert sort1.release.cohort_n == 34_557


def test_cli_demo_manifest_and_report_carry_both_numbers(tmp_path, capsys):
    """The offline demo reports the rows' N and the published cohort size,
    each labeled, and never the retired table value."""
    assert main(["--demo", "sort1_ukb_ppp_eur", "--output", str(tmp_path)]) == 0
    assert "| N = 33,628" in capsys.readouterr().out
    if (tmp_path / "manifest.yaml").exists():
        import yaml
        manifest = yaml.safe_load((tmp_path / "manifest.yaml").read_text())
    else:
        manifest = json.loads((tmp_path / "manifest.json").read_text())
    rel = manifest["release"]
    assert rel["n_samples"] == 33628
    assert rel["n_samples_min"] == 33628
    assert rel["n_samples_label"] == "33,628"
    assert rel["cohort_n"] == 34557
    assert "37794186" in rel["cohort_n_source"]
    report = (tmp_path / "report.md").read_text()
    assert "per-variant N" in report
    assert "33,628" in report
    assert "34,557" in report
    assert "46,673" not in report


def test_bundled_slice_takes_precedence_over_live_fetch(tmp_path, synthetic_tar):
    """Bundled slice wins over a fully-configured live downloader."""
    save_region_result_as_bundled_slice(
        _sort1_eur_test_result(), bundled_slices_dir=tmp_path,
    )

    lister = MagicMock()
    downloader = MagicMock()
    downloader.download_tar = MagicMock(return_value=synthetic_tar)
    client = UKBPPPClient(
        lister=lister, downloader=downloader,
        bundled_slices_dir=tmp_path,
    )
    loaded = client.fetch_region(
        protein_label="SORT1", ancestry="EUR",
        chromosome="1", start_bp=108_774_968, end_bp=109_774_968,
    )
    # The on-disk slice has 1 variant; the synthetic tar has 2 in-window
    # variants. If the live path ran, n_variants would be 2.
    assert loaded.n_variants == 1
    downloader.download_tar.assert_not_called()
    lister.list_folder.assert_not_called()


def test_no_bundle_no_live_raises_clear_error(tmp_path):
    """`allow_live_fetch=False` + no bundled slice = explicit error,
    NOT a silent Synapse hit."""
    client = UKBPPPClient(
        lister=MagicMock(), downloader=MagicMock(),
        bundled_slices_dir=tmp_path,   # empty
        allow_live_fetch=False,
    )
    with pytest.raises(UKBPPPAccessError, match="no bundled slice"):
        client.fetch_region(
            protein_label="SORT1", ancestry="EUR",
            chromosome="1", start_bp=1, end_bp=1000,
        )


# ----------------------- PAT-required UX -----------------------


def test_pat_required_message_is_loud_and_actionable():
    """Per user direction 2026-05-15 the error must be very explicit:
    name the Synapse signup URL, the PAT page, and the env var."""
    assert "https://www.synapse.org" in PAT_REQUIRED_MESSAGE
    assert "Profile:settings" in PAT_REQUIRED_MESSAGE
    assert "SYNAPSE_AUTH_TOKEN" in PAT_REQUIRED_MESSAGE
    assert "view" in PAT_REQUIRED_MESSAGE
    assert "download" in PAT_REQUIRED_MESSAGE
    # Clarifies the bundled-slice escape hatch and the no-UKB-application fact.
    assert "bundled" in PAT_REQUIRED_MESSAGE.lower()
    assert "NO UK Biobank Application" in PAT_REQUIRED_MESSAGE
