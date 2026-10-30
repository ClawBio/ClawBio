"""Live smoke tests: the real EBI FTP, gated by RUN_LIVE_TESTS=1 (or `pytest -m live`).

These exist because the unit tests mock pysam, and a skill whose only tests are mocked
can stay broken against its source for months with a green suite: that is what happened
when the catalogue retired its metadata API in September 2026 (eQTL-Catalogue-resources#59).
Network-dependent by design; not part of every-PR CI.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from eqtl_catalogue_region_fetch import (  # noqa: E402
    EQTLCatalogueClient,
    ftp_url_for,
    load_dataset_index,
)

pytestmark = pytest.mark.live


def test_bundled_demo_region_fetches_from_the_ftp():
    """The default demo (SORT1, GTEx minor salivary gland, QTD000276): resolution from the
    bundled index, then a real tabix range fetch. No result cache is involved here."""
    result = EQTLCatalogueClient().fetch_region(
        dataset_id="QTD000276", chromosome="1",
        start_bp=108_774_968, end_bp=109_774_968,
        molecular_trait_id="ENSG00000134243",
    )
    assert result.n_variants > 1000, result.n_variants
    assert result.release.dataset_release == "r7"
    assert result.release.study_label == "GTEx"
    v = result.variants[0]
    assert v.chromosome == "1" and v.beta is not None and v.se is not None


def test_bundled_index_agrees_with_the_ftp_for_the_odd_one_out():
    """QTD000584 (aptamer) is the one r7 dataset where the file the catalogue's table lists
    (`.all`) differs from what the old quant-method rule would pick (`.cc`). The FTP must
    serve the listed file's tabix index. It serves the `.cc` one too for this dataset
    (probed 2026-09-13 on 33 datasets: every `.all`-listed one also has `.cc`, no
    `.cc`-listed one has `.all`), so the point of reading the class from the table is not that the other file
    is always missing; it is that where both exist they hold different rows, and the table
    names the one the catalogue publishes as the dataset's per-variant sumstats."""
    import urllib.request

    row = load_dataset_index()["QTD000584"]
    assert row["file_class"] == "all"

    def status(url: str) -> int:
        req = urllib.request.Request(url, method="HEAD")
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                return r.status
        except urllib.error.HTTPError as e:
            return e.code

    listed = ftp_url_for(row["study_id"], "QTD000584", file_class="all") + ".tbi"
    assert status(listed) == 200, listed


def test_duplicate_rsid_rows_are_collapsed_on_the_sort1_window():
    """The catalogue writes one row per rsID, so a variant with several rsIDs comes
    back on several rows. Pins the measurement to the source: the SORT1 QTD000276
    1 Mb window read 2,833 rows for 2,740 (trait, variant) keys on 2026-10-22. If the
    counts move, the catalogue re-released this dataset; re-measure before editing."""
    result = EQTLCatalogueClient().fetch_region(
        dataset_id="QTD000276", chromosome="1",
        start_bp=108_774_968, end_bp=109_774_968,
        gene_id="ENSG00000134243",
    )
    assert result.n_rows_read == 2833
    assert result.n_variants == len(result.variants) == 2740
    assert len(result.collapsed) == 93
    by_id = {v.variant_id: v for v in result.variants}
    assert by_id["1_108778690_T_TTATA"].raw["rsids"] == "rs35423568;rs5776947"
