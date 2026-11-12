"""The pQTL panel title states the per-variant N, as a range when it varies."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

SKILL_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(SKILL_DIR))

import locuscompare_region_render  # noqa: E402  (puts the sibling skills on sys.path)
from locuscompare_region_render import (  # noqa: E402
    EXPOSURE_KIND_UKB_PPP,
    LocusCompareSpec,
    _build_short_labels_from_spec,
)
from ukb_ppp_region_fetch import UKBPPPRelease  # noqa: E402


def _release(n_max, n_min):
    return UKBPPPRelease(
        release_label="UKB-PPP r1 2023 (Sun 2023)",
        fetched_at_utc="2026-05-15T14:30:00Z",
        protein_hgnc="SORT1",
        protein_uniprot="Q99523",
        olink_reagent_id="OID20684",
        olink_panel="Cardiometabolic",
        ancestry="EUR",
        ancestry_label="European (discovery)",
        n_samples=n_max,
        synapse_id="syn52362400",
        source_url="https://www.synapse.org/Synapse:syn52362400",
        n_samples_min=n_min,
        cohort_n=34_557,
    )


def _spec():
    return LocusCompareSpec(
        lead_variant_id="1_109274968_G_T", chromosome="1",
        lead_position_bp=109_274_968, window_bp=500_000,
        eqtl_dataset_id=None, molecular_trait_id=None,
        gwas_accession="GCST90239655", exposure_gene_symbol="SORT1",
        outcome_trait_label="LDL cholesterol",
        exposure_kind=EXPOSURE_KIND_UKB_PPP,
    )


@pytest.mark.parametrize(
    "n_max, n_min, expected",
    [
        (33_628, 33_628, "N = 33,628 "),
        (33_628, 33_000, "N = 33,000 to 33,628 (varies by variant)"),
    ],
)
def test_pqtl_title_states_the_n_range(n_max, n_min, expected):
    exposure, _ = _build_short_labels_from_spec(
        _spec(), pqtl_release=_release(n_max, n_min),
    )
    assert expected in exposure
