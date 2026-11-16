"""SuSiE-inf with est_tausq=True, checked against the reference implementation.

The expected numbers come from the gentropy SuSiE-inf port (opentargets/gentropy,
itself a port of FinucaneLab/fine-mapping-inf), as recorded by clawbio_bench
(biostochastics/clawbio_bench at 49a0f11, MIT) in
test_cases/finemapping/fm_{18,20,21}_*/derived/oracle_expected.json. The loci are
rebuilt here from the prose description in each case's ground_truth.txt; no bench
file is copied.
"""

import sys
from pathlib import Path

import numpy as np
import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

SKILL_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(SKILL_DIR))

from fine_mapping_core import susie_inf
from fine_mapping_core.susie_inf import run_susie_inf

BUMPS = (17, 123)


def _two_block_locus() -> tuple[np.ndarray, np.ndarray]:
    """200 variants in two LD blocks (r = 0.20 within a block, 0 between).

    z is +1.5 across block 1 and -1.5 across block 2, with one outlier per block
    (+3.5 at index 17, -3.5 at index 123): a diffuse signal that an infinitesimal
    component should absorb.
    """
    p, half, rho = 200, 100, 0.20
    R = np.zeros((p, p))
    R[:half, :half] = rho
    R[half:, half:] = rho
    np.fill_diagonal(R, 1.0)
    z = np.where(np.arange(p) < half, 1.5, -1.5)
    z[17], z[123] = 3.5, -3.5
    return z, R


def _fit_two_block(**kwargs) -> dict:
    z, R = _two_block_locus()
    return run_susie_inf(
        z=z, R=R, n=5000, L=5, ssq_init=0.04, null_weight=1 / 6, max_iter=200, **kwargs
    )


def test_est_tausq_estimate_is_applied_and_matches_the_reference():
    """tau^2 of order 1e-5 is a real estimate here, and the fit has to use it."""
    out = _fit_two_block(est_tausq=True)

    assert out["tausq"] == pytest.approx(5.464108383754903e-05, rel=1e-4)
    assert out["sigmasq"] == pytest.approx(0.9890775122196939, abs=1e-6)
    for i in BUMPS:
        assert out["pip"][i] == pytest.approx(0.02167543, abs=1e-5)
    others = np.delete(out["pip"], BUMPS)
    assert others == pytest.approx(0.02061914, abs=1e-5)


def test_est_tausq_off_leaves_tausq_at_zero_and_matches_the_reference():
    """The same locus without the infinitesimal term: the outliers keep their PIP."""
    out = _fit_two_block(est_tausq=False)

    assert out["tausq"] == 0.0
    assert out["sigmasq"] == pytest.approx(0.9988407506195481, abs=1e-6)
    for i in BUMPS:
        assert out["pip"][i] == pytest.approx(0.15426454, abs=1e-5)
    others = np.delete(out["pip"], BUMPS)
    assert others == pytest.approx(0.01956324, abs=1e-5)


def test_est_tausq_changes_the_fit_on_a_diffuse_locus():
    on, off = _fit_two_block(est_tausq=True), _fit_two_block(est_tausq=False)

    assert on["tausq"] > 0.0
    for i in BUMPS:
        assert on["pip"][i] < off["pip"][i] / 5


def test_null_locus_with_est_tausq_falls_back_to_zero_and_keeps_every_effect():
    """A non-positive moment solution means tau^2 = 0, and the PIP still sums over all L effects.

    Five independent variants with |z| <= 0.12, L = 3, null weight 0.25: each effect
    stays near its prior, so every PIP is close to 1 - (1 - 0.75 / 5)^3 = 0.386.
    """
    z = np.array([0.1, -0.05, 0.08, -0.12, 0.03])
    out = run_susie_inf(
        z=z,
        R=np.eye(5),
        n=5000,
        L=3,
        ssq_init=0.04,
        null_weight=0.25,
        max_iter=200,
        est_tausq=True,
    )

    assert out["tausq"] == 0.0
    assert out["sigmasq"] == pytest.approx(1.0000128033958564, abs=1e-6)
    expected = [0.38469885, 0.38466425, 0.38468224, 0.38471915, 0.38465687]
    assert out["pip"] == pytest.approx(expected, abs=1e-5)


@settings(max_examples=60, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(
    z=st.lists(
        st.floats(min_value=-6.0, max_value=6.0, allow_nan=False),
        min_size=3,
        max_size=10,
    ),
    rho=st.floats(min_value=0.0, max_value=0.9),
    n_effects=st.integers(min_value=1, max_value=4),
    est_tausq=st.booleans(),
)
def test_pip_is_aggregated_over_every_effect(z, rho, n_effects, est_tausq):
    """PIP_j = 1 - prod_l (1 - alpha_jl) over all L effects (Wang et al. 2020),
    whatever the locus looks like and whether or not tau^2 is estimated."""
    p = len(z)
    R = np.full((p, p), rho)
    np.fill_diagonal(R, 1.0)

    out = run_susie_inf(z=np.array(z), R=R, n=5000, L=n_effects, est_tausq=est_tausq)

    expected = 1.0 - np.prod(1.0 - out["alpha"], axis=1)
    assert np.all(np.isfinite(out["pip"]))
    assert out["pip"] == pytest.approx(np.clip(expected, 0.0, 1.0), abs=1e-12)
    assert out["tausq"] >= 0.0


NULL_Z = np.array([0.1, -0.05, 0.08, -0.12, 0.03])


def test_non_positive_moment_solution_resets_a_positive_starting_tausq():
    """Started from tau^2 > 0 on a locus with no diffuse signal, the fit ends at tau^2 = 0.

    The reference sets tau^2 to 0 whenever the moment solution is not positive; keeping
    the starting value would report an infinitesimal component the data never supported.
    """
    out = run_susie_inf(
        z=NULL_Z,
        R=np.eye(5),
        n=5000,
        L=3,
        ssq_init=0.04,
        null_weight=0.25,
        max_iter=200,
        est_tausq=True,
        tausq=1e-4,
    )

    assert out["tausq"] == 0.0
    assert out["pip"] == pytest.approx(0.3847, abs=1e-3)


def test_unsolvable_moment_system_falls_back_to_zero_tausq(monkeypatch):
    """If the 2x2 moment system cannot be solved, the fit also ends at tau^2 = 0."""

    def refuse(*args, **kwargs):
        raise np.linalg.LinAlgError("singular matrix")

    monkeypatch.setattr(susie_inf.scipy.linalg, "solve", refuse)

    out = run_susie_inf(z=NULL_Z, R=np.eye(5), n=5000, L=2, est_tausq=True, tausq=1e-4)

    assert out["tausq"] == 0.0
    assert np.all(np.isfinite(out["pip"]))
