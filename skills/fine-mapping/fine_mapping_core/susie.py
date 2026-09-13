"""
susie.py — SuSiE fine-mapping backed by the sushie package.

Delegates the Sum of Single Effects model (Wang et al. 2020, JRSS-B
doi:10.1111/rssb.12388) to sushie (mancusolab/sushie), a published,
maintained JAX implementation, via its summary-statistics interface
``infer_sushie_ss`` run with a single ancestry. This replaced the previous
hand-rolled numpy IBSS implementation.

Requires the ``fine-mapping`` extra: ``uv sync --extra fine-mapping``.
"""

from __future__ import annotations

import warnings
from importlib.metadata import PackageNotFoundError, version

import numpy as np


def run_susie(
    z: np.ndarray,
    R: np.ndarray,
    n: int,
    L: int = 10,
    w: float = 0.04,
    max_iter: int = 500,
    tol: float = 1e-4,
    min_purity: float = 0.5,
    coverage: float = 0.95,
) -> dict:
    """Run SuSiE fine-mapping on one locus via sushie.

    The keyword set is a contract with the external clawbio_bench harness
    (biostochastics/clawbio_bench, drivers/finemapping_driver.py), which
    imports this function directly and always passes ``w``, ``max_iter``,
    ``tol`` and ``min_purity``. Do not remove or rename them.

    Parameters
    ----------
    z : (p,) z-score vector
    R : (p, p) LD correlation matrix
    n : effective sample size
    L : maximum number of causal signals
    w : prior effect-size variance. Seeds sushie's ``effect_var``; sushie
        then re-estimates it by EM each iteration, so this is a starting
        point rather than a fixed prior as in ABF.
    max_iter : maximum optimization iterations
    tol : ELBO convergence tolerance
    min_purity : minimum pairwise |r| for a credible set to be kept.
        Forwarded to sushie's ``purity`` pruning, so it governs which
        signals survive, not only how they are flagged downstream. sushie
        requires 0 < min_purity < 1.
    coverage : credible-set coverage. Forwarded to sushie's ``threshold``
        so fit-time pruning and downstream set construction use the same
        value. sushie requires 0 < coverage < 1.

    Returns
    -------
    dict with keys:
        alpha     : (k, p) posterior weight matrix over the k kept signals
                    (each row sums to 1)
        mu        : (k, p) conditional posterior mean effect per signal,
                    E[b | included], from sushie ``post_mean``
        mu2       : (k, p) conditional posterior second moment, E[b^2]
                    (not the variance), from sushie ``post_mean_sq``
                    Units: mu/mu2 are on sushie's standardised effect-size
                    scale, NOT susieR ``susie_rss`` z-units — do not compare
                    them to ``r * z`` shrinkage formulas.
        pip       : (p,) posterior inclusion probabilities
        elbo      : list of ELBO values per iteration
        converged : bool
        n_iter    : int
        engine    : "sushie"
        engine_version : installed sushie version string

    Raises
    ------
    ValueError : if n <= 0 or z contains NaN. sushie itself does not reject
        NaN z-scores — it silently returns degenerate PIPs — so we must.
    ImportError : if sushie is not installed (fine-mapping extra missing).
    """
    if n <= 0:
        raise ValueError("Sample size n must be positive, got %d" % n)
    z = np.asarray(z, dtype=float)
    R = np.asarray(R, dtype=float)
    if np.any(np.isnan(z)):
        raise ValueError("z-score vector contains NaN values")
    if R.shape != (len(z), len(z)):
        raise ValueError(
            f"LD matrix shape {R.shape} does not match z length {len(z)}"
        )
    if w <= 0:
        raise ValueError(f"Prior variance w must be positive, got {w}")
    if not 0 < min_purity < 1:
        raise ValueError(
            f"min_purity must satisfy 0 < min_purity < 1 (sushie constraint), got {min_purity}"
        )
    if not 0 < coverage < 1:
        raise ValueError(
            f"coverage must satisfy 0 < coverage < 1 (sushie constraint), got {coverage}"
        )

    try:
        import jax
        from sushie.infer_ss import infer_sushie_ss
    except ImportError as exc:
        raise ImportError(
            "The sushie package is required for SuSiE fine-mapping. "
            "Install it with: uv sync --extra fine-mapping"
        ) from exc

    # float32 (jax default) produces NaN ELBOs on some loci; sushie's own
    # log message recommends enabling x64.
    jax.config.update("jax_enable_x64", True)

    p = len(z)
    result = infer_sushie_ss(
        lds=[R],
        ns=np.array([float(n)]),
        zs=[z],
        L=L,
        # one ancestry → one prior effect variance; EM-updated from here
        effect_var=[float(w)],
        max_iter=max_iter,
        min_tol=tol,
        threshold=coverage,
        purity=min_purity,
        # sushie's guard defaults to 100 common SNPs; keep it for large loci
        # but allow small test/demo loci through.
        min_snps=min(p, 100),
    )

    # Keep only the signals sushie retained as credible sets (coverage
    # threshold + purity pruning). Its alpha rows are ordered so that kept
    # credible set k corresponds to row k-1; the remaining rows are inactive
    # near-uniform effects that would otherwise surface as phantom credible
    # sets downstream. pip_cs is the PIP over kept signals only.
    alpha_all = np.asarray(result.posteriors.alpha, dtype=float)
    kept = sorted(set(result.cs["CSIndex"].to_list())) if result.cs.height else []
    rows = [k - 1 for k in kept]
    alpha = alpha_all[rows] if kept else np.zeros((0, p))
    # Single ancestry: post_mean is (L, p, 1) and post_mean_sq is
    # (L, p, 1, 1); drop the ancestry axes. mu2 is E[b^2], not Var(b).
    mu_all = np.asarray(result.posteriors.post_mean, dtype=float)[..., 0]
    mu2_all = np.asarray(result.posteriors.post_mean_sq, dtype=float)[..., 0, 0]
    mu = mu_all[rows] if kept else np.zeros((0, p))
    mu2 = mu2_all[rows] if kept else np.zeros((0, p))
    pip = np.clip(np.asarray(result.pip_cs, dtype=float), 0.0, 1.0)
    # sushie seeds the ELBO history with a -inf sentinel before iteration 1;
    # strip it so elbo holds only real iterations.
    elbo = [float(e) for e in np.atleast_1d(np.asarray(result.elbo))
            if np.isfinite(e)]
    n_iter = len(elbo)
    converged = bool(result.elbo_increase) and (
        n_iter < max_iter
        or (n_iter >= 2 and abs(elbo[-1] - elbo[-2]) < tol)
    )
    if not converged:
        # Mirror susieR: warn AND flag. PIPs are still returned so the
        # caller can inspect them, but `converged` is the authority.
        warnings.warn(
            f"SuSiE (sushie) did not converge in {n_iter} iterations "
            f"(max_iter={max_iter}, tol={tol}); treat PIPs as provisional.",
            RuntimeWarning,
            stacklevel=2,
        )

    try:
        engine_version = version("sushie")
    except PackageNotFoundError:
        engine_version = "unknown"

    return {
        "alpha": alpha,
        "mu": mu,
        "mu2": mu2,
        "pip": pip,
        "elbo": elbo,
        "converged": converged,
        "n_iter": n_iter,
        "engine": "sushie",
        "engine_version": engine_version,
    }
