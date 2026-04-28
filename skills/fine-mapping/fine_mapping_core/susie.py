"""
susie.py — Pure-Python SuSiE (Sum of Single Effects) fine-mapping.

Implements the Iterative Bayesian Stepwise Selection (IBSS) algorithm from:
    Wang et al. (2020) JRSS-B doi:10.1111/rssb.12388

This is a pure-numpy implementation that requires no R or external SuSiE
package. It matches the core algorithm but omits some advanced features
(e.g. SuSiE-inf, intercept estimation).
"""

from __future__ import annotations

import warnings

import numpy as np
from .abf import _log_abf, DEFAULT_W


def run_susie(
    z: np.ndarray,
    R: np.ndarray,
    n: int,
    L: int = 10,
    w: float = DEFAULT_W,
    max_iter: int = 100,
    tol: float = 1e-3,
    min_purity: float = 0.5,
    null_weight: float | None = None,
) -> dict:
    """Run SuSiE fine-mapping.

    Parameters
    ----------
    z : (p,) z-score vector
    R : (p, p) LD correlation matrix
    n : effective sample size
    L : maximum number of causal signals
    w : prior variance on each single effect (Wakefield W)
    max_iter : maximum IBSS iterations
    tol : ELBO convergence tolerance
    min_purity : minimum pairwise |r| within a credible set (Wang 2020 section 3.2)
    null_weight : prior weight on the null hypothesis (no effect) for each
        single-effect regression. When > 0, the model can assign posterior mass
        to "no effect at this locus", preventing phantom PIPs on null loci.
        The susieR reference implementation uses null_weight to mitigate
        forced signal assignment. Default: 1/(L+1), giving equal prior to
        "no effect" as to each of L possible effects. Set to 0 to disable.

    Returns
    -------
    dict with keys:
        alpha   : (L, p) posterior weight matrix (each row sums to 1)
        mu      : (L, p) posterior mean effect sizes
        mu2     : (L, p) posterior second moments
        pip     : (p,) posterior inclusion probabilities
        null_weight_used : float, the null_weight that was applied
        elbo    : list of ELBO values per iteration
        converged : bool
        n_iter  : int

    Raises
    ------
    ValueError : if n <= 0, w <= 0, or z contains NaN
    """
    if n <= 0:
        raise ValueError("Sample size n must be positive, got %d" % n)
    if w <= 0:
        raise ValueError("Prior variance w must be positive, got %s" % w)
    if np.any(np.isnan(z)):
        raise ValueError("z-score vector contains NaN values")

    # Default null_weight: 1/(L+1) gives equal prior to "no effect" as to
    # each of L possible effects. Prevents forced signal on null loci.
    if null_weight is None:
        null_weight = 1.0 / (L + 1)

    p = len(z)
    z = z.astype(float)
    R = R.astype(float)

    # Variance of z-scores ≈ 1/n (used to derive V_i = 1/n for all variants)
    V = np.full(p, 1.0 / n)

    # Null component: log prior odds of "no effect" vs uniform prior on variants
    # When null_weight > 0, each single-effect regression includes a null
    # hypothesis that competes with the p variant hypotheses.
    use_null = null_weight > 0
    if use_null:
        # log prior: log(null_weight) for null, log((1 - null_weight)/p) for each variant
        log_prior_null = np.log(null_weight)
        log_prior_variant = np.log((1.0 - null_weight) / p)
    else:
        log_prior_variant = -np.log(p)  # uniform 1/p

    # Initialise
    alpha = np.ones((L, p)) / p       # posterior weights (uniform init)
    mu    = np.zeros((L, p))           # alpha-weighted posterior means (for IBSS updates)
    mu2   = np.zeros((L, p))           # alpha-weighted posterior second moments
    cond_mu = np.zeros((L, p))        # conditional posterior means (pure, not alpha-weighted)

    elbo_history = []
    converged = False

    for iteration in range(max_iter):
        alpha_prev = alpha.copy()

        # Precompute total fitted effect for residual updates
        fitted_all = (alpha * mu).sum(axis=0)  # (p,)

        for l in range(L):
            # Residual z-score: remove all other effects from z
            # r_l = z - R @ sum_{l' != l} alpha_{l'} * mu_{l'}
            other = fitted_all - alpha[l] * mu[l]  # shape (p,)
            r_l = z - R @ other  # shape (p,)

            # Single-effect regression: compute log ABF for each variant
            # treating r_l as observed z-score with variance V
            log_bf = _log_abf(r_l, V, w)

            # Add prior: log_bf + log_prior_variant for each variant
            log_posterior = log_bf + log_prior_variant

            if use_null:
                # Null component: BF = 1 (log BF = 0), prior = null_weight
                log_null_posterior = log_prior_null  # log(null_weight) + 0

                # Normalise across p variants + 1 null
                all_log_posts = np.append(log_posterior, log_null_posterior)
                max_lp = np.max(all_log_posts)
                all_log_posts_shifted = all_log_posts - max_lp
                all_weights = np.exp(all_log_posts_shifted)
                total = all_weights.sum()

                # alpha[l] gets the variant weights (excluding null)
                alpha[l] = all_weights[:p] / total
                # null_alpha is all_weights[p] / total (not stored, just absorbed)
            else:
                # Original behaviour: normalise across variants only
                log_bf_shifted = log_bf - np.max(log_bf)
                alpha[l] = np.exp(log_bf_shifted) / np.exp(log_bf_shifted).sum()

            # Posterior mean and second moment (Gaussian single-effect)
            # mu_l_j  = w / (V_j + w) * r_l_j    (scalar approximation per variant)
            post_mean_j = (w / (V + w)) * r_l          # conditional posterior mean per variant
            post_var_j  = w * V / (V + w)              # conditional posterior variance per variant

            cond_mu[l] = post_mean_j                   # pure conditional posterior mean
            mu[l]  = alpha[l] * post_mean_j            # alpha-weighted (for IBSS fitted values)
            mu2[l] = alpha[l] * (post_mean_j**2 + post_var_j)

            # Keep fitted_all in sync for the next effect's residual
            fitted_all = (alpha * mu).sum(axis=0)

        # ELBO (approximate): use KL between current and previous alpha
        elbo = _compute_elbo(z, R, alpha, mu, mu2, V, w)
        elbo_history.append(elbo)

        if iteration > 0 and abs(elbo_history[-1] - elbo_history[-2]) < tol:
            converged = True
            break

    # When null component is active, prune null effects: if an effect row's
    # maximum alpha is below the uniform prior (1/p), the null hypothesis
    # dominates and the effect should not contribute to PIPs. This prevents
    # phantom PIP accumulation from multiple null effects on a null locus.
    if use_null:
        uniform_prior = 1.0 / p
        active_mask = alpha.max(axis=1) > uniform_prior
        alpha_active = alpha[active_mask] if active_mask.any() else np.zeros((0, p))
    else:
        alpha_active = alpha

    # Compute PIPs: PIP_i = 1 - prod_l (1 - alpha_{l,i})  [active effects only]
    if alpha_active.shape[0] > 0:
        pip = 1.0 - np.prod(1.0 - alpha_active, axis=0)
    else:
        pip = np.zeros(p)
    pip = np.clip(pip, 0.0, 1.0)

    # Non-convergent results: raise ValueError to prevent downstream use of
    # unreliable estimates. This matches susieR behavior of warning AND
    # setting $converged=FALSE. Raising ensures callers cannot silently
    # consume non-converged PIPs.
    if not converged:
        raise ValueError(
            f"SuSiE IBSS did not converge after {iteration + 1} iterations "
            f"(tol={tol}). Increase max_iter or check input data quality."
        )

    return {
        "alpha": alpha,
        "mu": cond_mu,          # conditional posterior mean (pure, per Wang et al. eq. 4)
        "mu_weighted": mu,      # alpha-weighted posterior mean (used in IBSS fitted values)
        "mu2": mu2,
        "pip": pip,
        "null_weight_used": null_weight,
        "elbo": elbo_history,
        "converged": converged,
        "n_iter": iteration + 1,
    }


def _compute_elbo(
    z: np.ndarray,
    R: np.ndarray,
    alpha: np.ndarray,
    mu: np.ndarray,
    mu2: np.ndarray,
    V: np.ndarray,
    w: float,
) -> float:
    """Approximate ELBO for convergence monitoring.

    Uses the expected log-likelihood minus KL divergence.
    This is a simplified scalar approximation sufficient for convergence checks.
    """
    L, p = alpha.shape
    n = 1.0 / V[0]  # approximate

    # Expected fitted values
    fitted = (alpha * mu).sum(axis=0)  # (p,)
    residual = z - R @ fitted
    ell = -0.5 * n * float(residual @ residual)

    # KL: sum_l sum_j alpha[l,j] * (log alpha[l,j] - log(1/p))
    with np.errstate(divide="ignore", invalid="ignore"):
        log_alpha = np.where(alpha > 0, np.log(alpha), 0.0)
    kl = float(np.sum(alpha * (log_alpha - np.log(1.0 / p))))

    return ell - kl
