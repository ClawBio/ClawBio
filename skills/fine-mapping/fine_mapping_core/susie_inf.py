"""
susie_inf.py — SuSiE with infinitesimal background effects.

Implements the SuSiE-inf algorithm from:
    Cui et al. (2024) "Improving fine-mapping by modeling infinitesimal effects"
    Nature Genetics 56, 162–169. https://doi.org/10.1038/s41588-023-01597-3

SuSiE-inf extends standard SuSiE (Wang et al. 2020) by adding a
polygenic infinitesimal component τ² that captures diffuse background
signal not modelled by the L sparse effects. This improves calibration
when the genetic architecture is not fully sparse.

Key difference from susie.py
-----------------------------
Standard SuSiE assumes residual variance Ω⁻¹ = σ² I in z-score space.
SuSiE-inf uses Ω⁻¹ = (τ² · D² + σ² · I) in the LD eigenbasis, where
D² are the eigenvalues of X'X. When τ²→0 the two models are equivalent.

Algorithm
---------
1. Eigendecompose LD:  LD = V diag(d²/n) V'
2. Precompute  Dsq = n · eigvals,  VtXty = V' · (√n · z)
3. Main IBSS loop (L effects):
   a. Compute weighted residual  XtΩr = V · (V'b / var)  where var = τ²Dsq + σ²
   b. Update per-effect prior variance s² via bounded scalar optimisation
   c. Update α (posterior weights), μ (posterior means), ω (precisions)
4. Update σ², τ² by method-of-moments every iteration
5. Convergence on max |ΔPIP| < tol
6. Final PIP = 1 − ∏_l (1 − Σ_j α_{lj})  (collapsed across effects)
"""

from __future__ import annotations

import numpy as np
import scipy.linalg
import scipy.special
from scipy.optimize import minimize_scalar


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def run_susie_inf(
    z: np.ndarray,
    R: np.ndarray,
    n: int,
    L: int = 10,
    meansq: float = 1.0,
    ssq_init: float = 0.2,
    ssq_range: tuple[float, float] = (0.0, 1.0),
    est_ssq: bool = True,
    est_sigmasq: bool = True,
    est_tausq: bool = False,
    sigmasq: float = 1.0,
    tausq: float = 0.0,
    max_iter: int = 100,
    tol: float = 1e-3,
    null_weight: float | None = None,
) -> dict:
    """Run SuSiE-inf fine-mapping.

    Parameters
    ----------
    z : (p,) z-score vector
    R : (p, p) LD correlation matrix  (equal to X'X / n)
    n : effective sample size
    L : number of modelled sparse causal effects
    meansq : average squared magnitude of y  (||y||²/n), default 1
    ssq_init : initial prior effect-size variance for each of the L effects
    ssq_range : (lower, upper) bounds for per-effect s² optimisation
    est_ssq : estimate per-effect prior variances s² by MLE
    est_sigmasq : estimate residual variance σ² by method-of-moments
    est_tausq : estimate infinitesimal variance τ² by method-of-moments
    sigmasq : initial residual variance σ²
    tausq : initial infinitesimal variance τ² (0 = no infinitesimal term)
    max_iter : maximum IBSS iterations
    tol : convergence threshold on max |ΔPIP|

    Returns
    -------
    dict with keys:
        pip      : (p,) posterior inclusion probabilities (collapsed)
        alpha    : (p, L) per-effect posterior weights (each column sums to 1)
        mu       : (p, L) posterior means conditional on being causal
        omega    : (p, L) posterior precisions conditional on being causal
        ssq      : (L,) final per-effect prior variances
        sigmasq  : final residual variance σ²
        tausq    : final infinitesimal variance τ²
        n_iter   : number of iterations run
        converged: bool
    """
    z = np.asarray(z, dtype=float)
    R = np.asarray(R, dtype=float)
    p = len(z)

    # Eigendecompose LD once — O(p³) but done only once
    eigvals, V = scipy.linalg.eigh(R)            # ascending order
    Dsq = np.maximum(n * eigvals, 0.0)           # n · eigenvalues of X'X

    # Precompute projections
    Xty = np.sqrt(n) * z                          # X'y  (z-score convention)
    VtXty = V.T @ Xty                             # V'X'y
    yty = n * meansq                              # y'y

    # Initialise variance diagonal and weighted quantities
    var = tausq * Dsq + sigmasq                   # (p,) diagonal in eigenbasis
    diag_XtOX = np.sum(V ** 2 * (Dsq / var), axis=1)   # diag(X'ΩX)
    XtOy = V @ (VtXty / var)                     # X'Ωy

    # Initialise effects
    ssq = np.full(L, ssq_init, dtype=float)
    alpha = np.ones((p, L)) / p                  # uniform init, (p, L)
    mu = np.zeros((p, L))
    omega = diag_XtOX[:, None] + 1.0 / ssq       # (p, L)
    lbf_variable = np.zeros((p, L))

    # Null weight: prior probability that each single effect is absent.
    # Default 1/(L+1) gives equal prior to "no effect" as to each of L
    # possible effects, consistent with standard SuSiE null handling.
    if null_weight is None:
        null_weight = 1.0 / (L + 1)
    use_null = null_weight > 0
    if use_null:
        log_prior_null = np.log(null_weight)
        log_prior_variant = np.log((1.0 - null_weight) / p)
        logpi0 = np.full(p, log_prior_variant)
    else:
        logpi0 = np.full(p, -np.log(p))              # log uniform prior

    converged = False
    n_iter = 0

    for iteration in range(max_iter):
        alpha_prev = alpha.copy()

        for l in range(L):
            # Weighted residual: remove all effects except l
            b_others = (alpha * mu).sum(axis=1) - alpha[:, l] * mu[:, l]
            XtOXb = V @ (V.T @ b_others * Dsq / var)
            XtOr = XtOy - XtOXb                  # X'Ω r_l

            if est_ssq:
                # Maximise single-effect ELBO over s²_l (bounded scalar)
                def neg_elbo(s: float) -> float:
                    return -scipy.special.logsumexp(
                        XtOr ** 2 / (2.0 * (diag_XtOX + 1.0 / s))
                        - 0.5 * np.log((diag_XtOX + 1.0 / s) * s)
                        + logpi0
                    )

                res = minimize_scalar(neg_elbo, bounds=ssq_range, method="bounded")
                if res.success:
                    ssq[l] = res.x

            # Update posterior quantities for effect l
            omega[:, l] = diag_XtOX + 1.0 / ssq[l]
            mu[:, l] = XtOr / omega[:, l]
            lbf_variable[:, l] = (
                XtOr ** 2 / (2.0 * omega[:, l])
                - 0.5 * np.log(omega[:, l] * ssq[l])
            )
            log_alpha_l = lbf_variable[:, l] + logpi0
            if use_null:
                # Include null hypothesis in normalisation
                all_log = np.append(log_alpha_l, log_prior_null)
                log_norm = scipy.special.logsumexp(all_log)
                alpha[:, l] = np.exp(log_alpha_l - log_norm)
            else:
                alpha[:, l] = np.exp(log_alpha_l - scipy.special.logsumexp(log_alpha_l))

        # Method-of-moments update for σ² (and optionally τ²)
        if est_sigmasq or est_tausq:
            sigmasq, tausq = _mom_update(
                alpha, mu, omega, sigmasq, tausq,
                n, V, Dsq, VtXty, Xty, yty,
                est_sigmasq=est_sigmasq, est_tausq=est_tausq,
            )
            # Only apply tausq to the variance structure if it exceeds a
            # minimum threshold. Tiny tausq estimates (< 1e-3) spread
            # posterior weight across eigenvectors without meaningfully
            # modelling a polygenic background, degrading SER quality.
            effective_tausq = tausq if tausq >= 1e-3 else 0.0
            var = effective_tausq * Dsq + sigmasq
            diag_XtOX = np.sum(V ** 2 * (Dsq / var), axis=1)
            XtOy = V @ (VtXty / var)

        n_iter = iteration + 1
        if np.max(np.abs(alpha - alpha_prev)) < tol:
            converged = True
            break

    # Prune inactive effects: exclude effects with near-zero ssq AND
    # near-uniform alpha from PIP computation. When null_weight is active,
    # the null component in alpha normalization already down-weights
    # inactive effects, so pruning is less critical. When est_tausq is
    # active, use aggressive pruning to prevent background inflation.
    active_mask = np.ones(L, dtype=bool)
    if est_tausq:
        # With infinitesimal component, aggressively prune effects whose
        # ssq has been shrunk below threshold by the optimizer.
        uniform_alpha = 1.0 / p
        ssq_threshold = max(ssq_init * 0.01, 1e-6)
        for l in range(L):
            alpha_max = alpha[:, l].max()
            is_uniform = alpha_max < uniform_alpha * 3.0
            is_tiny_ssq = ssq[l] < ssq_threshold
            if is_uniform and is_tiny_ssq:
                active_mask[l] = False

    # Collapse per-effect alpha into marginal PIPs (active effects only).
    # When no effects are active (all pruned), assign uniform PIPs (1/p)
    # representing maximum uncertainty, consistent with the prior.
    if active_mask.any():
        alpha_active = alpha[:, active_mask]
        pip = 1.0 - np.prod(1.0 - alpha_active, axis=1)
    else:
        pip = np.full(p, 1.0 / p)
    pip = np.clip(pip, 0.0, 1.0)

    return {
        "pip": pip,
        "alpha": alpha,
        "mu": mu,
        "omega": omega,
        "ssq": ssq,
        "sigmasq": float(sigmasq),
        "tausq": float(tausq),
        "n_iter": n_iter,
        "converged": converged,
    }


def _mom_update(
    alpha: np.ndarray,
    mu: np.ndarray,
    omega: np.ndarray,
    sigmasq: float,
    tausq: float,
    n: int,
    V: np.ndarray,
    Dsq: np.ndarray,
    VtXty: np.ndarray,
    Xty: np.ndarray,
    yty: float,
    est_sigmasq: bool,
    est_tausq: bool,
) -> tuple[float, float]:
    """Update σ² (and optionally τ²) via method-of-moments.

    Solves the linear system  A [σ², τ²]' = x  derived from the
    second-order moments of the residual.  If only σ² is estimated,
    uses the first equation directly.
    """
    p, L = mu.shape

    # Build the 2×2 moment matrix A
    A = np.array([[float(n), float(Dsq.sum())],
                  [float(Dsq.sum()), float((Dsq ** 2).sum())]])

    # Compute diag(V' M V) where M is the posterior second-moment matrix
    b = (alpha * mu).sum(axis=1)               # total fitted effect (p,)
    Vtb = V.T @ b
    diag_VtMV = Vtb ** 2                       # start with (V'b)² contribution

    tmpD = np.zeros(p)
    for l in range(L):
        bl = alpha[:, l] * mu[:, l]
        Vtbl = V.T @ bl
        diag_VtMV -= Vtbl ** 2                 # subtract per-effect rank-1
        tmpD += alpha[:, l] * (mu[:, l] ** 2 + 1.0 / omega[:, l])

    # Add posterior variance contribution
    diag_VtMV += np.sum((V.T) ** 2 * tmpD, axis=1)

    # Moment vector x: use Dsq-weighted formulation for correct
    # polygenic background subtraction (Cui et al. 2024, Supplementary eq. 12)
    x = np.zeros(2)
    x[0] = yty - 2.0 * float(b @ Xty) + float(Dsq @ diag_VtMV)
    x[1] = float(VtXty ** 2 @ np.ones(p)) - 2.0 * float(Vtb * Dsq @ VtXty) + float((Dsq ** 2) @ diag_VtMV)

    if est_tausq:
        try:
            sol = scipy.linalg.solve(A, x)
            if sol[0] > 0 and sol[1] > 0:
                sigmasq, tausq = float(sol[0]), float(sol[1])
            else:
                sigmasq = float(x[0] / n)
        except np.linalg.LinAlgError:
            sigmasq = float(x[0] / n)
    elif est_sigmasq:
        new_sigma = (x[0] - A[0, 1] * tausq) / n
        if new_sigma > 0:
            sigmasq = float(new_sigma)

    return sigmasq, tausq


def cred_inf(
    alpha: np.ndarray,
    R: np.ndarray,
    coverage: float = 0.95,
    purity: float = 0.5,
    dedup: bool = True,
) -> list[list[int]]:
    """Build credible sets from the per-effect posterior weight matrix.

    Parameters
    ----------
    alpha : (p, L) per-effect posterior weights from run_susie_inf
    R : (p, p) LD correlation matrix (used for purity filtering)
    coverage : minimum cumulative PIP to include in each credible set
    purity : minimum |r| between all pairs within a credible set;
             sets that fail this threshold are excluded (not returned)
    dedup : deduplicate identical credible sets across effects

    Returns
    -------
    List of credible sets, each a sorted list of 0-based variant indices.
    """
    p, L = alpha.shape
    csets: list[list[int]] = []

    for l in range(L):
        weights = alpha[:, l]
        # Greedy: accumulate from highest weight until coverage reached
        order = np.argsort(-weights)
        cumulative = np.cumsum(weights[order])
        cutoff = int(np.searchsorted(cumulative, coverage)) + 1
        cs_idx = sorted(order[:cutoff].tolist())

        # Purity filter: skip sets with any pair below threshold
        if len(cs_idx) > 1 and purity > 0.0:
            sub = R[np.ix_(cs_idx, cs_idx)]
            # Use a random subsample of 100 when the set is large
            if len(cs_idx) > 100:
                rng = np.random.default_rng(123)
                sample = rng.choice(len(cs_idx), size=100, replace=False)
                sub = sub[np.ix_(sample, sample)]
            min_r = np.min(np.abs(sub - np.diag(np.diag(sub))))
            if min_r < purity:
                continue

        csets.append([int(i) for i in cs_idx])

    if dedup:
        seen: list[tuple[int, ...]] = []
        unique: list[list[int]] = []
        for cs in csets:
            t = tuple(cs)
            if t not in seen:
                seen.append(t)
                unique.append(cs)
        csets = unique

    return csets
