"""Discrete-MLE model comparison for P(s), the avalanche mass distribution.

Why MLE and not curve-fitting the binned PDF (readme's route):
  * the data ARE an exact integer histogram, so the likelihood is exact:
        ln Lik = sum_s  count[s] * ln p(s)
    with p normalized on the discrete support s = smin..CAP.  No binning choice,
    no weighting choice, no arbitrary residual metric.
  * an explicitly NORMALIZED model kills the degeneracy the readme hit.  In a
    free 4-param fit A*s^-tau*exp(-(s/s0)^b), A absorbs whatever tau and b do,
    so the tau<->b ridge is flat and the fit runs to tau=-9.57.  Under MLE the
    amplitude is not free -- it IS the normalization -- so the ridge is lifted.
"""
import numpy as np
from scipy.optimize import minimize
from scipy.special import gammaln

CAP = 200000          # discrete support ceiling; >> s_max (~10^3), model tails die long before


def _norm_logp(logw, s_all):
    """Normalize an unnormalized log-weight over the support; return log p."""
    mx = logw.max()
    Z = np.exp(logw - mx).sum()
    return logw - mx - np.log(Z)


# --- models: each maps params -> log p(s) on the full support s_all ----------
# parameters are passed in unconstrained (log) form where positivity is needed.

def m_pl(th, s):
    tau, = th
    return -tau * np.log(s)


def m_plexp(th, s):
    tau, lsc = th
    return -tau * np.log(s) - s / np.exp(lsc)


def m_stretch(th, s):
    """pure stretched exponential, no power-law prefactor"""
    ls0, lb = th
    return -(s / np.exp(ls0)) ** np.exp(lb)


def m_weibull(th, s):
    """readme's constrained form: tau = 1-k and stretch exponent both = k"""
    ls0, lk = th
    k = np.exp(lk)
    return (k - 1) * np.log(s) - (s / np.exp(ls0)) ** k


def m_plstretch(th, s):
    """the full 3-param family; nests every model above"""
    tau, ls0, lb = th
    return -tau * np.log(s) - (s / np.exp(ls0)) ** np.exp(lb)


def m_lognorm(th, s):
    mu, lsig = th
    sig = np.exp(lsig)
    return -np.log(s) - (np.log(s) - mu) ** 2 / (2 * sig ** 2)


MODELS = {
    "power law                s^-t":              (m_pl,        ["tau"],              [2.0]),
    "power law x exp          s^-t e^-s/sc":      (m_plexp,     ["tau", "ln sc"],     [1.5, 4.0]),
    "stretched exp            e^-(s/s0)^b":       (m_stretch,   ["ln s0", "ln b"],    [2.0, -1.0]),
    "Weibull (readme, t=1-k)  s^k-1 e^-(s/s0)^k": (m_weibull,   ["ln s0", "ln k"],    [2.0, -1.0]),
    "power law x stretched    s^-t e^-(s/s0)^b":  (m_plstretch, ["tau", "ln s0", "ln b"], [1.5, 2.0, -1.0]),
    "lognormal":                                  (m_lognorm,   ["mu", "ln sig"],     [1.5, 1.0]),
}


def fit(counts_by_s, smin, model, x0):
    """counts_by_s: dict/array giving count at each integer s. Returns (params, lnL, n)."""
    s_all = np.arange(smin, CAP + 1, dtype=float)
    s_obs = np.array(sorted(k for k in counts_by_s if k >= smin), dtype=float)
    c_obs = np.array([counts_by_s[int(k)] for k in s_obs], dtype=float)
    idx = (s_obs - smin).astype(int)
    n = c_obs.sum()

    def nll(th):
        if not np.all(np.isfinite(th)):
            return 1e18
        lw = model(th, s_all)
        if not np.all(np.isfinite(lw)):
            return 1e18
        lp = _norm_logp(lw, s_all)
        return -(c_obs * lp[idx]).sum()

    best = None
    for jit in range(6):
        start = np.array(x0, dtype=float) * (1 + 0.25 * jit * np.random.RandomState(jit).randn(len(x0)))
        r = minimize(nll, start, method="Nelder-Mead",
                     options=dict(maxiter=40000, maxfev=40000, xatol=1e-9, fatol=1e-6))
        r = minimize(nll, r.x, method="Nelder-Mead",
                     options=dict(maxiter=40000, maxfev=40000, xatol=1e-10, fatol=1e-8))
        if best is None or r.fun < best.fun:
            best = r
    return best.x, -best.fun, n


def worst_factor(counts_by_s, smin, model, th):
    """max multiplicative discrepancy between model and data, on log-binned P(s)."""
    s_all = np.arange(smin, CAP + 1, dtype=float)
    lp = _norm_logp(model(th, s_all), s_all)
    p = np.exp(lp)
    s_obs = np.array(sorted(k for k in counts_by_s if k >= smin), dtype=float)
    c_obs = np.array([counts_by_s[int(k)] for k in s_obs], dtype=float)
    n = c_obs.sum()
    # log bins so tail bins hold real counts
    edges = np.unique(np.round(np.logspace(np.log10(smin), np.log10(s_obs.max() + 1), 25)).astype(int))
    worst, where = 1.0, None
    for a, b in zip(edges[:-1], edges[1:]):
        sel = (s_obs >= a) & (s_obs < b)
        if not sel.any():
            continue
        obs = c_obs[sel].sum() / n
        exp = p[(s_all >= a) & (s_all < b)].sum()
        if obs <= 0 or exp <= 0:
            continue
        r = max(obs / exp, exp / obs)
        if r > worst:
            worst, where = r, (a, b)
    return worst, where
