"""P(s), the avalanche-mass distribution at integer N, is ONE fixed shape whose
only N-dependence is a cutoff scale.

    P(s) ~ s^{-tau} exp[-(s/s0)^b],   tau and b FIXED for every N,
                                      only the cutoff s0 = s0(N) moves.

Fractional N is a different system (the partial species is a self-poisoning
impurity that fragments clusters -- see readme's sawtooth section), so this
family is integer N only: N = 6..20 at L=1024, ~1024 sims each (up to 7e9
avalanches at N=6).

METHOD.  Discrete maximum likelihood on the raw integer histogram,

    ln Lik = sum_s count[s] * ln p(s),     p normalized on s = smin..CAP,

is exact (the data ARE an integer histogram) and, with p normalized, keeps the
amplitude from soaking up the exponents.  tau and b are fitted ONCE, shared
across all N (each N contributes only its own scale s0), then frozen; s0(N) is
read off per N.

WHY BOTH EXPONENTS ARE FROZEN.  A per-N free fit of (tau_N, s0_N, b_N) makes
tau AND b appear to rise with N (b even crossing 1).  That is a ridge artifact:
tau, b and s0 trade off inside the abundant bins just above smin, and with ~1e9
counts the fit exploits sub-percent curvature there.  On the GENUINE TAIL the
motion vanishes -- freezing both exponents and moving only s0 costs a negligible
ln-likelihood that collapses as smin is raised (dlnL vs the 3-param-per-N fit:
-5498 at smin=30 -> -206 at 50 -> -37 at 70 -> -11 at 100), while the tail is
tracked to ~1.1x throughout.  So the physical statement is: the tail shape is
universal (a permanent mild stretch, b~0.6); N only shifts the cutoff mass s0
downward -- more species -> matches rarer -> cascades truncate at smaller mass.
The discarded per-N and single-exponent comparisons are in the git history.

Below s ~ 18 a strong even/odd parity oscillation (the elementary cluster is a
pair) means no smooth form applies; fits start at smin=50, past both the parity
band and the near-smin ridge.
"""
import json
import os

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize, minimize_scalar

from common import PLOTS, logbin, pooled

plt.rcParams.update({
    "font.size": 17, "axes.labelsize": 21, "xtick.labelsize": 15,
    "ytick.labelsize": 15, "legend.fontsize": 13, "lines.linewidth": 2.1,
    "figure.autolayout": True,
})

L = 1024
NS = [6, 7, 8, 9, 10, 11, 12, 14, 16, 20]
SMIN = 50                 # past the even/odd parity band and the near-smin ridge
CAP = 6000                # discrete support ceiling; s_max <= 1394, model tail dies far below
SUB = "integerN"
HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, ".integerN_cache.npz")


def save(fig, name):
    d = f"{PLOTS}/{SUB}"
    os.makedirs(d, exist_ok=True)
    fig.savefig(f"{d}/{name}.svg")
    fig.savefig(f"{d}/{name}.png", dpi=150)
    plt.close(fig)
    print(f"saved {SUB}/{name}")


# --------------------------------------------------------------------- data ---
def load_pooled():
    """Pool the s-histograms once (1024 files/N is slow) and cache to disk."""
    if os.path.exists(CACHE):
        d = np.load(CACHE)
        return {N: (d[f"v{N}"], d[f"m{N}"]) for N in NS}
    out, blob = {}, {}
    for N in NS:
        v, h, _ = pooled(L, N)
        out[N] = (v.astype(np.int64), h["mass"].astype(np.int64))
        blob[f"v{N}"], blob[f"m{N}"] = out[N]
        print(f"pooled N={N}")
    np.savez(CACHE, **blob)
    return out


def prep(v, m, smin):
    """Likelihood arrays for one N: support, observed masses, their counts."""
    s_all = np.arange(smin, CAP + 1, dtype=float)
    keep = (m > 0) & (v >= smin)
    s_obs = v[keep].astype(float)
    c_obs = m[keep].astype(float)
    idx = (s_obs - smin).astype(int)
    return dict(s_all=s_all, ls_all=np.log(s_all), s_obs=s_obs, c_obs=c_obs,
                idx=idx, n=c_obs.sum(), smin=smin)


# ----------------------------------------------------------- fixed-shape fit ---
def _lnL(tau, b, s0, D):
    """Normalized discrete log-likelihood of one N at the given shape and scale."""
    w = -tau * D["ls_all"] - (D["s_all"] / s0) ** b
    mx = w.max()
    lp = w - mx - np.log(np.exp(w - mx).sum())
    return (D["c_obs"] * lp[D["idx"]]).sum()


def best_s0(tau, b, D):
    """Per-N cutoff scale that maximizes the likelihood at fixed exponents."""
    r = minimize_scalar(lambda ls0: -_lnL(tau, b, np.exp(ls0), D),
                        bounds=(0.0, 9.0), method="bounded")
    return np.exp(r.x)


def fit(Ds):
    """Fit ONE shared (tau, b) over all N, each N free only in its scale s0."""
    def negtot(x):
        tau, b = x
        if b <= 0:
            return 1e18
        return -sum(_lnL(tau, b, best_s0(tau, b, Ds[N]), Ds[N]) for N in NS)
    r = minimize(negtot, [2.0, 0.6], method="Nelder-Mead",
                 options=dict(xatol=1e-4, fatol=1e-2, maxiter=4000))
    tau, b = r.x
    s0 = {N: best_s0(tau, b, Ds[N]) for N in NS}
    return tau, b, s0


def worst_factor(tau, b, s0, D, mincount=200):
    """Max multiplicative model/data gap on log-binned P(s), tail bins only."""
    s = D["s_all"]
    w = -tau * np.log(s) - (s / s0) ** b
    p = np.exp(w - w.max()); p /= p.sum()
    edges = np.unique(np.round(np.logspace(np.log10(D["smin"]),
                                           np.log10(D["s_obs"].max() + 1), 25)).astype(int))
    worst = 1.0
    for a, bb in zip(edges[:-1], edges[1:]):
        sel = (D["s_obs"] >= a) & (D["s_obs"] < bb)
        cnt = D["c_obs"][sel].sum()
        if cnt < mincount:
            continue
        obs = cnt / D["n"]
        exp = p[(s >= a) & (s < bb)].sum()
        if obs > 0 and exp > 0:
            worst = max(worst, obs / exp, exp / obs)
    return worst


# ------------------------------------------------------------------ figures ---
def fig_pdf(P, Ds, tau, b, s0):
    """P(s) at every integer N with the single frozen-exponent model over the tail."""
    fig, ax = plt.subplots(figsize=(8, 6.5))
    show = [6, 7, 8, 10, 12, 16, 20]
    cols = plt.cm.viridis(np.linspace(0, 0.9, len(show)))
    for N, c in zip(show, cols):
        v, m = P[N]
        x, y = logbin(v, m, nbins=45)
        ax.plot(x, y, "o", ms=4.5, color=c, alpha=0.75,
                label=rf"$N={N}$, $s_0={s0[N]:.0f}$")
        s = Ds[N]["s_all"]
        w = -tau * np.log(s) - (s / s0[N]) ** b
        p = np.exp(w - w.max()); p /= p.sum()
        frac = Ds[N]["n"] / m.sum()
        sel = s <= v[m > 0].max()
        ax.plot(s[sel], p[sel] * frac, "-", color=c, lw=1.6)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"avalanche mass $s$")
    ax.set_ylabel(r"$P(s)$")
    ax.set_ylim(1e-10, 3)
    ax.text(0.03, 0.03, rf"$P(s)\sim s^{{-{tau:.2f}}}\,\exp[-(s/s_0)^{{{b:.2f}}}]$"
            "\n" r"$\tau,b$ fixed; only $s_0(N)$ varies",
            transform=ax.transAxes, fontsize=14)
    ax.legend(frameon=False, ncol=2, loc="upper right")
    save(fig, f"avalanche_pdf_integerN_L_{L}_steps_32768")


def fig_collapse(P, Ds, tau, b, s0):
    """The decisive picture: P(s) s^tau vs s/s0(N) collapses every N onto the one
    master curve exp[-(s/s0)^b]. Only the scale s0 differs between N."""
    fig, ax = plt.subplots(figsize=(8.2, 6.4))
    cols = plt.cm.viridis(np.linspace(0, 0.9, len(NS)))
    for N, c in zip(NS, cols):
        v, m = P[N]
        keep = (m > 0) & (v >= SMIN)
        s = v[keep].astype(float)
        pv = m[keep].astype(float) / m.sum()
        x = s / s0[N]
        y = pv * s ** tau
        ax.plot(x, y / y[0] * np.exp(-x[0] ** b), ".", ms=3.5, color=c,
                alpha=0.5, label=rf"$N={N}$")
    xx = np.logspace(np.log10(0.3), np.log10(30), 300)
    ax.plot(xx, np.exp(-xx ** b), "k-", lw=2.4, label=rf"$\exp[-(s/s_0)^{{{b:.2f}}}]$")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_ylim(1e-9, 5)
    ax.set_xlabel(r"$s / s_0(N)$")
    ax.set_ylabel(r"$P(s)\, s^{\tau}\ $ (rescaled)")
    ax.text(0.04, 0.06, rf"$\tau={tau:.2f}$, $b={b:.2f}$ fixed for all $N$"
            "\n" rf"only $s_0(N)$ varies: ${s0[NS[0]]:.0f}\to{s0[NS[-1]]:.1f}$",
            transform=ax.transAxes, fontsize=13)
    ax.legend(frameon=False, ncol=2, loc="upper right", handletextpad=0.3,
              columnspacing=0.8)
    save(fig, f"tail_collapse_fixedshape_L_{L}_steps_32768")


if __name__ == "__main__":
    P = load_pooled()
    Ds = {N: prep(*P[N], SMIN) for N in NS}
    tau, b, s0 = fit(Ds)

    print(f"\nL={L}, integer N, tail s>={SMIN}.  ONE fixed shape, only s0(N) varies:")
    print(f"  P(s) ~ s^-{tau:.3f} exp[-(s/s0)^{b:.3f}]\n")
    print("  N        s0     <s2>/<s>   worst-factor   n_tail")
    fits = {"tau": float(tau), "b": float(b), "s0": {}}
    for N in NS:
        v, m = P[N]
        vf = v.astype(float)
        cut = (vf ** 2 * m).sum() / (vf * m).sum()
        wf = worst_factor(tau, b, s0[N], Ds[N])
        print("%3d   %7.2f    %7.2f      %5.2f       %.2e"
              % (N, s0[N], cut, wf, Ds[N]["n"]))
        fits["s0"][str(N)] = float(s0[N])
    json.dump(fits, open(os.path.join(HERE, ".integerN_fits.json"), "w"))

    fig_pdf(P, Ds, tau, b, s0)
    fig_collapse(P, Ds, tau, b, s0)