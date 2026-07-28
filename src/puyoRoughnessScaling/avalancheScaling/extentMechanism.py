"""WHY P(s) is s^-tau exp[-(s/s0)^b] with b~0.6, and why only s0 moves with N.

The mass distribution is not fundamental -- it is inherited from the cascade's
spatial EXTENT w (the number of distinct columns it eliminates in), through two
facts that this script measures:

  (1) GEOMETRY.  Mass grows with extent as a power, <s|w> ~ w^d, with d ~ 1.72
      independent of N.  (An extent-w cascade digs a depth ~ w^{d-1}, so it is a
      compact 2D-ish region of the surface, not a thin sliver.)

  (2) SPREADING.  The extent is a TRUNCATED POWER LAW, P(w) ~ w^{-a} exp(-w/w*)
      (a ~ 3.0->2.2, w* ~ 4.9->1.8).  This is where the power law is data-forced,
      not assumed: `test_pw_form` shows a plain exponential is off by 1e4-1e5, and
      a stretched exponential or a lognormal in w is rejected too (worst factor
      2-5 vs 1.1).  A power-law extent means the lateral spread is (weakly
      near-critical) BRANCHING -- an eliminated column topples its neighbours, a
      scale-free process -- cut off at a finite width w* set by how often a
      neighbour matches (~1/N).  Only w* changes with N.

Change variables s = C w^d in P(w):

      P(s) ~ s^{-tau} exp[-(s/s0)^b],   b = 1/d,   s0 = C w*^d,
                                        tau = (a + d - 1)/d.

So the STRETCH EXPONENT is the reciprocal of the mass-extent exponent, b = 1/d ~
0.58 (fitted b = 0.62), N-independent because the geometry d is; the whole
N-dependence is the cutoff s0 = C w*^d.

HONESTY, both ways (the probe is the tail worst-factor, not raw lnL, which on
~1e7-1e9 counts calls everything significant):
  * P(s) does NOT by itself require the power-law prefactor.  `test_ps_prefactor`
    shows a pure stretched exp exp[-(s/s0)^b] fits the P(s) tail about as well
    (worst factor 1.1-1.3 vs 1.0-1.1); the power law buys only dlnL ~ 30-270,
    the same order this repo dismisses when freezing exponents.  The prefactor's
    justification is P(w) (above), not P(s).  tau is likewise inherited from the
    extent power a and is not pinned by P(s) alone.
  * The b = 1/d composition law was earlier thought REFUTED, on the grounds that
    a per-N free fit made b appear to rise while d stayed flat.  That b(N) rise
    was the fitting artifact (integerN.py); with b correctly fixed, flat 1/d =
    flat b holds -- consistent with, though not proven by, P(s).
"""
import os

import numpy as np
import matplotlib.pyplot as plt

import glob
from scipy.optimize import minimize
from common import slope_resolved, PLOTS, SLOPE, nfmt, STEPS
from integerN import load_pooled

plt.rcParams.update({
    "font.size": 16, "axes.labelsize": 21, "xtick.labelsize": 14,
    "ytick.labelsize": 14, "legend.fontsize": 12, "lines.linewidth": 2.0,
    "figure.autolayout": True,
})

L = 1024
NS = [6, 7, 8, 9, 10, 11, 12, 14, 16, 20]
B_FIT, TAU_FIT = 0.62, 2.03      # the fitted P(s) shape (integerN.py)


def measure():
    """d from <s|w>~w^d and (a, w*) from P(w)~w^-a e^-w/w*, per N."""
    out = {}
    for N in NS:
        _, _, _, wk, Pw, sw = slope_resolved(L, N)
        if wk is None:
            continue
        good = (sw > 0) & (wk >= 3) & (wk <= 60)
        d = np.polyfit(np.log(wk[good]), np.log(sw[good]), 1)[0]
        tail = (Pw > 0) & (wk >= 4)
        A = np.vstack([np.ones(tail.sum()), -np.log(wk[tail]), -wk[tail]]).T
        c, a, inv = np.linalg.lstsq(A, np.log(Pw[tail]), rcond=None)[0]
        out[N] = dict(d=d, a=a, wstar=1.0 / inv if inv > 0 else np.inf,
                      wk=wk, Pw=Pw, sw=sw)
    return out


# --------------------------------------------------------------------------
# Honesty checks: is the power law actually REQUIRED, at either level?  The
# probe is the tail worst-factor (max multiplicative model/data gap on log bins
# with >=200 counts) -- not raw lnL, which on ~1e7-1e9 counts calls everything
# "significant".  A component is needed only if dropping it visibly hurts the
# worst-factor, the same standard integerN.py uses to freeze the exponents.
# --------------------------------------------------------------------------
def _norm_logp(lw):
    mx = lw.max()
    return lw - mx - np.log(np.exp(lw - mx).sum())


def _mle(models, xmin, x_obs, c_obs, cap):
    """Fit each model (name -> (logp_fn, [starts])) by discrete MLE; return
    {name: (lnL, params, worst_factor)}."""
    x_all = np.arange(xmin, cap + 1, dtype=float)
    lx = np.log(x_all)
    idx = (x_obs - xmin).astype(int)
    n = c_obs.sum()
    edges = np.unique(np.round(np.logspace(np.log10(xmin),
                      np.log10(x_obs.max() + 1), 22)).astype(int))
    out = {}
    for name, (fn, starts) in models.items():
        best = None
        for x0 in starts:
            def nll(th, fn=fn):
                lw = fn(th, x_all, lx)
                if not np.all(np.isfinite(lw)):
                    return 1e18
                return -(c_obs * _norm_logp(lw)[idx]).sum()
            r = minimize(nll, x0, method="Nelder-Mead",
                         options=dict(xatol=1e-8, fatol=1e-6, maxiter=20000))
            if best is None or r.fun < best.fun:
                best = r
        p = np.exp(_norm_logp(fn(best.x, x_all, lx)))
        wf = 1.0
        for a, b in zip(edges[:-1], edges[1:]):
            sel = (x_obs >= a) & (x_obs < b)
            cnt = c_obs[sel].sum()
            if cnt < 200:
                continue
            ob, ex = cnt / n, p[(x_all >= a) & (x_all < b)].sum()
            if ob > 0 and ex > 0:
                wf = max(wf, ob / ex, ex / ob)
        out[name] = (-best.fun, best.x, wf)
    return out


def test_ps_prefactor(P, smin=50, cap=6000):
    """Is the power-law prefactor of P(s) required?  Pure stretched vs pl*stretched."""
    print("\n=== P(s): pure stretched vs power-law x stretched (smin=%d) ===" % smin)
    print(" N | pure: b   (wf) | pl*str: tau  b    (wf) | dlnL(pl-pure)")
    models = {
        "pure": (lambda th, s, ls: -(s / np.exp(th[0])) ** np.exp(th[1]),
                 [[3.0, -0.5], [2.0, -0.3]]),
        "pl":   (lambda th, s, ls: -th[0] * ls - (s / np.exp(th[1])) ** np.exp(th[2]),
                 [[2.0, 3.0, -0.5], [1.0, 2.5, -0.3], [0.4, 2.0, -0.6]]),
    }
    for N in NS:
        v, m = P[N]
        keep = (m > 0) & (v >= smin)
        r = _mle(models, smin, v[keep].astype(float), m[keep].astype(float), cap)
        (lp, xp, wp), (ll, xl, wl) = r["pure"], r["pl"]
        print(" %2d | %.2f (%.2f) | %4.2f %.2f (%.2f) | %+7.1f"
              % (N, np.exp(xp[1]), wp, xl[0], np.exp(xl[2]), wl, ll - lp))


def _pool_w(N, cap=3000):
    acc = {}
    for f in sorted(glob.glob(f"{SLOPE}/L_{L}_N_{nfmt(N)}_steps_{STEPS}_sim_*.tsv")):
        for line in open(f).read().partition("# extent")[2].strip().split("\n"):
            if not line or line.startswith("w\t"):
                continue
            p = line.split()
            acc[int(p[0])] = acc.get(int(p[0]), 0) + int(p[1])
    w = np.array(sorted(k for k in acc if k <= cap))
    return w, np.array([acc[k] for k in w], float)


def test_pw_form(wmin=4, cap=3000):
    """Is P(w) a POWER LAW x exp, or does another sub-exponential fit as well?"""
    print("\n=== P(w): pow*exp vs stretched vs lognormal (wmin=%d) ===" % wmin)
    print(" N | pow*exp: a  w*  (wf) | stretch (wf,dlnL) | lognorm (wf,dlnL)")
    models = {
        "pow": (lambda th, w, lw: -th[0] * lw - w / np.exp(th[1]), [[2.5, 1.5], [3, 1]]),
        "str": (lambda th, w, lw: -(w / np.exp(th[0])) ** np.exp(th[1]), [[1.0, -0.3], [0.5, 0.0]]),
        "logn": (lambda th, w, lw: -lw - (lw - th[0]) ** 2 / (2 * np.exp(th[1]) ** 2), [[1.5, 0.0]]),
    }
    for N in NS:
        w, c = _pool_w(N, cap)
        keep = (c > 0) & (w >= wmin)
        r = _mle(models, wmin, w[keep].astype(float), c[keep], cap)
        (lpo, xpo, wpo), (_, _, wst), (_, _, wln) = r["pow"], r["str"], r["logn"]
        print(" %2d | %4.2f %4.1f (%.2f) | (%.2f, %+8.1f) | (%.2f, %+8.1f)"
              % (N, xpo[0], np.exp(xpo[1]), wpo,
                 wst, r["str"][0] - lpo, wln, r["logn"][0] - lpo))


def figure(M):
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.6))
    cols = plt.cm.viridis(np.linspace(0, 0.9, len(NS)))

    # (a) <s|w> ~ w^d : one slope for all N -> b = 1/d is N-independent
    ax = axes[0]
    for N, c in zip(NS, cols):
        wk, sw = M[N]["wk"], M[N]["sw"]
        sel = (sw > 0) & (wk >= 2)
        ax.plot(wk[sel], sw[sel], "o", ms=4, color=c, alpha=0.7, label=rf"$N={N}$")
    dbar = np.mean([M[N]["d"] for N in NS])
    ww = np.array([2, 80], float)
    ax.plot(ww, sw[sel][0] / (wk[sel][0]) ** dbar * ww ** dbar * 1.0, "k--", lw=2,
            label=rf"$w^{{{dbar:.2f}}}$")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"cascade extent $w$ (columns)")
    ax.set_ylabel(r"$\langle s \,|\, w\rangle$")
    ax.text(0.05, 0.9, rf"$d={dbar:.2f}\Rightarrow b=1/d={1/dbar:.2f}$"
            "\n" rf"(fitted $b={B_FIT}$)", transform=ax.transAxes, fontsize=14)
    ax.legend(frameon=False, ncol=2, loc="lower right", handletextpad=0.3,
              columnspacing=0.8)

    # (b) P(w): a truncated power law w^-a e^-w/w* -- a straight power-law
    # stretch in log-log, then the w*-cutoff. (A pure exponential, straight in
    # semi-log, is excluded by test_pw_form: off by 1e4-1e5.)
    ax = axes[1]
    for N, c in zip(NS, cols):
        wk, Pw = M[N]["wk"], M[N]["Pw"]
        sel = Pw > 0
        ax.plot(wk[sel], Pw[sel], "o", ms=4, color=c, alpha=0.7,
                label=rf"$N={N}$, $w^*={M[N]['wstar']:.1f}$")
    abar = np.mean([M[N]["a"] for N in NS])
    ww = np.array([2.0, 12.0])
    ax.plot(ww, 0.5 * ww ** (-abar), "k--", lw=2, label=rf"$w^{{-{abar:.1f}}}$")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(1, 60)
    ax.set_ylim(1e-8, 1)
    ax.set_xlabel(r"cascade extent $w$ (columns)")
    ax.set_ylabel(r"$P(w)$")
    ax.text(0.97, 0.95, r"$P(w)\sim w^{-a}e^{-w/w^*}$"
            "\n" r"(power law forced; not $e^{-w/w^*}$)",
            transform=ax.transAxes, fontsize=13, ha="right", va="top")
    ax.legend(frameon=False, ncol=2, loc="lower left", fontsize=9.5,
              handletextpad=0.3, columnspacing=0.8)

    d = f"{PLOTS}/integerN"; os.makedirs(d, exist_ok=True)
    fig.savefig(f"{d}/extent_mechanism_L_{L}_steps_32768.svg")
    fig.savefig(f"{d}/extent_mechanism_L_{L}_steps_32768.png", dpi=150)
    plt.close(fig)
    print(f"saved integerN/extent_mechanism_L_{L}_steps_32768")


if __name__ == "__main__":
    M = measure()
    print(" N      d      1/d      a      w*     tau_pred=(a+d-1)/d")
    for N in NS:
        m = M[N]
        print("%3d   %.3f   %.3f   %.2f   %5.1f    %.2f"
              % (N, m["d"], 1 / m["d"], m["a"], m["wstar"], (m["a"] + m["d"] - 1) / m["d"]))
    print(f"\nfitted P(s): tau={TAU_FIT}, b={B_FIT};  mean 1/d = "
          f"{np.mean([1/M[N]['d'] for N in NS]):.3f}")

    # the honesty checks the readme leans on
    test_ps_prefactor(load_pooled())
    test_pw_form()
    figure(M)