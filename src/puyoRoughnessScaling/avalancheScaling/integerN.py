"""The shape of P(s) at INTEGER N, and how it changes with N.

Fractional N is a different system -- the partial species is a self-poisoning
impurity that fragments clusters (readme's sawtooth section) -- so this family
is integer N only: N = 6,7,8,9,10,11,12,14,16,20 at L=1024, 1024 sims each
(up to 7e9 active avalanches at N=6).

METHOD.  P(s) is fitted by discrete MAXIMUM LIKELIHOOD on the raw integer
histogram, not by least squares on a binned PDF:

    ln Lik = sum_s count[s] * ln p(s),    p normalized on s = smin..CAP

This is what makes the 3-parameter form P(s) ~ s^-tau exp[-(s/s0)^b]
IDENTIFIABLE.  The readme's least-squares fit of A s^-tau exp[-(s/s0)^b] ran
away to tau=-9.57, b=0.081, because a FREE amplitude A absorbs whatever tau and
b do, leaving a flat tau<->b ridge.  Under MLE the amplitude is not free -- it
IS the normalization -- so the ridge is lifted and the fit is sharp.

RESULT.
  * At N=6 the tail is a power law times a STRETCHED exponential,
        P(s) ~ s^-2.0 exp[-(s/17)^0.62],   s >~ 20,
    good to ~25% over ~8 decades.  A pure power law is off by 97x; b=1 (plain
    exponential cutoff) is excluded by dlnL~300; forcing N=6 to N=20's b=1.76
    costs dlnL=21000 and a 9.7x misfit.  So the stretch at low N is decisive.
  * Both tau AND b RISE monotonically with N.  b passes 1 near N~11: below it
    the cutoff is STRETCHED (b<1, a fat sub-exponential tail), above it
    COMPRESSED (b>1, a Gaussian-like sharp cutoff).  As species multiply, matches
    get rarer, cascades die faster, and the tail is cut off ever more sharply.
  * This is NOT the readme's b = 1/d composition law: the extent exponent d is
    flat in N (1.77 -> 1.85), so 1/d ~ 0.55 predicts NO N dependence.  It lands
    near the measured b only at N=6 by coincidence; the whole trend is unexplained
    by that mechanism.

CAVEAT.  b is well-pinned at low N (s0 ~ 17 sits at the bottom of the tail, so
the cutoff is heavily sampled) but its systematic error grows with N (at N=20,
s0 ~ 105 sits near s_max, so the cutoff region is thin).  The MONOTONIC RISE and
the crossing of b=1 survive every smin choice; the precise b at N=20 does not.

Below s ~ 18 a strong even/odd parity oscillation (even masses favoured -- the
elementary cluster is a pair) means no smooth form applies, so fits start at
smin=30.
"""
import json
import os

import numpy as np
import matplotlib.pyplot as plt

import mle
import taueff
from common import PLOTS, logbin, pooled

plt.rcParams.update({
    "font.size": 17, "axes.labelsize": 21, "xtick.labelsize": 15,
    "ytick.labelsize": 15, "legend.fontsize": 13, "lines.linewidth": 2.1,
    "figure.autolayout": True,
})

L = 1024
NS = [6, 7, 8, 9, 10, 11, 12, 14, 16, 20]
SMIN = 30
SMINS = [20, 30, 40]      # smin=60 is unreliable at high N (tail too short)
SUB = "integerN"
CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)), ".integerN_cache.npz")
mle.CAP = 6000            # s_max <= 1394; the model tail is negligible far below


def save(fig, name):
    d = f"{PLOTS}/{SUB}"
    os.makedirs(d, exist_ok=True)
    fig.savefig(f"{d}/{name}.svg")
    fig.savefig(f"{d}/{name}.png", dpi=150)
    plt.close(fig)
    print(f"saved {SUB}/{name}")


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


def counts_dict(v, m):
    return {int(s): int(c) for s, c in zip(v, m) if c > 0}


def fit_all(P):
    """MLE fit at each N, with the smin spread carried as the systematic error."""
    out = {}
    for N in NS:
        v, m = P[N]
        cb = counts_dict(v, m)
        th, lnL, n = mle.fit(cb, SMIN, mle.m_plstretch, [2.5, 3.0, -0.3])
        spread = [mle.fit(cb, s, mle.m_plstretch, [2.5, 3.0, -0.3])[0] for s in SMINS]
        tau_s = [p[0] for p in spread]
        b_s = [np.exp(p[2]) for p in spread]
        vf = v.astype(float)
        out[N] = dict(tau=th[0], b=np.exp(th[2]), s0=np.exp(th[1]), th=th, n=n,
                      tau_err=(max(tau_s) - min(tau_s)) / 2,
                      b_err=(max(b_s) - min(b_s)) / 2,
                      cut=(vf ** 2 * m).sum() / (vf * m).sum(),
                      smax=int(v[m > 0].max()))
    return out


def fig_pdf(P, F):
    """P(s) at every integer N, with the MLE curves over the fitted range."""
    fig, ax = plt.subplots(figsize=(8, 6.5))
    show = [6, 7, 8, 10, 12, 16, 20]
    cols = plt.cm.viridis(np.linspace(0, 0.9, len(show)))
    for N, c in zip(show, cols):
        v, m = P[N]
        cb = counts_dict(v, m)
        x, y = logbin(v, m, nbins=45)
        ax.plot(x, y, "o", ms=4.5, color=c, alpha=0.75,
                label=rf"$N={N}$, $b={F[N]['b']:.2f}$")
        th = F[N]["th"]
        s = np.arange(SMIN, mle.CAP + 1, dtype=float)
        p = np.exp(mle._norm_logp(mle.m_plstretch(th, s), s))
        frac = sum(cc for ss, cc in cb.items() if ss >= SMIN) / m.sum()
        sel = s <= v[m > 0].max()
        ax.plot(s[sel], p[sel] * frac, "-", color=c, lw=1.6)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"avalanche mass $s$")
    ax.set_ylabel(r"$P(s)$")
    ax.set_ylim(1e-10, 3)
    ax.legend(frameon=False, ncol=2)
    save(fig, f"avalanche_pdf_integerN_L_{L}_steps_32768")


def fig_exponents(F):
    """b(N) and tau(N).  b crosses 1 near N~11: stretched below, compressed above."""
    fig, axes = plt.subplots(2, 1, figsize=(7.6, 8.6), sharex=True)
    ns = np.array(NS, float)

    ax = axes[0]
    b = np.array([F[N]["b"] for N in NS])
    be = np.array([F[N]["b_err"] for N in NS])
    ax.errorbar(ns, b, yerr=be, fmt="o-", ms=8, capsize=4, color="#1f4e79")
    ax.axhline(1.0, color="k", ls=":", lw=1.6)
    ax.axhline(0.55, color="gray", ls="-.", lw=1.4)
    ax.set_ylim(0.4, 2.15)
    ax.text(6.3, 1.08, "plain exponential", fontsize=13)
    ax.text(19.7, 0.62, r"readme $1/d\approx0.55$", fontsize=12, color="gray", ha="right")
    ax.text(9.5, 1.85, "stretched\n(fat tail)", fontsize=13, style="italic", ha="center")
    ax.text(17.5, 1.05, "compressed\n(sharp cutoff)", fontsize=13, style="italic",
            ha="center")
    ax.set_ylabel(r"stretch exponent $b$")

    ax = axes[1]
    t = np.array([F[N]["tau"] for N in NS])
    te = np.array([F[N]["tau_err"] for N in NS])
    ax.errorbar(ns, t, yerr=te, fmt="o-", ms=8, capsize=4, color="#c1440e")
    ax.axhline(2.0, color="gray", ls="-.", lw=1.4)
    ax.set_ylabel(r"apparent exponent $\tau$")
    ax.set_xlabel(r"number of species $N$")
    save(fig, f"stretch_exponent_vs_N_L_{L}_steps_32768")


def fig_localslope(P, F):
    """Model-independent view: the local log-slope tau_eff(s) = -dlnP/dlns.
    For P ~ s^-tau exp[-(s/s0)^b] it rises from a tau plateau; the rise steepens
    with N.  This just displays the data's slope -- no fit is imposed."""
    fig, ax = plt.subplots(figsize=(8, 6.2))
    show = [6, 8, 10, 12, 16, 20]
    cols = plt.cm.viridis(np.linspace(0, 0.9, len(show)))
    for N, c in zip(show, cols):
        v, m = P[N]
        s, te, nn = taueff.tau_eff(v, m, nbins=30)
        sel = (nn > 800) & (s > 18)      # below ~18: even/odd parity noise, no smooth slope
        ax.plot(s[sel], te[sel], "o-", ms=5, color=c, label=rf"$N={N}$")
    ax.axhline(2.0, color="gray", ls=":", lw=1.3)
    ax.set_xscale("log")
    ax.set_ylim(1.5, 8.5)
    ax.set_xlabel(r"avalanche mass $s$")
    ax.set_ylabel(r"local slope $\tau_{\rm eff}(s)$")
    ax.legend(frameon=False, ncol=2, loc="upper left")
    save(fig, f"local_slope_integerN_L_{L}_steps_32768")


if __name__ == "__main__":
    P = load_pooled()
    F = fit_all(P)
    print()
    print(f"L={L}, 1024 sims/N, tail s>={SMIN}.   P(s) ~ s^-tau exp[-(s/s0)^b]")
    print("  N       tau            b            s0     <s2>/<s>   n_tail")
    for N in NS:
        f = F[N]
        print("%3d   %.3f(%.0f)     %.3f(%.0f)     %6.2f    %6.3f    %.2e"
              % (N, f["tau"], f["tau_err"] * 1000, f["b"], f["b_err"] * 1000,
                 f["s0"], f["cut"], f["n"]))
    json.dump({str(N): {k: float(v) for k, v in F[N].items() if k != "th"}
               for N in NS},
              open(os.path.join(os.path.dirname(CACHE), ".integerN_fits.json"), "w"))
    fig_pdf(P, F)
    fig_exponents(F)
    fig_localslope(P, F)
