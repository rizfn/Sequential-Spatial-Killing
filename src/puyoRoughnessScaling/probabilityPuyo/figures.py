"""Figures for 'what variable linearises v(N)?'.

  staircase     v(N) is a staircase: inside each interval it jumps at f=0+ and
                then flattens.  This is the thing to be explained
  collapse      the answer: v vs sum(p.rho) is straight across ALL N, and
                reduces to 1/N at integer N.  Entropy (1/N, sum p^2) does not
  entropy_scan  the negative: NO Renyi/Hill order q linearises v
  impurity      the mechanism: the rare species self-enriches (rho ~ f^b), so
                the pile composition is not the deposition composition
"""
import os

import numpy as np
import matplotlib.pyplot as plt

import common
from common import PLOTS, hill, weights, linearity

plt.rcParams.update({
    "font.size": 17, "axes.labelsize": 21, "xtick.labelsize": 15,
    "ytick.labelsize": 15, "legend.fontsize": 13, "lines.linewidth": 2.1,
    "figure.autolayout": True,
})


def save(fig, name):
    os.makedirs(PLOTS, exist_ok=True)
    fig.savefig(f"{PLOTS}/{name}.svg")
    fig.savefig(f"{PLOTS}/{name}.png", dpi=150)
    plt.close(fig)
    print(f"saved {name}")


T = common.table()
NS = sorted(T)
V = np.array([T[N]["v"] for N in NS])
INT = [N for N in NS if abs(N - round(N)) < 1e-9]


# --------------------------------------------------------------------------
def staircase():
    fig, ax = plt.subplots(figsize=(8.5, 6))
    ax.plot(NS, V, "o-", ms=7)
    ax.plot(INT, [T[N]["v"] for N in INT], "s", ms=13, color="crimson",
            zorder=5, label="integer $N$")
    for N in INT:
        ax.axvline(N, color="grey", ls=":", lw=1.2)
    ax.set(xlabel=r"$N$", ylabel=r"$v = \mathrm{d}(M/L)/\mathrm{d}t$")
    ax.legend()
    save(fig, f"staircase_L_{common.L}_steps_{common.STEPS}")


# --------------------------------------------------------------------------
def collapse():
    """v against the naive variable and against the right one."""
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13.5, 5.8))
    for ax, key, lab in [(a1, "1/N", r"$1/N$"),
                         (a2, "sum p.rho", r"$\sum_i p_i \rho_i$")]:
        x = np.array([T[N]["cands"][key] for N in NS])
        ax.plot(x, V, "o", ms=8, label="all $N$")
        xi = np.array([T[N]["cands"][key] for N in INT])
        ax.plot(xi, [T[N]["v"] for N in INT], "s", ms=13, color="crimson",
                zorder=5, label="integer $N$")
        rms, _, p = linearity(x, V)
        xx = np.linspace(x.min(), x.max(), 10)
        ax.plot(xx, np.polyval(p, xx), "k--", lw=2,
                label=f"line fit, RMS={rms:.4f}")
        ax.set(xlabel=lab, ylabel=r"$v$")
        ax.legend(loc="upper right")
    save(fig, f"collapse_L_{common.L}_steps_{common.STEPS}")


# --------------------------------------------------------------------------
def entropy_scan():
    """Is it information entropy?  Scan every Renyi order.  No: the best order
    is an unprincipled q~0.7 and is still far worse than sum(p.rho).
    """
    qs = np.arange(0.2, 3.01, 0.05)
    rms = []
    for q in qs:
        x = np.array([1 / hill(weights(N), q) for N in NS])
        rms.append(linearity(x, V)[0])
    rms = np.array(rms)

    noise = float(np.mean([T[N]["err"] for N in NS]))
    r_int = linearity(np.array([1 / N for N in INT]),
                      np.array([T[N]["v"] for N in INT]))[0]
    r_pr = linearity(np.array([T[N]["cands"]["sum p.rho"] for N in NS]), V)[0]

    fig, ax = plt.subplots(figsize=(8.5, 6))
    ax.plot(qs, rms, "-", lw=2.5, label=r"$1/{}^{q}\!D$  (Hill order $q$)")
    ax.axhline(r_pr, color="crimson", ls="--", lw=2.4,
               label=r"$\sum_i p_i\rho_i$")
    ax.axhline(r_int, color="grey", ls=":", lw=2.4,
               label=r"$1/N$, integer $N$ only (benchmark)")
    ax.axhline(noise, color="k", ls="-", lw=1.5, label="noise floor")
    for q, lab in [(1.0, "Shannon"), (2.0, "collision")]:
        ax.axvline(q, color="grey", lw=1, ls=":")
        ax.text(q + 0.03, rms.max() * 0.92, lab, fontsize=13, rotation=90)
    ax.set(xlabel=r"Rényi order $q$", ylabel="RMS residual of linear fit",
           yscale="log")
    ax.legend(fontsize=12, loc="lower right")
    save(fig, "entropy_scan")


# --------------------------------------------------------------------------
def impurity():
    """Why every function of p alone must fail: the pile is not the deposit."""
    fig, ax = plt.subplots(figsize=(8.5, 6))
    for lo, c in zip([6, 7, 8, 9], plt.cm.viridis(np.linspace(0, .85, 4))):
        fs = [N - lo for N in NS if lo < N <= lo + 1]
        en = [T[N]["rho"][-1] / ((N - lo) / N) for N in NS if lo < N <= lo + 1]
        b = np.polyfit(np.log(fs), np.log([T[N]["rho"][-1] for N in NS
                                           if lo < N <= lo + 1]), 1)[0]
        ax.plot(fs, en, "o-", ms=8, color=c,
                label=rf"$N\in[{lo},{lo+1}]$:  $\rho\sim f^{{{b:.2f}}}$")
    ax.axhline(1, color="k", ls="--", lw=2)
    ax.set(xlabel=r"fractional weight $f = N - \lfloor N \rfloor$",
           ylabel="enrichment  $\\rho_{\\rm imp} / (f/N)$")
    ax.legend(fontsize=12)
    save(fig, "impurity_enrichment")


if __name__ == "__main__":
    staircase()
    collapse()
    entropy_scan()
    impurity()
