"""Figures for the growth-arrest transition at N_c.

  unbinding       THE result: the arrest height h* diverges as N -> N_c and stays
                  proportional to the interface width W, so the transition is an
                  interface UNBINDING transition and h* is the diverging length
  finite_size     h*(L) at N=5.075: converges for L >= 4096 but GROWS at small L.
                  This is why a single-L estimate of N_c is biased, and it is the
                  evidence that the arrest survives the thermodynamic limit
  psi_fit         h* ~ (N_c-N)^-psi with N_c BRACKETED (not fitted) -> psi=0.58(3)
  scale_audit     the honest limit: h* reaches only ~5 UV cutoffs l*(N), so the
                  exponents are effective, not asymptotic
"""
import os

import numpy as np
import matplotlib.pyplot as plt

import common
from common import PLOTS

plt.rcParams.update({
    "font.size": 17, "axes.labelsize": 21, "xtick.labelsize": 15,
    "ytick.labelsize": 15, "legend.fontsize": 13, "lines.linewidth": 2.1,
    "figure.autolayout": True,
})

# N_c is BRACKETED, not fitted: at L=4096 every N <= 5.075 arrests and every
# N >= 5.078 grows linearly.  Pinning it this way is what makes psi identifiable.
N_C = 5.0765
N_C_ERR = 0.0015
L_MAIN, STEPS_MAIN = 4096, 2000000
L_UV = 20.0     # UV cutoff l*(N) near N_c, from roughness-not-kpz (l_G ~ 19.8 at N=6)

N_ALL = [5.00, 5.02, 5.04, 5.05, 5.06, 5.065, 5.07, 5.072, 5.075, 5.078, 5.080, 5.085, 5.090]
N_BOUND = [5.00, 5.02, 5.04, 5.05, 5.06, 5.065, 5.070, 5.072, 5.075]

# h* at N=5.075 must come from L=8192: the finite_size panel shows L=4096 is
# still drifting there, because xi_par has grown to a sizeable fraction of L.
HSTAR_L = {N: L_MAIN for N in N_BOUND}
HSTAR_L[5.075] = 8192


def save(fig, name):
    os.makedirs(PLOTS, exist_ok=True)
    fig.savefig(f"{PLOTS}/{name}.svg")
    fig.savefig(f"{PLOTS}/{name}.png", dpi=150)
    plt.close(fig)
    print(f"saved {name}")


def plateau(N, L, steps=STEPS_MAIN, key="h"):
    """Late-time plateau of <h> or W, averaged over the last stretch of samples.
    Returns None if the run is absent or was corrupted by a box-ceiling hit --
    a ceiling caps the pile and fakes an arrest, so those must never be used.
    """
    if common.ceiling_hits(L, N, steps) > 0:
        return None
    d = common.load(L, N, steps)
    return None if d is None else float(d[key][-40:].mean())


# --------------------------------------------------------------------------
def unbinding():
    Ns = [N for N in N_BOUND if plateau(N, HSTAR_L[N]) is not None]
    h = np.array([plateau(N, HSTAR_L[N]) for N in Ns])
    w = np.array([plateau(N, HSTAR_L[N], key="w") for N in Ns])
    d = np.array([N_C - N for N in Ns])

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13.5, 5.8))
    a1.plot(d, h, "o-", ms=9, label=r"$h^*$  (binding length)")
    a1.plot(d, w, "s-", ms=9, label=r"$W^*$  (interface width)")
    a1.axhline(L_UV, color="k", ls=":", lw=2)
    a1.text(0.055, L_UV * 1.12, r"UV cutoff $\ell^*$", fontsize=13)
    a1.set(xscale="log", yscale="log", xlabel=r"$N_c - N$", ylabel="length")
    a1.invert_xaxis()
    a1.legend(loc="lower left")

    a2.plot(Ns, h / w, "o-", ms=9, color="crimson")
    a2.axvline(N_C, color="k", ls="--")
    a2.text(N_C - 0.001, 0.6, r"$N_c$", fontsize=18, ha="right")
    a2.set(xlabel=r"$N$", ylabel=r"$h^* / W^*$", ylim=(0, 4.5))
    save(fig, f"unbinding_steps_{STEPS_MAIN}")


# --------------------------------------------------------------------------
def finite_size():
    """At N=5.075 the pile GROWS at small L but ARRESTS at large L.  So the
    arrest is a thermodynamic-limit property, and any N_c read off a single
    small L is biased.  The old N_c=5.075(10) came from L=512.
    """
    Ls, hs = [], []
    for L in [512, 1024, 2048, 4096, 8192]:
        v = plateau(5.075, L)
        if v is not None:
            Ls.append(L); hs.append(v)

    fig, ax = plt.subplots(figsize=(7.8, 6))
    ax.plot(Ls, hs, "o-", ms=12)
    ax.axhline(97.2, color="crimson", ls="--", lw=2,
               label=r"$h^* \simeq 97$ (converged)")
    ax.set(xscale="log", yscale="log", xlabel=r"$L$",
           ylabel=r"$\langle h \rangle$ at $t=2\times10^6$,  $N=5.075$")
    ax.legend()
    save(fig, "finite_size_N_5.075")


# --------------------------------------------------------------------------
def psi_fit():
    """h* ~ (N_c-N)^-psi.  Only points with h* > l* are fitted: below the UV
    cutoff there is no continuum and hence nothing to scale.
    """
    Ns = [N for N in N_BOUND
          if plateau(N, HSTAR_L[N]) is not None and plateau(N, HSTAR_L[N]) > L_UV]
    h = np.array([plateau(N, HSTAR_L[N]) for N in Ns])

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13.5, 5.8))
    d = np.array([N_C - N for N in Ns])
    p = np.polyfit(np.log(d), np.log(h), 1)
    a1.plot(d, h, "o", ms=11)
    xx = np.logspace(np.log10(d.min()), np.log10(d.max()), 20)
    a1.plot(xx, np.exp(np.polyval(p, np.log(xx))), "-",
            label=rf"$\psi = {-p[0]:.2f}$,  $N_c={N_C}$")
    a1.set(xscale="log", yscale="log", xlabel=r"$N_c - N$", ylabel=r"$h^*$")
    a1.invert_xaxis()
    a1.legend()

    # psi against the assumed N_c.  N_c must stay above the largest bound N,
    # or N_c-N goes negative for that point and the fit is undefined.
    ncs = np.linspace(max(Ns) + 0.0003, N_C + 2 * N_C_ERR, 25)
    psis = [-np.polyfit(np.log([nc - N for N in Ns]), np.log(h), 1)[0] for nc in ncs]
    a2.plot(ncs, psis, "-", lw=2.5, color="crimson")
    a2.axvspan(N_C - N_C_ERR, N_C + N_C_ERR, color="grey", alpha=0.25,
               label="bracketed $N_c$")
    a2.axhline(1 / 3, color="k", ls="--", lw=2)
    a2.text(ncs[0], 0.35, r"bEW: $\psi=1/3$", fontsize=13)
    a2.set(xlabel=r"assumed $N_c$", ylabel=r"fitted $\psi$")
    a2.legend(loc="lower right")
    save(fig, "psi_fit")


# --------------------------------------------------------------------------
def scale_audit():
    Ns = [N for N in N_BOUND if plateau(N, HSTAR_L[N]) is not None]
    h = np.array([plateau(N, HSTAR_L[N]) for N in Ns])

    fig, ax = plt.subplots(figsize=(7.8, 6))
    ax.plot([N_C - N for N in Ns], h / L_UV, "o-", ms=11)
    ax.axhline(1, color="k", ls="--", lw=2)
    ax.axhspan(0.1, 1.0, color="crimson", alpha=0.15)
    ax.text(0.06, 1.15, r"$h^*=\ell^*$: no continuum below", fontsize=13)
    ax.text(0.06, 0.45, "no scaling regime", fontsize=13, color="darkred")
    ax.set(xscale="log", yscale="log", xlabel=r"$N_c - N$",
           ylabel=r"$h^*/\ell^*$   (binding length in UV cutoffs)")
    ax.invert_xaxis()
    ax.set_ylim(0.2, 12)
    save(fig, "scale_audit")


if __name__ == "__main__":
    unbinding()
    finite_size()
    psi_fit()
    scale_audit()
