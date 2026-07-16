"""Regenerate all roughness-scaling figures into plots/.

  beta_eff_vs_t   effective growth exponent d ln W / d ln t   (L=64, family of N)
  slope_collapse  steady-state slope PDF, P(m) and collapse P(m) vs m/N
  roughness_vs_L  W_sat ~ L^alpha per N (fit over saturated L)
  alpha_vs_N      the roughness exponent alpha as a function of N
  collapse_N6     Family-Vicsek collapse at N=6 vs the KPZ prediction
  family_vicsek   test of the identity  alpha = beta * z
  growth_anatomy  W(t) log-log showing the transient is not a single power law
"""
import glob
import numpy as np
import matplotlib.pyplot as plt

from common import (mean_W, stack, parse_slopes, pattern, steps_of,
                    PLOTS, FSWEEP, FIXED_L, all_N_at_fixed_L)
from exponents import measure_all, measure_fv, saturated_Ls, SAT_THRESH

plt.rcParams.update({
    "font.size": 17, "axes.labelsize": 21, "xtick.labelsize": 15,
    "ytick.labelsize": 15, "legend.fontsize": 13, "lines.linewidth": 2.1,
    "figure.autolayout": True,
})


def save(fig, name):
    fig.savefig(f"{PLOTS}/{name}.svg")
    fig.savefig(f"{PLOTS}/{name}.png", dpi=150)
    plt.close(fig)
    print("saved", name)


# --------------------------------------------------------------------------
def fig_beta_eff():
    Ns = all_N_at_fixed_L()
    cmap = plt.get_cmap("viridis", len(Ns))
    fig, ax = plt.subplots(figsize=(9, 6))
    for i, N in enumerate(Ns):
        t, W = mean_W(FIXED_L, N)
        if t is None:
            continue
        m = t > 0
        lt, lW = np.log(t[m]), np.log(W[m])
        beta = np.convolve(np.gradient(lW, lt), np.ones(9) / 9, mode="same")
        show = t[m] <= 200
        ax.plot(t[m][show], beta[show], color=cmap(i), label=f"N={N}")
    ax.axhline(0.5, ls="--", color="grey", lw=1.5)
    ax.axhline(1 / 3, ls="--", color="black", lw=1.5)
    ax.text(1.3, 0.505, "random deposition (1/2)", color="grey", va="bottom", fontsize=12)
    ax.text(1.3, 0.335, "KPZ (1/3)", color="black", va="bottom", fontsize=12)
    ax.set_xscale("log")
    ax.set_xlabel("time $t$")
    ax.set_ylabel(r"$\beta_{\rm eff}=d\ln W/d\ln t$")
    ax.set_ylim(0.15, 0.6); ax.set_xlim(1, 200)
    ax.legend(ncol=2, title=f"$L={FIXED_L}$", loc="lower right")
    save(fig, "beta_eff_vs_t")


# --------------------------------------------------------------------------
def _steady_slopes(N):
    thr = 0.7 * steps_of(N)
    out = []
    for f in sorted(glob.glob(pattern(FIXED_L, N))):
        with open(f) as fh:
            next(fh)
            for line in fh:
                p = line.split("\t")
                if len(p) < 3 or not p[2].strip() or float(p[0]) < thr:
                    continue
                out.append(parse_slopes(p[2].strip()))
    return np.concatenate(out) if out else np.array([])


def fig_slope_collapse():
    """Collapse the steady-state slope PDFs by their own Laplace scale
    lambda(N) = <|m|> (NOT by N: lambda grows super-linearly ~ N^1.5)."""
    Ns = all_N_at_fixed_L()
    cmap = plt.get_cmap("viridis", len(Ns))
    fig, (axa, axb) = plt.subplots(1, 2, figsize=(15, 6))
    for i, N in enumerate(Ns):
        m = _steady_slopes(N)
        if m.size == 0:
            continue
        lam = np.abs(m).mean()                       # MLE Laplace scale
        bins = np.arange(np.floor(m.min()) - 0.5, np.ceil(m.max()) + 1.5, 1.0)
        pdf, edges = np.histogram(m, bins=bins, density=True)
        ctr = 0.5 * (edges[:-1] + edges[1:])
        g = pdf > 0
        axa.plot(ctr[g], pdf[g], color=cmap(i), label=f"N={N}")
        axb.plot(ctr[g] / lam, lam * pdf[g], color=cmap(i), label=f"N={N}")
    xx = np.linspace(-9, 9, 200)
    axb.plot(xx, 0.5 * np.exp(-np.abs(xx)), "k--", lw=1.6,
             label=r"Laplace $\frac{1}{2} e^{-|x|}$")
    for ax in (axa, axb):
        ax.set_yscale("log")
    axa.set_xlabel(r"local slope $m=\nabla h$"); axa.set_ylabel("$P(m)$")
    axa.legend(title=f"$L={FIXED_L}$", ncol=2)
    axb.set_xlabel(r"rescaled slope $m/\lambda(N)$,  $\lambda=\langle|m|\rangle$")
    axb.set_ylabel(r"$\lambda\,P(m)$")
    axb.set_xlim(-9, 9); axb.set_ylim(1e-6, 1)
    axb.legend(fontsize=11)
    save(fig, "slope_collapse")


# --------------------------------------------------------------------------
def fig_slope_scale():
    """The Laplace scale lambda(N)=<|m|> of the slope PDF vs N: grows FASTER than
    N (so an m/N rescaling cannot collapse the family)."""
    Ns = all_N_at_fixed_L()
    lam = np.array([np.abs(_steady_slopes(N)).mean() for N in Ns])
    Na = np.array(Ns, float)
    # fit over the mid-range (N<=20, cleanest / most saturated)
    mid = Na <= 20
    lx, ly = np.log(Na[mid]), np.log(lam[mid])
    a, b = np.polyfit(lx, ly, 1)             # free exponent
    c1 = np.mean(ly - lx)                     # fixed exponent 1, best amplitude
    ssr_free = np.sum((ly - (a * lx + b)) ** 2)
    ssr_one = np.sum((ly - (lx + c1)) ** 2)
    print(f"lambda(N), N<=20:  free  a={a:.2f}  SSR={ssr_free:.4f}   |   "
          f"fixed a=1  SSR={ssr_one:.4f}  ({ssr_one/ssr_free:.0f}x worse)")
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(Na, lam, "o", ms=10, color="#2a9d8f", label=r"$\lambda=\langle|m|\rangle$")
    xx = np.linspace(Na.min(), Na.max(), 50)
    ax.plot(xx, np.exp(b) * xx ** a, "-", color="#2a9d8f",
            label=fr"free fit: $\lambda\sim N^{{{a:.2f}}}$")
    ax.plot(xx, np.exp(c1) * xx, "k--", lw=1.5, label=r"fixed fit: $\lambda\propto N$")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("number of colours $N$")
    ax.set_ylabel(r"slope scale $\lambda(N)=\langle|m|\rangle$")
    ax.legend()
    save(fig, "slope_scale")


def fig_roughness_vs_L(res):
    fig, ax = plt.subplots(figsize=(8.5, 6.5))
    Ns = list(res)
    cmap = plt.get_cmap("plasma", len(Ns) + 1)
    for j, N in enumerate(Ns):
        r = res[N]
        Ls = np.array(r["Ls"], float)
        Wsat = r["Wsat"]
        sat = np.array([r["ratios"][L] < SAT_THRESH for L in r["Ls"]])
        a = np.nanmean(r["alpha"])
        ax.plot(Ls[sat], Wsat[sat], "o-", color=cmap(j), label=fr"$N={N}$: $\alpha={a:.2f}$")
        if (~sat).any():                       # unsaturated L shown hollow, excluded from fit
            ax.plot(Ls[~sat], Wsat[~sat], "o", mfc="white", mec=cmap(j), ms=7)
    ax.plot([16, 256], [0.9 * 16 ** 0.5, 0.9 * 256 ** 0.5], "k--", lw=1.4,
            label=r"KPZ/EW $\alpha=1/2$")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("system size $L$"); ax.set_ylabel(r"saturated roughness $W_{\rm sat}$")
    ax.legend(ncol=2, fontsize=11)
    save(fig, "roughness_vs_L")


def fig_alpha_vs_N(res):
    Ns = np.array(list(res), float)
    a = np.array([np.nanmean(res[N]["alpha"]) for N in res])
    da = np.array([np.nanstd(res[N]["alpha"]) for N in res])
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.errorbar(Ns, a, yerr=da, fmt="o-", ms=9, capsize=4, color="#b5179e")
    ax.axhline(0.5, ls="--", color="grey", lw=1.4)
    ax.text(Ns.max() * 0.6, 0.51, "KPZ/EW  (1/2)", color="grey", fontsize=12)
    ax.set_xlabel("number of colours $N$")
    ax.set_ylabel(r"roughness exponent $\alpha$")
    ax.set_ylim(0, 0.55)
    save(fig, "alpha_vs_N")


# --------------------------------------------------------------------------
def _collapse_cost(data, alpha, z):
    curves = []
    for L, (t, W) in data.items():
        m = t > 0
        curves.append((np.log(t[m] / L**z), np.log(W[m] / L**alpha)))
    lo = max(c[0][0] for c in curves); hi = min(c[0][-1] for c in curves)
    if hi <= lo:
        return np.inf
    grid = np.linspace(lo, hi, 80)
    ys = np.array([np.interp(grid, c[0], c[1]) for c in curves])
    return np.mean(np.var(ys, axis=0))


def fig_collapse_N6():
    N = 6
    stacks = {L: stack(L, N) for L in FSWEEP[N]["Ls"]}
    satLs, _ = saturated_Ls(stacks)
    data = {L: (stacks[L][0], stacks[L][1].mean(0)) for L in satLs}
    best = (np.inf, None, None)
    for a in np.arange(0.18, 0.45, 0.01):
        for z in np.arange(0.5, 1.8, 0.02):
            c = _collapse_cost(data, a, z)
            if c < best[0]:
                best = (c, a, z)
    _, ab, zb = best
    print(f"N=6 best collapse alpha={ab:.2f} z={zb:.2f} "
          f"(cost {best[0]:.4f} vs KPZ {_collapse_cost(data,0.5,1.5):.4f})")
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(14, 5.5))
    cmap = plt.get_cmap("viridis", len(satLs))
    for i, L in enumerate(satLs):
        t, W = data[L]; m = t > 0
        a1.plot(t[m] / L**zb, W[m] / L**ab, color=cmap(i), label=f"L={L}")
        a2.plot(t[m] / L**1.5, W[m] / L**0.5, color=cmap(i), label=f"L={L}")
    for ax, ttl in ((a1, fr"best fit $\alpha={ab:.2f},\ z={zb:.2f}$"),
                    (a2, r"KPZ $\alpha=1/2,\ z=3/2$")):
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel(r"$t/L^{z}$"); ax.set_ylabel(r"$W/L^{\alpha}$")
        ax.set_title(ttl); ax.legend()
    fig.suptitle(r"$N=6$ (just above $N_c\approx5.24$)", fontsize=17)
    save(fig, "collapse_N6")


# --------------------------------------------------------------------------
def fig_family_vicsek(res):
    fig, ax = plt.subplots(figsize=(7.5, 7))
    Ns = list(res)
    cmap = plt.get_cmap("plasma", len(Ns) + 1)
    lo, hi = 0.10, 0.40
    ax.plot([lo, hi], [lo, hi], "k--", lw=1.5, label=r"$\alpha=\beta z$")
    for j, N in enumerate(Ns):
        a_s, bz_s = measure_fv(N)
        if np.isnan(a_s).all():
            continue
        ax.errorbar(np.nanmean(bz_s), np.nanmean(a_s),
                    xerr=np.nanstd(bz_s), yerr=np.nanstd(a_s),
                    fmt="o", ms=10, capsize=4, color=cmap(j), label=f"$N={N}$")
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi); ax.set_aspect("equal")
    ax.set_xlabel(r"$\beta\,z$  (growth $\times$ dynamic)")
    ax.set_ylabel(r"$\alpha$  (roughness)")
    ax.legend(title="Family–Vicsek check", ncol=2, fontsize=11)
    save(fig, "family_vicsek")


# --------------------------------------------------------------------------
def fig_growth_anatomy():
    fig, ax = plt.subplots(figsize=(8.5, 6.5))
    colours = {6: "#0d0887", 10: "#d8576b", 20: "#f0a202"}
    for N, c in colours.items():
        t, W = mean_W(128, N)
        if t is None:
            continue
        m = t > 0
        ax.plot(t[m], W[m], color=c, label=f"N={N}, L=128")
    tt = np.array([2, 40]); ax.plot(tt, 1.35 * tt**0.5, "k--", lw=1.4, label=r"$t^{1/2}$ (random deposition)")
    tt2 = np.array([3, 400]); ax.plot(tt2, 2.1 * tt2**0.30, "k:", lw=1.4, label=r"$t^{0.30}$")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("time $t$"); ax.set_ylabel("roughness $W$")
    ax.set_xlim(1, 200000); ax.set_ylim(1, 40)
    ax.legend()
    save(fig, "growth_anatomy")


# Large-L runs (beyond the crossover l*). L -> steps, per N.
# l* = Gaussianisation length, measured independently from the increment kurtosis.
LARGE_L_SPEC = {
    6:  ({16: 20000, 32: 20000, 64: 20000, 128: 20000, 256: 20000,
          512: 40000, 1024: 80000, 2048: 250000},
         20.0, [16, 32, 64, 128, 256], [256, 512, 1024, 2048]),
    10: ({16: 60000, 32: 60000, 64: 60000, 128: 60000, 256: 150000,
          512: 250000, 1024: 1200000},
         31.0, [16, 32, 64, 128], [256, 512, 1024]),
}


def _meanW_spec(L, N, spec):
    import glob as _g
    from common import roughness_series, DATA
    fs = sorted(_g.glob(f"{DATA}/L_{L}_N_{N}_steps_{spec[L]}_sim_*.tsv"))
    rows, ref = [], None
    for f in fs:
        tt, ww = roughness_series(f)
        if ref is None:
            ref = tt
        rows.append(ww[:min(len(ref), len(ww))])
    n = min(len(r) for r in rows)
    return ref[:n], np.array([r[:n] for r in rows]).mean(0)


def _cost(data, Ls, a, z):
    cs = []
    for L in Ls:
        t, W = data[L]
        m = t > 0
        cs.append((np.log(t[m] / L**z), np.log(W[m] / L**a)))
    lo = max(c[0][0] for c in cs); hi = min(c[0][-1] for c in cs)
    if hi <= lo:
        return np.inf
    g = np.linspace(lo, hi, 120)
    ys = np.array([np.interp(g, c[0], c[1]) for c in cs])
    return float(np.mean(np.var(ys, axis=0)))


def fig_collapse_largeL(N=6):
    """Family-Vicsek collapse restricted to L >> l*, testing KPZ vs EW.
    Below l* the interface has no continuum description, so no asymptotic
    collapse should be expected there -- and indeed KPZ only works above l*."""
    spec, lstar, small, large = LARGE_L_SPEC[N]
    data = {L: _meanW_spec(L, N, spec) for L in spec}
    panels = [(small, 0.5, 1.5, r"KPZ, $L\lesssim\ell^*$  (fails)"),
              (large, 0.5, 1.5, r"KPZ, $L\gg\ell^*$  (collapses)"),
              (large, 0.5, 2.0, r"EW, $L\gg\ell^*$  (excluded)")]
    fig, axes = plt.subplots(1, 3, figsize=(19, 5.6))
    for ax, (Ls, a, z, title) in zip(axes, panels):
        cmap = plt.get_cmap("viridis", len(Ls))
        for i, L in enumerate(Ls):
            t, W = data[L]
            m = t > 0
            ax.plot(t[m] / L**z, W[m] / L**a, color=cmap(i), label=f"L={L}")
        c = _cost(data, Ls, a, z)
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel(r"$t/L^{z}$"); ax.set_ylabel(r"$W/L^{\alpha}$")
        ax.set_title(f"{title}\n" fr"$\alpha={a},\ z={z}$   cost$={c:.4f}$", fontsize=14)
        ax.legend(fontsize=11)
    fig.suptitle(fr"$N={N}$: the KPZ collapse only works beyond the crossover "
                 fr"$\ell^*\approx{lstar:.0f}$", fontsize=17)
    save(fig, f"collapse_N{N}_largeL")


def fig_universality_plane(res):
    """Exponents in 'universality-class space'. Family-Vicsek (beta = alpha/z) leaves
    only two of the three independent, so each class is a POINT in the (alpha, z)
    plane, and rays through the origin are lines of constant beta."""
    fig, ax = plt.subplots(figsize=(9, 7))
    # constant-beta rays: z = alpha / beta
    for beta, lab in [(0.5, r"$\beta=1/2$"), (1 / 3, r"$\beta=1/3$"), (0.25, r"$\beta=1/4$")]:
        aa = np.linspace(0, 0.62, 10)
        ax.plot(aa, aa / beta, ls=":", color="grey", lw=1.2, zorder=1)
        al = min(0.60, 2.15 * beta)
        ax.text(al, al / beta, lab, color="grey", fontsize=12, ha="left", va="bottom")
    # reference classes (1+1 dimensions)
    ax.plot(0.5, 1.5, "*", ms=24, color="black", zorder=5)
    ax.annotate("KPZ", (0.5, 1.5), xytext=(10, -4), textcoords="offset points", fontsize=15)
    ax.plot(0.5, 2.0, "*", ms=24, color="dimgrey", zorder=5)
    ax.annotate("EW", (0.5, 2.0), xytext=(10, -4), textcoords="offset points", fontsize=15)
    ax.plot(0, 0, "s", ms=12, color="crimson", zorder=5)
    ax.annotate("random deposition\n" r"($\alpha\to0,\ \beta\to1/2$)", (0, 0),
                xytext=(12, 8), textcoords="offset points", fontsize=12, color="crimson")
    # our model's trajectory as N varies
    Ns = list(res)
    a = np.array([np.nanmean(res[N]["alpha"]) for N in Ns])
    z = np.array([np.nanmean(res[N]["z"]) for N in Ns])
    da = np.array([np.nanstd(res[N]["alpha"]) for N in Ns])
    dz = np.array([np.nanstd(res[N]["z"]) for N in Ns])
    ax.errorbar(a, z, xerr=da, yerr=dz, fmt="none", ecolor="0.6", lw=1, zorder=2)
    ax.plot(a, z, "-", color="0.45", lw=1.3, zorder=3)
    sc = ax.scatter(a, z, c=Ns, cmap="plasma", s=130, zorder=4, edgecolor="k", linewidth=0.6)
    for N, x, y in zip(Ns, a, z):
        ax.annotate(f"{N}", (x, y), xytext=(7, 7), textcoords="offset points", fontsize=11)
    cb = fig.colorbar(sc, ax=ax)
    cb.set_label("number of colours $N$")
    ax.set_xlabel(r"roughness exponent $\alpha$")
    ax.set_ylabel(r"dynamic exponent $z$")
    ax.set_xlim(0, 0.65); ax.set_ylim(0, 2.25)

    # zoomed inset on our trajectory (it is a tight cluster on the full scale)
    axi = ax.inset_axes([0.06, 0.55, 0.40, 0.40])
    for beta in (0.5, 1 / 3, 0.25):
        aa = np.linspace(0, 0.62, 10)
        axi.plot(aa, aa / beta, ls=":", color="grey", lw=1.1)
    axi.errorbar(a, z, xerr=da, yerr=dz, fmt="none", ecolor="0.6", lw=1)
    axi.plot(a, z, "-", color="0.45", lw=1.2)
    axi.scatter(a, z, c=Ns, cmap="plasma", s=70, edgecolor="k", linewidth=0.5, zorder=4)
    for N, x, y in zip(Ns, a, z):
        axi.annotate(f"{N}", (x, y), xytext=(5, 4), textcoords="offset points", fontsize=9)
    axi.set_xlim(0.18, 0.35); axi.set_ylim(0.45, 1.10)
    axi.tick_params(labelsize=9)
    axi.set_title("zoom on our model", fontsize=11)
    ax.indicate_inset_zoom(axi, edgecolor="black")
    save(fig, "universality_plane")


def fig_exponent_flow(res):
    """alpha, beta, z as functions of N: how the three exponents flow together
    from the near-critical point toward the random-deposition limit."""
    Ns = np.array(list(res), float)
    fig, ax = plt.subplots(figsize=(8.5, 6))
    for key, col, lab in [("alpha", "#b5179e", r"$\alpha$ (roughness)"),
                          ("beta", "#4361ee", r"$\beta$ (growth)"),
                          ("z", "#2a9d8f", r"$z$ (dynamic)")]:
        m = np.array([np.nanmean(res[N][key]) for N in res])
        e = np.array([np.nanstd(res[N][key]) for N in res])
        ax.errorbar(Ns, m, yerr=e, fmt="o-", ms=8, capsize=3, color=col, label=lab)
    ax.axhline(0.5, ls="--", color="#4361ee", lw=1, alpha=0.6)
    ax.axhline(0.0, ls="--", color="#b5179e", lw=1, alpha=0.6)
    ax.text(Ns.max() * 0.72, 0.52, r"RD: $\beta\to1/2$", color="#4361ee", fontsize=12)
    ax.text(Ns.max() * 0.72, 0.02, r"RD: $\alpha\to0$", color="#b5179e", fontsize=12)
    ax.set_xlabel("number of colours $N$"); ax.set_ylabel("exponent")
    ax.set_ylim(-0.05, 1.05)
    ax.legend(loc="upper right")
    save(fig, "exponent_flow")


if __name__ == "__main__":
    res = measure_all()
    fig_universality_plane(res)
    fig_exponent_flow(res)
    fig_beta_eff()
    fig_slope_collapse()
    fig_slope_scale()
    fig_roughness_vs_L(res)
    fig_alpha_vs_N(res)
    fig_collapse_N6()
    fig_family_vicsek(res)
    fig_growth_anatomy()
