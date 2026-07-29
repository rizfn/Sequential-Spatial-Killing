"""Regenerate all avalanche-scaling figures into themed subfolders of plots/.

  mechanism_...        which variable drives it: NOT the local slope, but the
                       cascade extent w, via k = a/d (mechanism/)
  sawtooth_mechanism   the fractional species is a self-poisoning impurity
                       (speciesSawtooth/)
  transition           v(N): the order parameter of the growth-arrest
                       transition, vanishing at N_c = 5.075, + h ~ t^0.32 at N_c
  critical_test        THE decisive plot: <s^2>/<s> vs L *at* N_c.  Flat =>
                       avalanches are not scale-free even at the critical point
  mass_balance         v = 1 - <s> f_active, why arrest needs no divergence
  cutoff_vs_N          the avalanche scale vs N: sawtooth peaking at integer N
  steady_state         <s^2>/<s> vs t: the runs are equilibrated, so flatness is
                       physics and not a transient artifact
  avalanche_pdf_vs_L   P(s) for every L at N=6 -- curves lie on top of each other
  cutoff_vs_L          cutoff estimators vs L, against what SOC would predict
  clusters_duration    the companion cascade observables
"""
import os

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from common import (LS, MAIN_N, NCRIT, NCRIT_LS, NSAW, NSAW_L, NSCAN, NSCAN_LS,
                    NVEL, NVEL_L, N_C, PLOTS, PLOT_DIRS, STEPS, WARMUP_FRAC,
                    active_fraction, cutoff_probe, impurity_enrichment, logbin,
                    mean_and_moment, moments_vs_t, pdf, pooled, slope_resolved,
                    velocity)

plt.rcParams.update({
    "font.size": 17, "axes.labelsize": 21, "xtick.labelsize": 15,
    "ytick.labelsize": 15, "legend.fontsize": 13, "lines.linewidth": 2.1,
    "figure.autolayout": True,
})


def save(fig, name, sub):
    """Write into a themed subfolder of plots/ (sub must be in PLOT_DIRS)."""
    assert sub in PLOT_DIRS, sub
    d = f"{PLOTS}/{sub}"
    os.makedirs(d, exist_ok=True)
    fig.savefig(f"{d}/{name}.svg")
    fig.savefig(f"{d}/{name}.png", dpi=150)
    plt.close(fig)
    print(f"saved {sub}/{name}")


def quantile(value, counts, q):
    cdf = np.cumsum(counts.astype(float)) / counts.sum()
    return float(value[np.searchsorted(cdf, q)])


# a readable subset for the busy per-L PDF panels (16 sizes is too many curves)
LS_SHOW = [L for L in LS if L in (16, 32, 64, 128, 256, 512, 1024, 2048, 4096)]




# --------------------------------------------------------------------------
def fig_mechanism():
    """WHICH variable drives the stretched exponential?

    Left: the local slope m at the deposition site.  P(m) is a clean exponential
    -- but <s|m> SATURATES (~1.9), so the slope does NOT set the mass.  The
    naive "s ~ m^d" route is refuted directly: at m=12 it would predict s~600.

    Right: the cascade's spatial extent w = the COUNT of distinct columns it
    eliminates in (not the span!).  Here <s|w> ~ w^d with d ~ 1.7 (nearly compact
    in 2D, as geometry requires), and P(w) = A w^-t exp(-w/w0) -- a power law
    times a PURE exponential, correlation length w0 ~ 5 columns.  Carrying the
    prefactor through the change of variables gives tau = 1+(t-1)/d but leaves
    k = 1/d.  That predicts k=0.57 against 0.2-0.39 measured: the chain does not
    close, because s = A w^d is only a statement about the MEAN <s|w> while the
    change of variables needs it deterministic.  The conditional scatter fattens
    the tail of P(s).
    """
    L = 1024
    m, pm, s_of_m, w, pw, s_of_w = slope_resolved(L, MAIN_N)
    if m is None:
        print("skip mechanism: no slopeResolved data")
        return
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    ax = axes[0]
    sel = pm > 1e-6
    ax.semilogy(m[sel], pm[sel], "o-", color="navy", ms=5, label="$P(m)$")
    k = sel & (np.abs(m) >= 1) & (np.abs(m) <= 12)
    lam = -np.polyfit(np.abs(m[k]), np.log(pm[k]), 1)[0]
    ax.set_xlabel("local slope $m$ at deposition site")
    ax.set_ylabel("$P(m)$", color="navy")
    ax.set_xlim(-16, 16)
    ax2 = ax.twinx()
    ax2.plot(m[sel], s_of_m[sel], "s-", color="crimson", ms=5)
    ax2.set_ylabel(r"$\langle s\,|\,m\rangle$", color="crimson")
    ax2.set_ylim(0, 3)
    ax2.plot(np.arange(1, 13), 2.0 * np.arange(1, 13) ** 2.55 / 2.0 ** 2.55, ":",
             color="grey", lw=2)
    ax2.text(3.5, 2.5, r"$m^{2.55}$ (what $s\sim m^d$" "\n" r"would need) $\rightarrow$ off scale",
             fontsize=11, color="grey")
    ax.set_title(rf"slope: $P(m)\sim e^{{-{lam:.2f}|m|}}$ but $\langle s|m\rangle$ saturates",
                 fontsize=14)

    ax = axes[1]
    good = (pw > 0) & (w >= 1)
    ax.loglog(w[good], s_of_w[good], "o", color="crimson", ms=7, label=r"$\langle s\,|\,w\rangle$")
    kk = good & (w >= 2) & (w <= 30)
    d_fit, c0 = np.polyfit(np.log(w[kk]), np.log(s_of_w[kk]), 1)
    wf = np.logspace(0, np.log10(w[good].max()), 50)
    ax.plot(wf, np.exp(c0) * wf ** d_fit, "-", color="black", lw=1.8,
            label=rf"$s\sim w^{{{d_fit:.2f}}}$")
    ax.set_xlabel("cascade extent $w$ (columns)")
    ax.set_ylabel(r"$\langle s\,|\,w\rangle$", color="crimson")
    ax3 = ax.twinx()
    ax3.loglog(w[good], pw[good], "^-", color="navy", ms=4, alpha=0.8)
    ax3.set_ylabel("$P(w)$", color="navy")
    # P(w) is a power law TIMES a PURE exponential -- fitting a Weibull without
    # the prefactor fakes a stretch (it reads a~0.66); with the prefactor the
    # cutoff is exponential (a~1) and the correlation length is ~5 columns.
    gw = good & (w >= 3) & (pw * 1e9 > 50)
    (Aw, tw, bw), _ = curve_fit(lambda x, A, t, b: A - t * np.log(x) - b * x,
                                w[gw], np.log(pw[gw]), p0=[0., 1., 0.3], maxfev=200000)
    ax.set_title(rf"extent: $s\sim w^{{{d_fit:.2f}}}$;  $P(w)=Aw^{{-{tw:.2f}}}e^{{-w/{1/bw:.1f}}}$"
                 "\n" rf"$\Rightarrow k=1/d={1/d_fit:.2f}$ (measured $0.2$-$0.39$)", fontsize=14)
    ax.legend(loc="upper left", fontsize=12)
    save(fig, "mechanism_slope_vs_extent", "mechanism")


# --------------------------------------------------------------------------
def fig_sawtooth_mechanism():
    """Why the sawtooth: the fractional species is a self-poisoning impurity.

    At N = n + f, the last species carries weight f, so it is deposited with
    probability f/N.  Being rare, it seldom finds a same-species neighbour, so it
    is seldom eliminated -- and therefore ACCUMULATES.  Left: its pile fraction
    divided by its deposition fraction is up to ~2.5x, i.e. it is a frozen defect
    that fragments clusters and suppresses avalanches.  The effect must vanish at
    BOTH ends of a tooth (f=0 and f=1 are both uniform integer systems), which is
    exactly the sawtooth shape -- confirmed at right by resolving one tooth.
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    ax = axes[0]
    xs, ys = [], []
    for N in NSAW:
        e = impurity_enrichment(NSAW_L, N)
        if e is not None:
            xs.append(N); ys.append(e)
    ax.plot(xs, ys, "o-", color="crimson", ms=8)
    ax.axhline(1.0, ls="--", color="grey", lw=1.5)
    ax.text(6.55, 1.03, "1.0 = no enrichment", color="grey", fontsize=13)
    for n in (6, 7):
        ax.axvline(n, ls=":", color="black", lw=1.2)
    ax.set_xlabel("$N$")
    ax.set_ylabel("impurity enrichment\n(pile fraction / deposit fraction)")
    ax.set_title("the rare species accumulates", fontsize=15)

    ax = axes[1]
    xs, ys = [], []
    for N in NSAW:
        c = cutoff_probe(NSAW_L, N)
        if c is not None:
            xs.append(N); ys.append(c)
    ax.plot(xs, ys, "o-", color="navy", ms=8)
    for n in (6, 7):
        ax.axvline(n, ls=":", color="black", lw=1.2)
    ax.text(6.02, min(ys), "integer $N$:\nno impurity", fontsize=12)
    ax.set_xlabel("$N$")
    ax.set_ylabel(r"$\langle s^2\rangle/\langle s\rangle$")
    ax.set_title(f"one tooth resolved ($L={NSAW_L}$)", fontsize=15)
    save(fig, "sawtooth_mechanism", "speciesSawtooth")


# --------------------------------------------------------------------------
def fig_transition():
    """The growth-arrest transition is real: v(N) is a genuine order parameter.

    Left: v vanishes continuously at N_c.  Right: at N_c the pile still grows,
    but as a sublinear power law h ~ t^0.32 rather than h ~ vt -- the scaling
    signature of a critical point.  So this IS a critical point; the next figure
    shows the avalanches nevertheless do not know about it.
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    ax = axes[0]
    Ns = np.array(NVEL)
    v = np.array([velocity(NVEL_L, N) for N in NVEL], dtype=float)
    ax.plot(Ns, v, "o-", color="navy", ms=7)
    ax.axvline(N_C, ls="--", color="crimson", lw=2)
    ax.axhline(0, ls=":", color="grey", lw=1.5)
    ax.text(N_C + 0.004, max(v) * 0.55, f"$N_c={N_C}$", color="crimson", fontsize=15)
    ax.text(5.005, max(v) * 0.8, "arrested\n$v=0$", fontsize=14, color="grey")
    ax.set_xlabel("$N$")
    ax.set_ylabel(r"growth velocity $v=\mathrm{d}\langle h\rangle/\mathrm{d}t$")
    ax.set_title(f"$L={NVEL_L}$", fontsize=16)

    ax = axes[1]
    for N, col in ((NCRIT[0], "crimson"), (5.20, "darkorange"), (MAIN_N, "navy")):
        tg, m1, m21, af, ta, mh = moments_vs_t(4096 if N in NCRIT else 1024, N)
        if tg is None:
            continue
        m = ta > 0
        ax.plot(ta[m], mh[m], "-", color=col, lw=2.2, label=f"$N={N:.3f}$")
    # slope guides
    t = np.logspace(2, 4.3, 50)
    ax.plot(t, 0.9 * t ** 0.32, "--", color="black", lw=1.5)
    ax.text(2e3, 0.9 * 2e3 ** 0.32 * 0.45, r"$t^{0.32}$ (at $N_c$)", fontsize=13)
    ax.plot(t, 0.06 * t, ":", color="grey", lw=1.8)
    ax.text(3e3, 0.06 * 3e3 * 1.3, r"$t^{1}$ (free growth)", fontsize=13, color="grey")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("time $t$")
    ax.set_ylabel(r"mean height $\langle h\rangle$")
    ax.legend(loc="upper left")
    save(fig, "transition_vN", "criticality")


# --------------------------------------------------------------------------
def fig_critical_test():
    """THE test. At a critical point xi diverges, so the ONLY cutoff left is L
    and <s^2>/<s> ~ L^D must grow without bound.  It does not: flat from L=128
    to L=4096 at N_c.  The avalanches are not scale-free even at criticality.
    """
    fig, ax = plt.subplots(figsize=(9, 6.5))
    for N, col, mk in ((NCRIT[0], "crimson", "o"), (NCRIT[1], "darkorange", "s")):
        y = [cutoff_probe(L, N) for L in NCRIT_LS]
        ax.plot(NCRIT_LS, y, mk + "-", color=col, ms=8,
                label=f"$N={N:.3f}$" + ("  ($=N_c$)" if N == NCRIT[0] else "  (just above)"))
    y6 = [cutoff_probe(L, MAIN_N) for L in NCRIT_LS]
    ax.plot(NCRIT_LS, y6, "d-", color="navy", ms=7, label=f"$N={MAIN_N:.1f}$ (off-critical)")
    L = np.array(NCRIT_LS, float)
    ax.plot(L, 8.5 * (L / L[0]) ** 0.5, "--", color="black", lw=1.5)
    ax.text(L[2], 8.5 * (L[2] / L[0]) ** 0.5 * 1.1,
            r"$\propto L^{1/2}$: what criticality" "\n" r"would require", fontsize=13)
    ax.set_xscale("log", base=2); ax.set_yscale("log")
    ax.set_ylim(6, 60)
    ax.set_xlabel("system size $L$")
    ax.set_ylabel(r"$\langle s^2\rangle/\langle s\rangle$")
    ax.legend(loc="upper left")
    save(fig, "critical_test", "criticality")


# --------------------------------------------------------------------------
def fig_mass_balance():
    """Why arrest needs no diverging avalanche.

    Every deposition adds exactly 1 block and removes f_active*<s> on average, so
    v = 1 - <s> f_active exactly.  Arrest (v=0) is therefore a condition on the
    FIRST moment, which the swarm of s=2 events already dominates.  Nothing has
    to diverge.  (In the BTW sandpile, boundary-only dissipation forces
    <s> ~ L^2 instead -- which is what makes sandpiles critical.)
    """
    Ns = [n for n in sorted(set(NVEL) | set(NSCAN) | set(NCRIT)) if n <= 6.0]
    xs, pred, meas = [], [], []
    for N in Ns:
        L = NVEL_L
        c = cutoff_probe(L, N)
        if c is None:
            continue
        v_meas = velocity(L, N)
        v_p, h, _ = pooled(L, N)
        m1 = mean_and_moment(v_p, h["mass"].astype(float), 1)
        fa = active_fraction(L, N)
        if v_meas is None or fa is None:
            continue
        xs.append(N); pred.append(1 - m1 * fa); meas.append(v_meas)

    fig, ax = plt.subplots(figsize=(9, 6.5))
    ax.plot(xs, pred, "-", color="crimson", lw=2.5, label=r"$1-\langle s\rangle f_{\rm active}$")
    ax.plot(xs, meas, "o", color="navy", ms=8, label=r"$v$ measured from $\langle h\rangle(t)$")
    ax.axvline(N_C, ls="--", color="grey", lw=1.5)
    ax.axhline(0, ls=":", color="grey", lw=1.5)
    ax.text(N_C + 0.01, max(pred) * 0.6, f"$N_c={N_C}$\n$\\langle s\\rangle f_{{\\rm active}}=1$",
            fontsize=13, color="grey")
    ax.set_xlabel("$N$")
    ax.set_ylabel("growth velocity $v$")
    ax.legend(loc="upper left")
    save(fig, "mass_balance", "criticality")


# --------------------------------------------------------------------------
def fig_cutoff_vs_N():
    """The avalanche scale vs N: a sawtooth with maxima at INTEGER N.

    A non-integer N means a final species carrying only the fractional weight: a
    rare "impurity" block that seldom finds a partner and so fragments clusters.
    Crossing an integer introduces a fresh impurity and the scale drops.  The
    global max at N=6 is just the largest integer peak above N_c -- not a
    critical point, as fig_critical_test shows.
    """
    fig, ax = plt.subplots(figsize=(9, 6.5))
    Ns = sorted(set(NSCAN) | set(NCRIT) | {5.10, 5.15})
    cmap = plt.get_cmap("viridis", len(NSCAN_LS))
    for i, L in enumerate(NSCAN_LS):
        xs, ys = [], []
        for N in Ns:
            c = cutoff_probe(L, N)
            if c is not None:
                xs.append(N); ys.append(c)
        ax.plot(xs, ys, "o-", color=cmap(i), ms=6, label=f"$L={L}$")
    for n in (6, 7, 8):
        ax.axvline(n, ls=":", color="grey", lw=1.2)
    ax.axvline(N_C, ls="--", color="crimson", lw=2)
    ax.text(N_C + 0.05, 7.2, f"$N_c$", color="crimson", fontsize=15)
    ax.set_xlabel("$N$")
    ax.set_ylabel(r"$\langle s^2\rangle/\langle s\rangle$")
    ax.legend(title="dotted = integer $N$", loc="lower right")
    save(fig, "cutoff_vs_N", "speciesSawtooth")


# --------------------------------------------------------------------------
def fig_steady_state():
    """Show the runs are equilibrated: <s^2>/<s> plateaus long before WARMUP.

    Without this panel the flat cutoff_vs_L could be dismissed as every system
    being measured in a common transient.  The plateau is reached by t ~ L^z
    (z ~ 0.95, dotted marks), and the histogrammed window starts at WARMUP,
    far to the right of it.
    """
    fig, ax = plt.subplots(figsize=(9, 6.5))
    cmap = plt.get_cmap("viridis", len(LS_SHOW))
    for i, L in enumerate(LS_SHOW):
        t, m1, m21, af, ta, mh = moments_vs_t(L, MAIN_N)
        if t is None:
            continue
        ax.plot(t, m21, "-", color=cmap(i), label=f"$L={L}$")
        tsat = L ** 0.95
        ax.plot([tsat], [np.interp(tsat, t, m21)], "o", color=cmap(i), ms=7,
                mec="black", mew=0.8, zorder=5)
    ax.axvline(WARMUP_FRAC * STEPS, ls="--", color="crimson", lw=2)
    ax.text(WARMUP_FRAC * STEPS * 0.85, 2.4, "warmup: histograms start here  ",
            color="crimson", fontsize=13, va="bottom", ha="right")
    ax.annotate("", xy=(STEPS, 2.3), xytext=(WARMUP_FRAC * STEPS, 2.3),
                arrowprops=dict(arrowstyle="-|>", color="crimson", lw=2))
    ax.set_xscale("log")
    ax.set_xlabel("time $t$ (steps per lattice point)")
    ax.set_ylabel(r"$\langle s^2\rangle/\langle s\rangle$")
    ax.legend(ncol=2, title=f"$N={MAIN_N}$;  $\\bullet$ = $t_{{\\rm sat}}\\sim L^{{0.95}}$",
              loc="upper left")
    save(fig, f"steady_state_N_{MAIN_N:.1f}", "finiteSize")


# --------------------------------------------------------------------------
def fig_pdf_vs_L():
    """The headline: P(s) is independent of L."""
    fig, ax = plt.subplots(figsize=(9, 6.5))
    cmap = plt.get_cmap("viridis", len(LS_SHOW))
    for i, L in enumerate(LS_SHOW):
        v, h, _ = pooled(L, MAIN_N)
        if v is None:
            continue
        x, y = logbin(v, h["mass"])
        ax.plot(x, y, "o-", ms=4.5, color=cmap(i), label=f"$L={L}$")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("avalanche mass $s$")
    ax.set_ylabel("$P(s)$")
    ax.legend(ncol=2, title=f"$N={MAIN_N:.1f}$")
    save(fig, f"avalanche_pdf_vs_L_N_{MAIN_N:.1f}", "sizeDistribution")


# --------------------------------------------------------------------------
def fig_cutoff_vs_L():
    """Cutoff estimators vs L. Flat = no size scaling.

    <s^2>/<s> is the standard cutoff probe: for P ~ s^-tau f(s/s_c) it grows
    like s_c, and it is far less noisy than a high quantile or s_max.
    """
    Ls, m1, m21, q999, smax = [], [], [], [], []
    for L in LS:
        v, h, _ = pooled(L, MAIN_N)
        if v is None:
            continue
        c = h["mass"].astype(float)
        Ls.append(L)
        m1.append(mean_and_moment(v, c, 1))
        m21.append(mean_and_moment(v, c, 2) / mean_and_moment(v, c, 1))
        q999.append(quantile(v, c, 0.999))
        smax.append(float(v[c > 0][-1]))
    Ls = np.array(Ls, dtype=float)

    # the cutoff probe converges to a finite constant: fit m21 = A - B L^-x
    m21 = np.array(m21)
    fs = lambda L, A, B, x: A - B * L ** -x
    p, cov = curve_fit(fs, Ls, m21, p0=[10.6, 10, 1.0])
    err = np.sqrt(np.diag(cov))

    fig, ax = plt.subplots(figsize=(9, 6.5))
    ax.plot(Ls, m21, "o", color="crimson", ms=8, label=r"$\langle s^2\rangle/\langle s\rangle$")
    Lf = np.logspace(np.log10(Ls[0]), np.log10(Ls[-1] * 1.5), 100)
    ax.plot(Lf, fs(Lf, *p), "-", color="crimson", lw=1.6)
    ax.axhline(p[0], ls=":", color="crimson", lw=1.6)
    ax.text(Ls[0] * 1.1, p[0] * 1.06,
            rf"$A={p[0]:.2f}\pm{err[0]:.2f}$ (finite!),  corr. $\propto L^{{-{p[2]:.2f}}}$",
            color="crimson", fontsize=13)
    ax.plot(Ls, q999, "s-", color="darkorange", label=r"$s_{99.9\%}$")
    ax.plot(Ls, smax, "^-", color="grey", label=r"$s_{\max}$ (sample-size artifact)")
    ax.plot(Ls, m1, "d-", color="navy", label=r"$\langle s\rangle$")
    # what a size-limited (SOC-like) cutoff would look like, anchored at L=16
    ax.plot(Ls, m21[0] * (Ls / Ls[0]) ** 1.0, "--", color="black", lw=1.5)
    ax.text(Ls[-4], m21[0] * (Ls[-4] / Ls[0]) ** 1.0 * 1.25,
            r"$\propto L$ (SOC)", fontsize=13, rotation=30)
    ax.set_xscale("log", base=2); ax.set_yscale("log")
    ax.set_xlabel("system size $L$")
    ax.set_ylabel("avalanche mass scale")
    ax.set_ylim(2, 3e3)
    ax.legend(title=f"$N={MAIN_N:.1f}$", loc="upper left")
    save(fig, f"cutoff_vs_L_N_{MAIN_N:.1f}", "finiteSize")


# --------------------------------------------------------------------------
def fig_pdf_vs_N():
    Ns = sorted(set(NSCAN) | set(NCRIT))
    L = 1024
    fig, ax = plt.subplots(figsize=(9, 6.5))
    cmap = plt.get_cmap("plasma", len(Ns) + 1)
    for i, N in enumerate(Ns):
        v, h, _ = pooled(L, N)
        if v is None:
            continue
        x, y = logbin(v, h["mass"])
        ax.plot(x, y, "o-", ms=4.5, color=cmap(i), label=f"$N={N:.3f}$")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("avalanche mass $s$")
    ax.set_ylabel("$P(s)$")
    ax.legend(title=f"$L={L}$", ncol=2, fontsize=11)
    save(fig, f"avalanche_pdf_vs_N_L_{L}", "sizeDistribution")


# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
def fig_clusters_duration():
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    cmap = plt.get_cmap("viridis", len(LS))
    for ax, key, lab in zip(axes, ("clusters", "duration"),
                            ("clusters eliminated $n$", "cascade duration $d$")):
        for i, L in enumerate(LS):
            v, h, _ = pooled(L, MAIN_N)
            if v is None:
                continue
            x, y = pdf(v, h[key])
            ax.plot(x, y, "o-", ms=4.5, color=cmap(i), label=f"$L={L}$")
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel(lab)
        ax.set_ylabel("probability")
    axes[0].legend(ncol=2, title=f"$N={MAIN_N:.1f}$")
    save(fig, f"clusters_duration_N_{MAIN_N:.1f}", "sizeDistribution")


if __name__ == "__main__":
    fig_mechanism()
    fig_sawtooth_mechanism()
    fig_transition()
    fig_critical_test()
    fig_mass_balance()
    fig_cutoff_vs_N()
    fig_steady_state()
    fig_pdf_vs_L()
    fig_cutoff_vs_L()
    fig_pdf_vs_N()
    fig_clusters_duration()
