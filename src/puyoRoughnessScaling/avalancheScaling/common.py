"""Shared config + data-loading for the avalanche-size-distribution analysis.

avalancheDist.cpp writes two files per simulation.

outputs/avalancheDist/ : three histograms over a shared value axis,
    value  mass  clusters  duration
where for each deposition (after warmup) we record the whole triggered cascade:
    mass     s : blocks eliminated before the pile is quiet again
    clusters n : same-species components eliminated
    duration d : chain generations that eliminated at least one block

outputs/moments/ : the moments of s, and the mean pile height, in log-spaced
    time windows over ALL t (transient included).  Two jobs: verify steady state
    (<s^2>/<s> must plateau well before WARMUP), and measure the growth velocity
    v = d<h>/dt, the order parameter of the growth-arrest transition.

THE PHYSICS (see readme.md):
  * N is CONTINUOUS.  v(N) vanishes at N_c = 5.075(10): a genuine continuous
    transition, with h ~ t^0.32 at N_c instead of h ~ vt.
  * But avalanches are NOT scale-free, even at N_c: <s^2>/<s> is flat in L from
    128 to 4096.  The reason is a mass balance, exact per deposition:
        v = 1 - <s> * f_active
    Arrest is a condition on the FIRST moment, which is dominated by the tiny
    avalanches (P(s=2) ~ 0.74) and needs no divergence to reach 1.  Contrast the
    BTW sandpile, where dissipation is only at the boundary so <s> ~ L^2 is
    FORCED to diverge.  Puyo dissipates in the bulk, so nothing forces it.
"""
import glob
import os
import re
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = f"{HERE}/outputs/avalancheDist"
MOM = f"{HERE}/outputs/moments"
SLOPE = f"{HERE}/outputs/slopeResolved"
PLOTS = f"{HERE}/plots"

# plots/ is themed into subfolders; figures.py's save() takes one of these.
PLOT_DIRS = ("criticality", "sizeDistribution", "mechanism", "speciesSawtooth",
             "finiteSize", "integerN")

# ---------------------------------------------------------------------------
# Sweep definition (single source of truth for generation AND analysis).
#
# STEPS is fixed, not L-dependent.  It was sized against ../exponents.py's
# z ~ 0.95 at N=6 (t_sat ~ L^0.95 is only ~2600 steps even at L=4096, so 32768
# is >12x that), but do NOT lean on that number: those exponents are known to
# suffer a finite-size crossover.  The load-bearing check is empirical --
# outputs/moments/ shows <s^2>/<s> plateauing by t~7000 at every L and every N,
# INCLUDING at N_c where critical slowing down is a real worry.  WARMUP=8192 is
# past the measured plateau, which is the claim that actually matters.
# ---------------------------------------------------------------------------
STEPS = 32768
WARMUP_FRAC = 0.25

N_C = 5.075          # growth-arrest transition, v(N_c) = 0
MAIN_N = 6.000       # largest integer-N avalanche peak; the reference family

# family A: cutoff scaling in L at fixed off-critical N
LS = [16, 24, 32, 48, 64, 96, 128, 192, 256, 384, 512, 768, 1024, 1536, 2048, 4096]
NSIMS = {16: 64, 24: 64, 32: 64, 48: 64, 64: 64, 96: 64, 128: 64,
         192: 32, 256: 32, 384: 32, 512: 32,
         768: 16, 1024: 16, 1536: 16, 2048: 16, 4096: 8}

# family B: the decisive test -- L scaling AT the critical point.  If the
# avalanches were critical, <s^2>/<s> ~ L^D would grow here without bound.
NCRIT = [5.075, 5.080]
NCRIT_LS = [128, 256, 512, 1024, 2048, 4096]
NCRIT_SIMS = 16

# family C: v(N) near N_c, to locate the transition (order parameter)
NVEL = [5.00, 5.02, 5.04, 5.05, 5.06, 5.07, 5.08, 5.10, 5.12, 5.15, 5.20, 5.30]
NVEL_L = 512
NVEL_SIMS = 8

# family D: N dependence of the avalanche scale.  Note the sawtooth: local
# maxima sit at INTEGER N, because a fractional species is a rare "impurity"
# that seldom matches and fragments clusters.
NSCAN = [5.20, 5.50, 6.00, 6.25, 6.50, 7.00, 7.50, 8.00]
NSCAN_LS = [64, 128, 256, 512, 1024]
NSCAN_SIMS = 16

# family E: one sawtooth tooth resolved finely, N = 6 -> 7.  Both endpoints are
# integers (fractional weight f = 0 and f = 1 are the same uniform system), so
# any dip in between is caused by the impurity species and must vanish at both
# ends.  This is what separates "impurity damage" from "smooth decay in N".
NSAW = [6.00, 6.05, 6.10, 6.20, 6.30, 6.40, 6.50, 6.60, 6.70, 6.80, 6.90, 6.95, 7.00]
NSAW_L = 512
NSAW_SIMS = 16


def nfmt(N):
    """N appears in filenames with fixed precision; keep C++ and python in sync."""
    return f"{float(N):.3f}"


def steps_of(L=None):
    """Steps per lattice point -- the natural time unit (total drops = L*steps)."""
    return STEPS


def box_H(N, steps=STEPS):
    """Box height, sized from the measured growth velocity v(N) with margin.

    v is ~0 at N_c and rises steeply past N=6.5, so a single ratio would either
    waste memory at low N or hit the ceiling at high N.  Every run reports
    ceiling_hits, which must be 0.
    """
    if N <= 5.5:
        r = 0.10      # v <= 0.05 here
    elif N <= 6.5:
        r = 0.30      # v(6) ~ 0.061
    else:
        r = 0.55      # v(7) ~ 0.216
    return int(np.ceil(r * steps)) + 512


def jobs():
    """Every (L, N, steps, sim, H, warmup) in the sweep, biggest box first.

    Biggest-first keeps the tall/wide (high memory) jobs from all landing at the
    end, and lets xargs backfill with cheap ones.
    """
    out = []
    seen = set()

    def add(L, N, nsims):
        for sim in range(nsims):
            key = (L, nfmt(N), sim)
            if key in seen:
                continue
            seen.add(key)
            out.append((L, nfmt(N), STEPS, sim, box_H(N), int(WARMUP_FRAC * STEPS)))

    for L in LS:
        add(L, MAIN_N, NSIMS[L])
    for N in NCRIT:
        for L in NCRIT_LS:
            add(L, N, NCRIT_SIMS)
    for N in NVEL:
        add(NVEL_L, N, NVEL_SIMS)
    for N in NSCAN:
        for L in NSCAN_LS:
            add(L, N, NSCAN_SIMS)
    for N in NSAW:
        add(NSAW_L, N, NSAW_SIMS)
    out.sort(key=lambda j: -j[0] * j[4])
    return out


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------
_META = re.compile(r"(\w+)=(-?[\d.]+)")


def load_sim(path):
    """Return (meta dict, value array, {name: counts}) for one simulation file."""
    meta, rows = {}, []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                meta.update({k: float(v) for k, v in _META.findall(line)})
            elif line.startswith("value"):
                continue
            else:
                rows.append([int(x) for x in line.split()])
    a = np.array(rows, dtype=np.int64)
    if a.size == 0:
        return meta, np.array([]), {}
    return meta, a[:, 0], dict(mass=a[:, 1], clusters=a[:, 2], duration=a[:, 3])


def pooled(L, N):
    """Pool histograms over all sims of (L, N).

    Returns (value, {name: counts}, drops_counted) or (None, None, 0).
    Pooling counts (not averaging densities) is right here: every deposition is
    an independent sample, so summing bins across sims is just a longer run, and
    the tail bins get the statistics they need.
    """
    files = sorted(glob.glob(f"{DATA}/L_{L}_N_{nfmt(N)}_steps_{STEPS}_sim_*.tsv"))
    if not files:
        return None, None, 0
    vmax, loaded, drops = 0, [], 0
    for f in files:
        meta, v, h = load_sim(f)
        if v.size == 0:
            continue
        drops += int(meta.get("drops_counted", 0))
        vmax = max(vmax, int(v[-1]))
        loaded.append((v, h))
    if not loaded:
        return None, None, 0
    tot = {k: np.zeros(vmax + 1, dtype=np.int64) for k in ("mass", "clusters", "duration")}
    for v, h in loaded:
        for k in tot:
            tot[k][v] += h[k]
    value = np.arange(vmax + 1)
    return value[1:], {k: c[1:] for k, c in tot.items()}, drops


def pdf(value, counts):
    """Normalized PDF over active avalanches, dropping empty bins."""
    m = counts > 0
    return value[m], counts[m] / counts.sum()


def logbin(value, counts, nbins=40):
    """Logarithmically binned PDF: essential for a clean tail.

    Linear bins in the tail hold 0 or 1 count each and the log-log plot becomes a
    fan of noise; log bins keep counts per bin roughly constant.
    Density = (counts in bin) / (total counts) / (bin width).
    """
    m = counts > 0
    v, c = value[m].astype(float), counts[m].astype(float)
    if v.size == 0:
        return np.array([]), np.array([])
    edges = np.unique(np.round(np.logspace(0, np.log10(v.max() + 1), nbins)).astype(int))
    edges = edges[edges >= 1]
    if edges.size < 2:
        return np.array([]), np.array([])
    idx = np.digitize(v, edges) - 1
    tot = c.sum()
    x, y = [], []
    for b in range(len(edges) - 1):
        sel = idx == b
        if not sel.any():
            continue
        w = edges[b + 1] - edges[b]
        x.append(np.exp(np.log(v[sel]).mean()))   # geometric centre of occupied values
        y.append(c[sel].sum() / tot / w)
    return np.array(x), np.array(y)


def mean_and_moment(value, counts, k=1):
    """<s^k> over active avalanches."""
    c = counts.astype(float)
    return (value.astype(float) ** k * c).sum() / c.sum()


def cutoff_probe(L, N):
    """<s^2>/<s>: the cutoff scale, free of any tau fit. None if no data.

    For P(s) ~ s^-tau f(s/s_c), <s^k> ~ s_c^(k-tau+1), so this ratio ~ s_c with
    the tau dependence cancelling.  <s> is useless (dominated by s=2) and s_max
    is an extreme-value/sample-size artifact.
    """
    v, h, _ = pooled(L, N)
    if v is None:
        return None
    c = h["mass"].astype(float)
    return mean_and_moment(v, c, 2) / mean_and_moment(v, c, 1)


def active_fraction(L, N):
    """Fraction of depositions that eliminate anything."""
    v, h, drops = pooled(L, N)
    if v is None or drops == 0:
        return None
    return h["mass"].sum() / drops


# ---------------------------------------------------------------------------
# Time-resolved: steady state + growth velocity
# ---------------------------------------------------------------------------
def moments_vs_t(L, N):
    """Pool outputs/moments over sims for (L, N).

    Returns (t_geom, <s>, <s^2>/<s>, active_frac, t_arith, mean_h).

    Two time axes on purpose.  The avalanche moments are stationary, so the
    geometric window centre is fine (and right for log plots).  But mean_h grows
    with t, and a window reports its *average* h; for a linear h(t) that average
    sits at the ARITHMETIC centre.  Using the geometric centre for h biases late
    wide windows and fakes a decreasing velocity.
    """
    files = sorted(glob.glob(f"{MOM}/L_{L}_N_{nfmt(N)}_steps_{STEPS}_sim_*.tsv"))
    if not files:
        return (None,) * 6
    acc = {}
    for f in files:
        for line in open(f):
            if line.startswith("#") or line.startswith("step_lo"):
                continue
            p = line.split()
            lo, hi = int(p[0]), int(p[1])
            a = acc.setdefault((lo, hi), [0, 0, 0, 0, 0.0, 0])
            a[0] += int(p[2]); a[1] += int(p[3]); a[2] += int(p[4]); a[3] += int(p[5])
            a[4] += float(p[7]); a[5] += 1
    keys = sorted(acc)
    tg, m1, m21, af, ta, mh = [], [], [], [], [], []
    for lo, hi in keys:
        drops, active, s1, s2, hsum, n = acc[(lo, hi)]
        ta.append((lo + hi) / 2.0)
        mh.append(hsum / n)
        if active == 0 or s1 == 0:
            tg.append(np.sqrt(max(lo, 1) * max(hi, 1))); m1.append(np.nan)
            m21.append(np.nan); af.append(0.0)
            continue
        tg.append(np.sqrt(max(lo, 1) * max(hi, 1)))
        m1.append(s1 / active); m21.append(s2 / s1); af.append(active / drops)
    return (np.array(tg), np.array(m1), np.array(m21), np.array(af),
            np.array(ta), np.array(mh))


# ---------------------------------------------------------------------------
# Mechanism: what variable drives the stretched exponential, and the sawtooth
# ---------------------------------------------------------------------------
def slope_resolved(L, N):
    """Pool outputs/slopeResolved for (L, N).

    Returns (m, P(m), <s|m>, w, P(w), <s|w>).
      m : local slope h[c+1]-h[c] at the deposition column, read pre-landing
      w : distinct columns the whole cascade eliminates in (its spatial extent)
    <s|m> is over ALL depositions (inactive ones contribute s=0), so it is the
    honest "does slope predict mass" answer.
    """
    files = sorted(glob.glob(f"{SLOPE}/L_{L}_N_{nfmt(N)}_steps_{STEPS}_sim_*.tsv"))
    if not files:
        return (None,) * 6
    macc, wacc = {}, {}
    for f in files:
        head, _, tail = open(f).read().partition("# extent")
        for line in head.strip().split("\n"):
            if line.startswith("#") or line.startswith("m\t"):
                continue
            p = [int(x) for x in line.split()]
            a = macc.setdefault(p[0], [0, 0, 0])
            a[0] += p[1]; a[1] += p[2]; a[2] += p[3]
        for line in tail.strip().split("\n"):
            if not line or line.startswith("w\t"):
                continue
            p = [int(x) for x in line.split()]
            a = wacc.setdefault(p[0], [0, 0])
            a[0] += p[1]; a[1] += p[2]
    mk = sorted(macc)
    md = np.array([macc[k][0] for k in mk], float)
    ms = np.array([macc[k][2] for k in mk], float)
    wk = sorted(wacc)
    wc = np.array([wacc[k][0] for k in wk], float)
    ws = np.array([wacc[k][1] for k in wk], float)
    return (np.array(mk, float), md / md.sum(), np.divide(ms, md, out=np.zeros_like(ms), where=md > 0),
            np.array(wk, float), wc / wc.sum(), np.divide(ws, wc, out=np.zeros_like(ws), where=wc > 0))


def composition(L, N):
    """Return (deposition probability, pile fraction) per species.

    Tests the sawtooth mechanism: a fractional species is deposited with
    probability f/N but, being rare, seldom finds a partner -- so if it is
    ENRICHED in the pile it is acting as a frozen defect.
    """
    files = sorted(glob.glob(f"{SLOPE}/L_{L}_N_{nfmt(N)}_steps_{STEPS}_sim_*.tsv"))
    if not files:
        return None, None
    tot = None
    for f in files:
        for line in open(f):
            if line.startswith("# pile_by_species="):
                pass
            if "pile_by_species=" in line:
                c = np.array([int(x) for x in line.split("pile_by_species=")[1].split(",")], float)
                tot = c if tot is None else tot + c
                break
    if tot is None:
        return None, None
    n_int = int(np.floor(N))
    frac = N - n_int
    w = [1.0] * n_int + ([frac] if frac > 0 else [])
    dep = np.array(w) / sum(w)
    return dep, tot / tot.sum()


def impurity_enrichment(L, N):
    """Pile fraction / deposition fraction for the fractional species.

    1.0 = no effect (integer N has no fractional species and returns None).
    >1 = the impurity accumulates: it is a frozen defect.
    """
    n_int = int(np.floor(N))
    if N - n_int <= 0:
        return None
    dep, pile = composition(L, N)
    if dep is None:
        return None
    return float(pile[-1] / dep[-1])


def weibull_slope(value, counts, lo, hi):
    """Slope k of the Weibull plot ln(-ln S) vs ln(value) over [lo, hi].

    P(x) Weibull  <=>  S(x)=exp(-(x/a)^k)  <=>  ln(-ln S) linear in ln x.
    Uses the survival function, so no binning and far less noise than the PDF.
    Returns (k, R^2).
    """
    p = counts.astype(float) / counts.sum()
    S = 1.0 - np.cumsum(p)
    m = (S > 1e-12) & (value >= lo) & (value <= hi) & (counts > 0)
    if m.sum() < 4:
        return np.nan, np.nan
    x = np.log(value[m].astype(float))
    y = np.log(-np.log(S[m]))
    k = np.polyfit(x, y, 1)[0]
    return float(k), float(np.corrcoef(x, y)[0, 1] ** 2)


def velocity(L, N, tmin_frac=0.35):
    """Growth velocity v = d<h>/dt from a late-time linear fit. The order parameter.

    Fit only t > tmin_frac*STEPS: h(t) is linear only once the surface has
    equilibrated, and at N_c it is a sublinear power law instead.
    """
    tg, m1, m21, af, ta, mh = moments_vs_t(L, N)
    if tg is None:
        return None
    m = ta > tmin_frac * STEPS
    if m.sum() < 3:
        return None
    return float(np.polyfit(ta[m], mh[m], 1)[0])
