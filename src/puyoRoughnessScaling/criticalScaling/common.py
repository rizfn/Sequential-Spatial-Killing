"""Loading and scaling analysis for the growth-arrest transition at N_c.

The object of study is the growth velocity v = d<h>/dt, which is the order
parameter: it is strictly positive for N > N_c and zero for N < N_c, where the
pile arrests at a finite thickness held up by the floor.

Everything here assumes the standard finite-size-scaling ansatz for a continuous
transition with correlation length xi ~ |N-N_c|^-nu and dynamic exponent z:

    <h>(t, N, L) = L^{a} F( t / L^z , (N-N_c) L^{1/nu} )

At N = N_c and t << L^z the system cannot know L, so <h> ~ t^beta with
beta = a/z; for t >> L^z it crosses over to linear growth at an L-dependent
velocity v(N_c, L) ~ L^{a-z} = L^{-theta/nu}.  Requiring the early-time regime to
be L-independent forces a = z*beta, hence the exponent relation

    theta/nu = z (1 - beta)

which is the main internal consistency check available to us: theta/nu comes from
v(N_c,L), beta from the t^beta regime, z from the crossover time, and the three
must close.
"""
import os
import re
import glob
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "outputs", "growth")
PLOTS = os.path.join(HERE, "plots")


def fname(L, N, steps, sim):
    return os.path.join(OUT, f"L_{L}_N_{N:.4f}_steps_{steps}_sim_{sim}.tsv")


def load(L, N, steps):
    """Average h(t), W(t) over all sims present.  Returns (t, h, W, nsims).

    Sims are averaged at the sample level, which is safe because the sample
    schedule is a deterministic function of `steps` alone -- every sim with the
    same `steps` samples at exactly the same times.
    """
    files = sorted(glob.glob(fname(L, N, steps, 0).replace("_sim_0", "_sim_*")))
    if not files:
        return None
    hs, ws, afs, ss = [], [], [], []
    t = None
    for f in files:
        # a run still in flight has a header but no rows yet
        if os.path.getsize(f) < 64:
            continue
        d = np.loadtxt(f, skiprows=2, ndmin=2)
        if d.size == 0:
            continue
        if t is None:
            t = d[:, 0]
        elif len(d) != len(t):
            continue
        hs.append(d[:, 1]); ws.append(d[:, 2]); afs.append(d[:, 3]); ss.append(d[:, 4])
    if not hs:
        return None
    return dict(t=t, h=np.mean(hs, axis=0), w=np.mean(ws, axis=0),
                active=np.mean(afs, axis=0), mean_s=np.mean(ss, axis=0),
                h_all=np.array(hs), nsims=len(hs))


def ceiling_hits(L, N, steps):
    tot = 0
    for f in glob.glob(fname(L, N, steps, 0).replace("_sim_0", "_sim_*")):
        with open(f) as fh:
            m = re.search(r"ceiling_hits=(\d+)", fh.readline())
            if m:
                tot += int(m.group(1))
    return tot


def local_slope(t, y, halfwidth=0.6):
    """d log y / d log t, from a sliding least-squares fit of half-width
    `halfwidth` decades (natural log).  Returns nan where the window is too thin.
    """
    m = (t > 0) & (y > 0)
    lt, ly = np.log(t[m]), np.log(y[m])
    out = np.full(len(lt), np.nan)
    for i in range(len(lt)):
        w = (lt > lt[i] - halfwidth) & (lt < lt[i] + halfwidth)
        if w.sum() > 3:
            out[i] = np.polyfit(lt[w], ly[w], 1)[0]
    return t[m], out


def velocity(L, N, steps, frac=0.5, err=False):
    """Late-time growth velocity, fitted over the last `frac` of the time axis.

    Fitting <h> vs t linearly over a late window measures v only if the run has
    actually crossed over to linear growth there; `linear_frac` reports whether
    it has, and callers MUST check it rather than trusting v blindly.
    """
    d = load(L, N, steps)
    if d is None:
        return None
    t, h = d["t"], d["h"]
    m = t >= t[-1] * frac
    if m.sum() < 4:
        return None
    v = np.polyfit(t[m], h[m], 1)[0]
    if not err:
        return v
    # error across sims: refit each sim independently
    vs = [np.polyfit(t[m], hh[m], 1)[0] for hh in d["h_all"]]
    return v, np.std(vs) / np.sqrt(len(vs))


def linear_frac(L, N, steps, frac=0.5):
    """Local log-slope of <h>(t) averaged over the late window used by
    `velocity`.  ~1 means the run really is in the linear regime and v is
    meaningful; well below 1 means it is still in the critical t^beta regime and
    the "velocity" is a transient artifact.
    """
    d = load(L, N, steps)
    if d is None:
        return None
    t, sl = local_slope(d["t"], d["h"])
    m = t >= t[-1] * frac
    return np.nanmean(sl[m])


def collapse_cost(xs, ys):
    """Spread of several curves about their common envelope, in log-log.

    Curves are compared on the overlap of their x-ranges by interpolating each
    onto a shared grid and taking the RMS deviation from the pointwise mean.
    Scale-free, so costs are comparable across different (a, z) trials.
    """
    lo = max(x.min() for x in xs)
    hi = min(x.max() for x in xs)
    if not (hi > lo):
        return np.inf
    grid = np.linspace(np.log(lo), np.log(hi), 60)
    curves = []
    for x, y in zip(xs, ys):
        o = np.argsort(x)
        curves.append(np.interp(grid, np.log(x[o]), np.log(y[o])))
    curves = np.array(curves)
    return float(np.sqrt(np.mean((curves - curves.mean(axis=0)) ** 2)))
