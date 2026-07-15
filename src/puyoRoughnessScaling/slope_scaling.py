"""Explore how the steady-state slope distribution P(m) scales with N.

The current figure rescales by N (assuming the Laplace scale lambda ~ N). Test
that assumption: measure lambda(N) directly and search for the scaling that
actually collapses the family.
"""
import glob
import numpy as np
from common import parse_slopes, pattern, steps_of, FIXED_L, all_N_at_fixed_L


def steady_slopes(N):
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
    return np.concatenate(out)


def pdf(m, scale):
    """Return (x, y) = (m/scale, scale*P(m)) so the area is preserved."""
    lo, hi = np.floor(m.min()) - 0.5, np.ceil(m.max()) + 1.5
    bins = np.arange(lo, hi, 1.0)
    p, e = np.histogram(m, bins=bins, density=True)
    c = 0.5 * (e[:-1] + e[1:])
    g = p > 0
    return c[g] / scale, p[g] * scale


def collapse_cost(curves, xmin, xmax):
    """Mean variance of log y across curves on a common x grid."""
    grid = np.linspace(xmin, xmax, 120)
    ys = []
    for x, y in curves:
        o = np.argsort(x)
        ys.append(np.interp(grid, x[o], np.log(y[o]), left=np.nan, right=np.nan))
    ys = np.array(ys)
    v = np.nanvar(ys, axis=0)
    return np.nanmean(v[np.isfinite(v)])


Ns = all_N_at_fixed_L()
print(f"{'N':>3} {'<|m|>=lam':>9} {'std':>7} {'P(0)':>7} {'std/lam':>7} {'exKurt':>7}  lam/N")
data = {}
for N in Ns:
    m = steady_slopes(N)
    lam = np.abs(m).mean()
    sd = m.std()
    p0 = np.mean(np.abs(m) < 0.5) / 1.0
    exk = ((m**4).mean() / sd**4) - 3
    data[N] = m
    print(f"{N:>3} {lam:9.2f} {sd:7.2f} {p0:7.3f} {sd/lam:7.3f} {exk:7.2f}  {lam/N:.3f}")

# lambda(N) power law
Na = np.array(Ns, float)
lam = np.array([np.abs(data[N]).mean() for N in Ns])
a, b = np.polyfit(np.log(Na), np.log(lam), 1)
print(f"\nlambda ~ N^{a:.3f}  (N^1 would justify the m/N collapse)")

# collapse costs for different scale choices (evaluate over the Laplace-tail region)
scales = {
    "m/N (current)": {N: N for N in Ns},
    "m/<|m|>": {N: np.abs(data[N]).mean() for N in Ns},
    "m/std": {N: data[N].std() for N in Ns},
    f"m/N^{a:.2f}": {N: Na_ ** a for N, Na_ in zip(Ns, Na)},
}
print(f"\n{'scaling':>16} : collapse cost (lower=better), core [-1.5,1.5] / tail [1.5,4]")
for name, sc in scales.items():
    curves = [pdf(data[N], sc[N]) for N in Ns]
    core = collapse_cost(curves, -1.5, 1.5)
    tail = collapse_cost([(np.abs(x), y) for x, y in curves], 1.5, 4.0)
    print(f"{name:>16} : core={core:.4f}  tail={tail:.4f}")
