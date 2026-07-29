"""Loading + candidate variables for v(N) across the fractional-species scan.

The question: v is (nearly) linear in 1/N at INTEGER N, but not at fractional N.
What variable restores linearity?

The candidates fall into two families, and the distinction is the whole point:

  functions of the DEPOSITION weights p alone -- 1/N, sum p^2, and every Renyi /
  Hill diversity of any order q.  These are what "is it information entropy?"
  is asking about.

  quantities involving the PILE composition rho -- sum p.rho, sum rho^2.  rho is
  a dynamical OUTPUT, not a control parameter, because the rare fractional
  species self-poisons (it can only be removed in pairs, so it accumulates).

At integer N symmetry forces rho = p = uniform and every candidate collapses to
1/N, which is exactly why only the interiors of the intervals can discriminate.
"""
import os
import re
import glob

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "outputs", "velocity")
PLOTS = os.path.join(HERE, "plots")

L = 1024
STEPS = 50000
NS = [round(x, 3) for x in np.arange(5.5, 10.0 + 1e-9, 0.1)]


def _files(N):
    return sorted(glob.glob(os.path.join(OUT, f"L_{L}_N_{N:.4f}_steps_{STEPS}_sim_*.tsv")))


HEADER_LINES = 5   # ceiling/pile / elim / weights / column names


def load(N):
    """Mean velocity, its standard error, the pile composition rho, and the
    per-species removal rate r.

    v is the slope of a straight-line fit of mass/L against t over the recorded
    (post-warmup) window, fitted per sim and then averaged, so the error bar is
    the spread across independent sims rather than of a single fit.

    r_i is eliminations of species i per deposition, over the same post-warmup
    window.  It is what closes the exact balance p_i - r_i = v*rho_i, so it is
    measured rather than guessed (see closure.py).
    """
    vs, rho, r, ceil = [], None, None, 0
    for f in _files(N):
        head = open(f).readlines()[:4]
        ceil += int(re.search(r"ceiling_hits=(\d+)", head[0]).group(1))
        d = np.loadtxt(f, skiprows=HEADER_LINES, ndmin=2)
        vs.append(np.polyfit(d[:, 0], d[:, 1], 1)[0])
        sp = np.array([float(x) for x in
                       re.search(r"pile_by_species=([\d,]+)", head[1]).group(1).split(",")])
        el = np.array([float(x) for x in
                       re.search(r"elim_by_species=([\d,]+)", head[2]).group(1).split(",")])
        nd = float(re.search(r"elim_drops=(\d+)", head[2]).group(1))
        rho = sp / sp.sum() if rho is None else rho + sp / sp.sum()
        r = el / nd if r is None else r + el / nd
    if not vs:
        return None
    return dict(v=float(np.mean(vs)), err=float(np.std(vs) / np.sqrt(len(vs))),
                rho=rho / len(vs), r=r / len(vs), nsims=len(vs), ceiling_hits=ceil)


def weights(N):
    """Deposition probabilities: floor(N) species of weight 1 plus a fractional
    one of weight f, all normalised by N.
    """
    n = int(np.floor(N))
    f = N - n
    w = np.array([1.0] * n + ([f] if f > 0 else []))
    return w / w.sum()


def hill(p, q):
    """Hill number of order q = exp(Renyi entropy of order q): the 'effective
    number of species'.  q=1 is Shannon, q=2 the collision/Simpson index.
    """
    p = p[p > 0]
    if abs(q - 1) < 1e-9:
        return float(np.exp(-(p * np.log(p)).sum()))
    return float((p ** q).sum() ** (1 / (1 - q)))


def candidates(N, rho):
    """Every candidate x for which we ask 'is v linear in x?'."""
    p = weights(N)
    return {
        "1/N": 1 / N,
        "sum p^2": float((p ** 2).sum()),
        "1/Hill_1": 1 / hill(p, 1.0),
        "sum p.rho": float((p * rho).sum()),
        "sum rho^2": float((rho ** 2).sum()),
    }


def linearity(xs, ys):
    """RMS residual of a straight-line fit -- the figure of merit throughout.

    Compared against the noise floor (the mean standard error on v) and against
    the integer-only 1/N fit, which is the best linearity the model itself
    offers and therefore the target, not zero.
    """
    p = np.polyfit(xs, ys, 1)
    r = ys - np.polyval(p, xs)
    return float(np.sqrt(np.mean(r ** 2))), float(np.max(np.abs(r))), p


def table():
    """N -> dict(v, err, rho, candidate values). Skips ceiling-contaminated N."""
    out = {}
    for N in NS:
        d = load(N)
        if d is None or d["ceiling_hits"] > 0:
            continue
        d["cands"] = candidates(N, d["rho"])
        out[N] = d
    return out
