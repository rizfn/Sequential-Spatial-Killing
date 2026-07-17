"""Fit-free reading of the stretch exponent b.

For P(s) = C s^-tau exp[-(s/s0)^b], the local log-slope is

    tau_eff(s) = -dlnP/dlns = tau + b*(s/s0)^b

so the EXCESS local slope above its small-s plateau grows as a pure power of s
with exponent b.  Plot ln(tau_eff - tau) vs ln s -> slope = b, with no fitting
of the distribution itself.  This is independent of the tau<->b ridge, because
it reads b off the *shape* of the local slope, not off a global normalization.
"""
import numpy as np


def tau_eff(value, counts, nbins=26, smin=8):
    """Local log-slope of P(s) from log-binned data (centred differences)."""
    m = counts > 0
    v, c = value[m].astype(float), counts[m].astype(float)
    edges = np.unique(np.round(np.logspace(np.log10(smin), np.log10(v.max()+1), nbins)).astype(int))
    xs, ys, ns = [], [], []
    tot = c.sum()
    for a, b in zip(edges[:-1], edges[1:]):
        sel = (v >= a) & (v < b)
        if not sel.any():
            continue
        cc = c[sel].sum()
        if cc < 50:            # need real statistics for a slope
            continue
        xs.append(np.exp(np.log(v[sel]).mean()))
        ys.append(cc / tot / (b - a))
        ns.append(cc)
    x, y, ns = np.log(np.array(xs)), np.log(np.array(ys)), np.array(ns)
    # centred difference of ln P wrt ln s
    te = -np.gradient(y, x)
    return np.exp(x), te, ns
