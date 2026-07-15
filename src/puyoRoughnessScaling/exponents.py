"""Finite-size scaling exponents for the Puyo deposition interface.

For each N, measure INDEPENDENTLY (bootstrap over sims):
  alpha : W_sat ~ L^alpha        (over the L that have SATURATED)
  beta  : W ~ t^beta             (early growth on the largest saturated L)
  z     : t_x  ~ L^z             (crossover time, W reaches 0.75 W_sat)
and test the Family-Vicsek identity  beta * z = alpha.

A per-(N,L) saturation check (sim-averaged late/mid W ratio) filters out sizes
that have not saturated, so they don't corrupt the alpha/z fits.

Run directly to print the table; import measure_all()/measure_fv() for figures.
"""
import numpy as np
from common import stack, FSWEEP

_rng = np.random.default_rng(1)
CFRAC = 0.75          # crossover threshold as fraction of W_sat
SAT_THRESH = 1.06     # late/mid W ratio below this => saturated


def _load(N):
    return {L: stack(L, N) for L in FSWEEP[N]["Ls"]}


def sat_ratio(t, Wm):
    T = t.max()
    return Wm[t >= 0.8 * T].mean() / Wm[(t >= 0.4 * T) & (t <= 0.6 * T)].mean()


def saturated_Ls(stacks):
    out = []
    ratios = {}
    for L, (t, Warr) in stacks.items():
        r = sat_ratio(t, Warr.mean(0))
        ratios[L] = r
        if r < SAT_THRESH:
            out.append(L)
    return sorted(out), ratios


def _bootstrap(N, stacks, satLs, nboot=300):
    gw = FSWEEP[N]["gw"]
    Lbeta = satLs[-1] if satLs else max(stacks)   # largest saturated L for beta
    La = np.array(satLs, float)
    a_s, b_s, z_s = [], [], []
    for _ in range(nboot):
        Wsat, tx = [], []
        for L in satLs:
            t, Warr = stacks[L]
            Wm = Warr[_rng.integers(0, Warr.shape[0], Warr.shape[0])].mean(0)
            ws = Wm[t >= 0.6 * t.max()].mean()
            Wsat.append(ws)
            k = np.argmax(Wm >= CFRAC * ws)
            tx.append(t[k] if Wm[k] >= CFRAC * ws else np.nan)
        # beta on largest saturated L
        t, Warr = stacks[Lbeta]
        Wm = Warr[_rng.integers(0, Warr.shape[0], Warr.shape[0])].mean(0)
        m = (t >= gw[0]) & (t <= gw[1]) & (Wm > 0)
        b_s.append(np.polyfit(np.log(t[m]), np.log(Wm[m]), 1)[0])
        Wsat = np.array(Wsat); tx = np.array(tx)
        a_s.append(np.polyfit(np.log(La), np.log(Wsat), 1)[0])
        mz = np.isfinite(tx) & (tx > 0)
        z_s.append(np.polyfit(np.log(La[mz]), np.log(tx[mz]), 1)[0] if mz.sum() >= 2 else np.nan)
    return np.array(a_s), np.array(b_s), np.array(z_s)


def measure_all(nboot=300):
    res = {}
    for N in FSWEEP:
        stacks = _load(N)
        satLs, ratios = saturated_Ls(stacks)
        Wsat = np.array([stacks[L][1].mean(0)[stacks[L][0] >= 0.6 * stacks[L][0].max()].mean()
                         for L in FSWEEP[N]["Ls"]])
        if len(satLs) >= 3:
            a_s, b_s, z_s = _bootstrap(N, stacks, satLs, nboot)
        else:
            a_s = b_s = z_s = np.array([np.nan])
        res[N] = dict(alpha=a_s, beta=b_s, z=z_s, Wsat=Wsat, ratios=ratios,
                      Ls=FSWEEP[N]["Ls"], satLs=satLs)
    return res


# --- Family-Vicsek inputs with SYSTEMATIC (window/threshold) uncertainty -------
_Z_THRESH = [0.60, 0.70, 0.80, 0.85]


def _beta_windows(N):
    lo, hi = FSWEEP[N]["gw"]
    return [(lo, hi), (lo, int(hi * 0.6)), (int(lo * 1.5) + 1, hi), (lo, int(hi * 1.5))]


def measure_fv(N, nboot=600):
    stacks = _load(N)
    satLs, _ = saturated_Ls(stacks)
    if len(satLs) < 3:
        return np.array([np.nan]), np.array([np.nan])
    Lbeta = satLs[-1]
    La = np.array(satLs, float)
    wins = _beta_windows(N)
    a_s, bz_s = [], []
    for _ in range(nboot):
        lo, hi = wins[_rng.integers(len(wins))]
        cf = _Z_THRESH[_rng.integers(len(_Z_THRESH))]
        Wsat, tx = [], []
        for L in satLs:
            t, Warr = stacks[L]
            Wm = Warr[_rng.integers(0, Warr.shape[0], Warr.shape[0])].mean(0)
            ws = Wm[t >= 0.6 * t.max()].mean()
            Wsat.append(ws)
            k = np.argmax(Wm >= cf * ws)
            tx.append(t[k] if Wm[k] >= cf * ws else np.nan)
        t, Warr = stacks[Lbeta]
        Wm = Warr[_rng.integers(0, Warr.shape[0], Warr.shape[0])].mean(0)
        m = (t >= lo) & (t <= hi) & (Wm > 0)
        beta = np.polyfit(np.log(t[m]), np.log(Wm[m]), 1)[0]
        Wsat = np.array(Wsat); tx = np.array(tx)
        alpha = np.polyfit(np.log(La), np.log(Wsat), 1)[0]
        mz = np.isfinite(tx) & (tx > 0)
        z = np.polyfit(np.log(La[mz]), np.log(tx[mz]), 1)[0] if mz.sum() >= 2 else np.nan
        a_s.append(alpha); bz_s.append(beta * z)
    return np.array(a_s), np.array(bz_s)


def _pm(x):
    return f"{np.nanmean(x):.3f}±{np.nanstd(x):.3f}"


if __name__ == "__main__":
    res = measure_all()
    hdr = f"{'N':>3} | {'alpha':>11} | {'beta':>11} | {'z':>11} | {'beta*z':>11} | saturated L / all"
    print(hdr); print("-" * len(hdr))
    for N, r in res.items():
        bz = r["beta"] * r["z"]
        rat = " ".join(f"{L}:{r['ratios'][L]:.2f}" for L in r["Ls"])
        print(f"{N:>3} | {_pm(r['alpha']):>11} | {_pm(r['beta']):>11} | {_pm(r['z']):>11} | "
              f"{_pm(bz):>11} | sat={r['satLs']}  [{rat}]")
