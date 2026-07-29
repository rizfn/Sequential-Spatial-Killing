"""Can sum_i p_i rho_i be computed from the drop probabilities p alone?

Exact steady-state balance (composition fixed while the pile grows at v per
deposition), per species:

    p_i - r_i = v rho_i        added - removed = what the growing pile carries

Summing gives v = 1 - sum_i r_i, the mass balance of ../avalancheScaling.  Exact,
but not closed: r_i needs the dynamics.

Closure: nothing happens until a block is DROPPED, so species i is removed when
species i is deposited (rate p_i) and lands next to a match (prob ~ rho_i):

    r_i = k p_i rho_i

Measured directly (R^2=0.998, 2.6% rel RMS) -- it beats the equilibrium-flavoured
guess r_i = k rho_i^2 (R^2=0.965, 27%), which is wrong here because pile blocks
do not move on their own.
That closes it -- rho, and hence sum p_i rho_i and v, follow from p alone with a
single constant k:

    p_i - k p_i rho_i = v rho_i     =>     rho_i = p_i / (v + k p_i)

and sum_i rho_i = 1 then fixes v.

VERDICT (readme §5): rho is reproduced to 1-2% -- the pile probabilities ARE
derivable from the drop probabilities.  But v is NOT: v = 1 - sum_i r_i is a
near-total cancellation (sum r = 0.94 at N=6), so the closure's 2.6% error in r
becomes a 42% error in v, diverging as v -> 0 near N_c.
"""
import numpy as np
from scipy.optimize import brentq, minimize_scalar

import common


def rho_of(v, k, p):
    """rho_i = p_i / (v + k p_i), from r_i = k p_i rho_i + the exact balance."""
    return p / (v + k * p)


def predict_v(N, k):
    """v from p alone. sum(rho) falls with v; if it is already <1 at v=0 there is
    no growing solution and the closure says the pile is arrested."""
    p = common.weights(N)
    if rho_of(0.0, k, p).sum() <= 1.0:
        return 0.0
    return brentq(lambda v: rho_of(v, k, p).sum() - 1.0, 0.0, 5.0)


def fit_k(T):
    ns = sorted(T)
    cost = lambda k: np.sqrt(np.mean([(predict_v(N, k) - T[N]["v"]) ** 2 for N in ns]))
    r = minimize_scalar(cost, bounds=(4.0, 7.0), method="bounded")
    return r.x, r.fun


if __name__ == "__main__":
    T = common.table()
    k, rms = fit_k(T)
    print(f"fitted k = {k:.3f}   RMS = {rms:.5f}  ({rms/2e-4:.0f}x noise floor)")
    print(f"closure N_c (integer family, where n_species = k) = {k:.3f}")
    print(f"  measured N_c = 5.0765 (../criticalScaling)  -> {abs(k-5.0765)/5.0765:.1%} off")
    print(f"\n{'N':>6} {'v meas':>9} {'v pred':>9} {'diff':>9}")
    for N in [5.5, 6.0, 6.1, 6.5, 7.0, 8.0, 9.0, 10.0]:
        pv = predict_v(N, k)
        print(f"{N:6.1f} {T[N]['v']:9.5f} {pv:9.5f} {pv-T[N]['v']:+9.5f}")
