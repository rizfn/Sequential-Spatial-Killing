"""Emit argument lines for the v(N) sweep.

    python run_sweep.py | xargs -P 16 -L 1 ./velocity > /dev/null

N is sampled finely THROUGH the integer intervals, because the whole question is
what happens between integers: at integer N the species are uniform and every
candidate variable coincides, so only the interiors of the intervals can
discriminate them.

Box height must cover v*steps with margin; a ceiling hit caps the pile and
corrupts v.  Every run reports ceiling_hits (must be 0).
"""
import numpy as np

L = 1024
STEPS = 50000
WARMUP = 10000
NSIMS = 16

# fine scan through the low intervals, where v is small and the curvature lives
NS = [round(x, 3) for x in np.arange(5.5, 10.0 + 1e-9, 0.1)]


def v_upper(N):
    """Generous upper bound on v, for sizing the box only."""
    return max(0.02, 1.0 - 4.9 / N)


def box_H(N):
    return int(v_upper(N) * STEPS * 1.35 + 500)


if __name__ == "__main__":
    for L_, N, s in sorted(((L, N, s) for N in NS for s in range(NSIMS)),
                           key=lambda j: -box_H(j[1])):
        print(f"{L_} {N} {STEPS} {s} {box_H(N)} {WARMUP}")
