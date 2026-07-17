"""Job list for the integer-N avalanche family.

Fractional N is a different system (the partial species is a self-poisoning
impurity -- see readme's sawtooth section), so this family is INTEGER N only.
Fixed L=1024 (P(s) is L-converged by ~512; avalanches have w_0 ~ 4.8 columns),
and buys tail statistics with sims instead of L: the tail shape is the question.

H is sized from v(N) measured at L=256 (probe): v rises from 0.38 at N=9 to
0.64 at N=16, so the r=0.55 of common.box_H would hit the ceiling above N~9.
"""
import sys

L = 1024
STEPS = 32768
WARMUP = STEPS // 4
NSIMS = 128

# measured v(N), L=256 probe -> box height with 25% margin
V = {6: 0.061, 7: 0.216, 8: 0.310, 9: 0.383, 10: 0.451, 12: 0.531, 16: 0.641}
NS = [6, 7, 8, 9, 10, 12, 16]


def box_H(N):
    # keep 6/7/8 identical to common.box_H so existing sims 0-15 pool cleanly
    if N == 6:
        return int(0.30 * STEPS) + 512
    if N in (7, 8):
        return int(0.55 * STEPS) + 512
    return int(1.25 * V[N] * STEPS) + 512


if __name__ == "__main__":
    jobs = [(L, N, STEPS, sim, box_H(N), WARMUP) for N in NS for sim in range(NSIMS)]
    jobs.sort(key=lambda j: -j[4])          # biggest box first: bounds peak memory
    for j in jobs:
        print(" ".join(str(x) for x in j))
