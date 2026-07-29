"""Emit the argument lines for the criticalScaling sweep.

    python run_sweep.py | xargs -P 12 -L 1 ./criticalScaling > /dev/null

Box height H is the one parameter that must be right: too small and the pile
hits the ceiling, which caps it and **fakes an arrest** -- the exact artifact
this study is trying to detect.  H is therefore sized from a measured v(N) with
a 1.8x margin, and every run reports ceiling_hits (must be 0; the analysis
refuses contaminated runs).

Parallelism is bounded by memory, not cores: a job holds 5*H*L bytes.
"""
STEPS = 2000000

# Measured late-time velocity at L=4096; ~0 in the bound phase (N < N_c=5.0765).
VEL = {5.00: 0, 5.02: 0, 5.04: 0, 5.05: 0, 5.06: 0, 5.065: 0,
       5.070: 1e-4, 5.072: 2e-4, 5.075: 5e-4,
       5.078: 1e-3, 5.080: 1.5e-3, 5.085: 3e-3, 5.090: 4.5e-3}

# main N scan at fixed large L
MAIN_L, MAIN_NSIMS = 4096, 4
# finite-size check at one N: the arrest must be shown to survive L -> infinity
FS_N, FS_LS, FS_NSIMS = 5.075, [512, 1024, 2048, 4096, 8192], 8


def box_H(N, L):
    return int(max(300, VEL[N] * STEPS * 1.8 + 800))


def jobs():
    for N in VEL:
        for s in range(MAIN_NSIMS):
            yield MAIN_L, N, STEPS, s, box_H(N, MAIN_L)
    for L in FS_LS:
        if L == MAIN_L:
            continue
        for s in range(FS_NSIMS):
            # small L stays above threshold and grows, so it needs a TALLER box
            # than the same N at large L -- this is what the first run got wrong
            yield L, FS_N, STEPS, s, int(2500 * 4096 / L) + 500


if __name__ == "__main__":
    # biggest boxes first, so the memory peak happens while the queue is full
    for L, N, steps, s, H in sorted(jobs(), key=lambda j: -j[4] * j[0]):
        print(f"{L} {N} {steps} {s} {H}")
