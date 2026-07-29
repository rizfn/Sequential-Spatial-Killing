"""Emit the sweep job list, one line per sim: "L N steps sim H warmup".

Single source of truth = common.jobs().  Piped into xargs for parallel
execution:

    python run_sweep.py | xargs -P 12 -L 1 ./avalancheDist > /dev/null

-P is bounded by memory, not cores: a job holds ~5 bytes * H * L (uint8 lattice
+ int32 BFS stamps), which is ~203 MB for the L=4096 jobs.
"""
from common import jobs


def main():
    print("\n".join(" ".join(str(x) for x in j) for j in jobs()))


if __name__ == "__main__":
    main()
