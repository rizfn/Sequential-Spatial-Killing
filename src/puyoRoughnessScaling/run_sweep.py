"""Emit the sweep job list (one line per sim: "L N steps sim H rec").

Single source of truth = common.FSWEEP / LEXTRA / NSIMS.  Piped into xargs for
execution (see the shell command that launches the sweep).
"""
from common import FSWEEP, LEXTRA, NSIMS, FIXED_L


def main():
    lines = []
    for N, cfg in FSWEEP.items():
        rec = max(1, cfg["steps"] // 2500)
        for L in cfg["Ls"]:
            for sim in range(NSIMS):
                lines.append(f"{L} {N} {cfg['steps']} {sim} {cfg['H']} {rec}")
    for N, cfg in LEXTRA.items():
        rec = max(1, cfg["steps"] // 2500)
        for sim in range(NSIMS):
            lines.append(f"{FIXED_L} {N} {cfg['steps']} {sim} {cfg['H']} {rec}")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
