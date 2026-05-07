import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

L = 128
STEPS = 2048*2
TRANSIENT_FRACTION = 0.5

BASE_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = BASE_DIR / "outputs" / "massVsTime2D"
PLOTS_DIR = BASE_DIR / "plots"


def load_mass_timeseries(l, steps):
    pattern = f"L_{l}_N_*_steps_{steps}_sim_*.tsv"
    rx = re.compile(r"L_(\d+)_N_([0-9eE+\-.]+)_steps_(\d+)_sim_(\d+)\.tsv$")

    series = []
    for path in OUTPUT_DIR.glob(pattern):
        m = rx.search(path.name)
        if not m:
            continue

        n_value = float(m.group(2))
        sim_no = int(m.group(4))
        step, mass, _height = np.loadtxt(path, delimiter="\t", comments="#", unpack=True)

        fit_start = int(len(step) * TRANSIENT_FRACTION)
        fit_start = min(max(fit_start, 0), len(step) - 2)
        slope = float(np.polyfit(step[fit_start:], mass[fit_start:], 1)[0])

        series.append((n_value, sim_no, step, mass, slope))

    # Sort by (N, sim_no) for consistent ordering
    series.sort(key=lambda x: (x[0], x[1]))
    return series


def plot_mass_vs_time(series, l, steps):
    fig, ax = plt.subplots(figsize=(12, 7))

    # Group by N value for coloring
    by_n = {}
    for n_value, sim_no, step, mass, _ in series:
        if n_value not in by_n:
            by_n[n_value] = []
        by_n[n_value].append((sim_no, step, mass))

    cmap = plt.get_cmap("rainbow", max(len(by_n), 2))
    for color_idx, (n_value, runs) in enumerate(sorted(by_n.items())):
        color = cmap(color_idx)
        for sim_no, step, mass in runs:
            ax.plot(step, mass, color=color, alpha=0.6, linewidth=1.0, label=f"N={n_value:g}, sim={sim_no}")

    ylim = ax.get_ylim()
    ax.plot([0, steps], [l, l * (steps + 1)], color="grey", linestyle="--", linewidth=2, label="Total mass added")
    ax.set_ylim(ylim)
    
    ax.set_xlabel("Time")
    ax.set_ylabel("Mass")
    ax.grid(alpha=0.3)
    ax.legend(ncols=2, fontsize=8)

    PLOTS_DIR.mkdir(parents=True, exist_ok=True)
    out = PLOTS_DIR / f"massVsTime_L_{l}_steps_{steps}.png"
    plt.tight_layout()
    plt.savefig(out, dpi=300)
    print(f"Saved: {out}")


def plot_drift_vs_inverse_n(series, l, steps):
    # Collect individual slopes by N value
    by_n = {}
    for n_value, sim_no, _step, _mass, slope in series:
        if n_value > 0:
            if n_value not in by_n:
                by_n[n_value] = []
            by_n[n_value].append((sim_no, slope))

    fig, ax = plt.subplots(figsize=(9, 6))

    cmap = plt.get_cmap("tab10", max(len(by_n), 2))
    for color_idx, (n_value, runs) in enumerate(sorted(by_n.items())):
        color = cmap(color_idx % 10)
        inv_n = 1.0 / n_value
        slopes = [slope for _, slope in runs]

        # Plot individual points
        for sim_no, slope in runs:
            ax.scatter([inv_n], [slope], color=color, s=50, alpha=0.6)

        # Plot mean with error bar
        mean_slope = np.mean(slopes)
        std_slope = np.std(slopes)
        ax.errorbar([inv_n], [mean_slope], yerr=std_slope, fmt='o', capsize=5,
                    color=color, markersize=8, linewidth=2, label=f"N={n_value:g}")

    ax.set_xlabel("Inverse N (1/N)")
    ax.set_ylabel("Drift (dM/dt)")
    ax.set_title("Drift vs Inverse N (individual runs + mean ± std)")
    ax.grid(alpha=0.3)
    ax.legend()

    PLOTS_DIR.mkdir(parents=True, exist_ok=True)
    out = PLOTS_DIR / f"driftVsInverseN_L_{l}_steps_{steps}.png"
    plt.tight_layout()
    plt.savefig(out, dpi=300)
    print(f"Saved: {out}")
    

def main():
    series = load_mass_timeseries(L, STEPS)
    if not series:
        raise FileNotFoundError(
            f"No files found in {OUTPUT_DIR} matching L={L}, steps={STEPS}."
        )

    plot_mass_vs_time(series, L, STEPS)
    plot_drift_vs_inverse_n(series, L, STEPS)


if __name__ == "__main__":
    main()
