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

    grouped = {}
    for path in OUTPUT_DIR.glob(pattern):
        m = rx.search(path.name)
        if not m:
            continue

        n_value = float(m.group(2))
        sim_no = int(m.group(4))
        step, mass, _height = np.loadtxt(path, delimiter="\t", skiprows=1, unpack=True)
        grouped.setdefault(n_value, []).append((sim_no, step, mass))

    series = []
    for n_value in sorted(grouped.keys()):
        runs = sorted(grouped[n_value], key=lambda x: x[0])
        min_len = min(len(step) for _, step, _ in runs)
        ref_step = runs[0][1][:min_len]
        masses = np.vstack([mass[:min_len] for _, _, mass in runs])

        mean_mass = masses.mean(axis=0)
        std_mass = masses.std(axis=0)

        fit_start = int(min_len * TRANSIENT_FRACTION)
        fit_start = min(max(fit_start, 0), min_len - 2)
        slopes = np.array(
            [np.polyfit(step[fit_start:min_len], mass[fit_start:min_len], 1)[0] for _, step, mass in runs]
        )
        mean_slope = float(np.mean(slopes))
        std_slope = float(np.std(slopes))

        series.append((n_value, ref_step, mean_mass, std_mass, mean_slope, std_slope))

    return series


def plot_mass_vs_time(series, l, steps):
    fig, ax = plt.subplots(figsize=(9, 6))

    cmap = plt.get_cmap("rainbow", max(len(series), 2))
    for i, (n_value, step, mean_mass, std_mass, _mean_slope, _std_slope) in enumerate(series):
        color = cmap(i)
        ax.plot(step, mean_mass, color=color, label=f"N={n_value:g}")
        ax.fill_between(step, mean_mass - std_mass, mean_mass + std_mass, color=color, alpha=0.18)

    ylim = ax.get_ylim()
    ax.plot([0, steps], [l, l * (steps + 1)], color="grey", linestyle="--", label="Total mass added")
    ax.set_ylim(ylim)
    
    ax.set_xlabel("Time")
    ax.set_ylabel("Mass")
    ax.grid()
    ax.legend(ncols=2)

    PLOTS_DIR.mkdir(parents=True, exist_ok=True)
    out = PLOTS_DIR / f"massVsTime_L_{l}_steps_{steps}.png"
    plt.tight_layout()
    plt.savefig(out, dpi=300)


def plot_drift_vs_inverse_n(series, l, steps):
    inv_n_values = []
    drifts = []
    drift_errs = []

    for n_value, _step, _mean_mass, _std_mass, mean_slope, std_slope in series:
        if n_value > 0:
            inv_n_values.append(1.0 / n_value)
            drifts.append(mean_slope)
            drift_errs.append(std_slope)

    fig, ax = plt.subplots(figsize=(7, 6))
    ax.errorbar(inv_n_values, drifts, yerr=drift_errs, fmt="o-", capsize=4)
    ax.set_xlabel("Inverse N (1/N)")
    ax.set_ylabel("Drift (dM/dt)")
    ax.set_title("Drift vs Inverse N")
    ax.grid()

    PLOTS_DIR.mkdir(parents=True, exist_ok=True)
    out = PLOTS_DIR / f"driftVsInverseN_L_{l}_steps_{steps}.png"
    plt.tight_layout()
    plt.savefig(out, dpi=300)
    

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
