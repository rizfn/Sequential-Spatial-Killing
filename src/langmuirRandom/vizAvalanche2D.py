import numpy as np
import matplotlib.pyplot as plt
from glob import glob
from pathlib import Path
import re

def load_individual_runs(data_glob, file_pattern):
    files = sorted(glob(data_glob))
    if not files:
        raise OSError(f"No files found for pattern: {data_glob}")

    runs = []
    for file_path in files:
        name = Path(file_path).name
        match = file_pattern.match(name)
        if not match:
            continue

        n_value = float(match.group("N"))
        sim_no = int(match.group("simNo"))

        try:
            data = np.loadtxt(file_path, delimiter="\t", comments="#")
        except ValueError:
            continue

        if data.size == 0:
            continue

        if data.ndim == 1:
            data = data.reshape(1, -1)

        if data.shape[1] < 3:
            continue

        runs.append({
            "N": n_value,
            "sim_no": sim_no,
            "avalanches": data[:, 1],
            "eliminated": data[:, 2],
            "file": file_path
        })

    # Sort by (N, sim_no) for consistent ordering
    runs.sort(key=lambda x: (x["N"], x["sim_no"]))
    return runs


def main():
    data_glob = "src/langmuirRandom/outputs/avalanche2D/*.tsv"
    file_pattern = re.compile(r"^L_(?P<L>\d+)_N_(?P<N>[^_]+)_(?P<simNo>\d+)\.tsv$")
    runs = load_individual_runs(data_glob, file_pattern)
    
    print(f"Loaded {len(runs)} individual runs")
    n_values = sorted(set(r["N"] for r in runs))
    print(f"Found N values: {', '.join(str(n) for n in n_values)}")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # Determine global bin range
    all_avalanches = np.concatenate([r["avalanches"] for r in runs])
    all_eliminated = np.concatenate([r["eliminated"] for r in runs])

    avalanche_bins = np.geomspace(1.0, max(2.0, float(np.max(all_avalanches))), 30)
    eliminated_bins = np.geomspace(1.0, max(2.0, float(np.max(all_eliminated))), 30)

    # Color by N value
    cmap = plt.get_cmap("tab10", max(len(n_values), 2))
    n_to_color = {n: cmap(i % 10) for i, n in enumerate(n_values)}

    # Plot each run individually
    for run in runs:
        n_value = run["N"]
        sim_no = run["sim_no"]
        avalanches = run["avalanches"]
        eliminated = run["eliminated"]
        color = n_to_color[n_value]

        avalanche_hist, _ = np.histogram(avalanches, bins=avalanche_bins)
        eliminated_hist, _ = np.histogram(eliminated, bins=eliminated_bins)

        avalanche_density = avalanche_hist / np.diff(avalanche_bins)
        eliminated_density = eliminated_hist / np.diff(eliminated_bins)

        label = f"N={n_value:g}, sim={sim_no}"
        ax1.plot(
            avalanche_bins[:-1][avalanche_density > 0],
            avalanche_density[avalanche_density > 0],
            marker="x",
            linewidth=0.8,
            alpha=0.7,
            label=label,
            color=color,
        )

        ax2.plot(
            eliminated_bins[:-1][eliminated_density > 0],
            eliminated_density[eliminated_density > 0],
            marker="x",
            linewidth=0.8,
            alpha=0.7,
            label=label,
            color=color,
        )

    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlabel("Avalanche size")
    ax1.set_ylabel("Frequency density")
    ax1.grid(alpha=0.3)
    ax1.legend(ncols=2, fontsize=8)

    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlabel("Total eliminated")
    ax2.set_ylabel("Frequency density")
    ax2.grid(alpha=0.3)
    ax2.legend(ncols=2, fontsize=8)

    out_path = Path("src/langmuirRandom/plots/avalanche/simple_combined_distributions.png")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()