import numpy as np
import matplotlib.pyplot as plt
from glob import glob
from pathlib import Path
import re

def load_grouped_by_n(data_glob, file_pattern):
    files = sorted(glob(data_glob))
    if not files:
        raise OSError(f"No files found for pattern: {data_glob}")

    grouped = {}
    for file_path in files:
        name = Path(file_path).name
        match = file_pattern.match(name)
        if not match:
            continue

        n_value = float(match.group("N"))
        grouped.setdefault(n_value, {"avalanches": [], "eliminated": [], "files": []})

        try:
            data = np.loadtxt(file_path, delimiter="\t", skiprows=1)
        except ValueError:
            continue

        if data.size == 0:
            continue

        if data.ndim == 1:
            data = data.reshape(1, -1)

        if data.shape[1] < 3:
            continue

        grouped[n_value]["avalanches"].append(data[:, 1])
        grouped[n_value]["eliminated"].append(data[:, 2])
        grouped[n_value]["files"].append(file_path)

    grouped = {
        n: {
            "avalanches": np.concatenate(values["avalanches"]),
            "eliminated": np.concatenate(values["eliminated"]),
            "files": values["files"],
        }
        for n, values in grouped.items()
        if values["avalanches"]
    }

    if not grouped:
        raise OSError("No usable rows in matched files")

    return grouped


def main():
    data_glob = "src/honeycombPuyo/outputs/avalanche2D/*.tsv"
    file_pattern = re.compile(r"^L_(?P<L>\d+)_N_(?P<N>[^_]+)_(?P<simNo>\d+)\.tsv$")
    grouped = load_grouped_by_n(data_glob, file_pattern)
    total_files = sum(len(values["files"]) for values in grouped.values())
    print(f"Loaded {total_files} files from {data_glob}")
    print(f"Found N values: {', '.join(str(n) for n in sorted(grouped.keys()))}")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    all_avalanches = np.concatenate([values["avalanches"] for values in grouped.values()])
    all_eliminated = np.concatenate([values["eliminated"] for values in grouped.values()])

    avalanche_bins = np.geomspace(1.0, max(2.0, float(np.max(all_avalanches))), 30)
    eliminated_bins = np.geomspace(1.0, max(2.0, float(np.max(all_eliminated))), 30)

    cmap = plt.get_cmap("rainbow", len(grouped))

    for i, n_value in enumerate(sorted(grouped.keys())):
        avalanches = grouped[n_value]["avalanches"]
        eliminated = grouped[n_value]["eliminated"]
        n_sims = len(grouped[n_value]["files"])
        color = cmap(i / max(1, len(grouped) - 1))

        avalanche_hist, _ = np.histogram(avalanches, bins=avalanche_bins)
        eliminated_hist, _ = np.histogram(eliminated, bins=eliminated_bins)

        avalanche_density = avalanche_hist / np.diff(avalanche_bins)
        eliminated_density = eliminated_hist / np.diff(eliminated_bins)

        ax1.plot(
            avalanche_bins[:-1][avalanche_density > 0],
            avalanche_density[avalanche_density > 0],
            marker="x",
            linewidth=1.2,
            label=f"N={n_value:g} ({n_sims} sims)",
            color=color,
        )

        ax2.plot(
            eliminated_bins[:-1][eliminated_density > 0],
            eliminated_density[eliminated_density > 0],
            marker="x",
            linewidth=1.2,
            label=f"N={n_value:g} ({n_sims} sims)",
            color=color,
        )

    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlabel("Avalanche size")
    ax1.set_ylabel("Frequency density")
    ax1.grid()
    ax1.legend(ncols=2)

    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlabel("Total eliminated")
    ax2.set_ylabel("Frequency density")
    ax2.grid()
    ax2.legend(ncols=2)

    out_path = Path("src/honeycombPuyo/plots/avalanche/simple_combined_distributions.png")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_path, dpi=300)


if __name__ == "__main__":
    main()