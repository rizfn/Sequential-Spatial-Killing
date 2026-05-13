import numpy as np
import matplotlib.pyplot as plt
from glob import glob
from pathlib import Path
import re

def load_grouped_by_rho(data_glob, file_pattern, transient_time=0.1):
    files = sorted(glob(data_glob))
    if not files:
        raise OSError(f"No files found for pattern: {data_glob}")

    grouped = {}
    l_values = set()
    n_values = set()
    final_times = []
    for file_path in files:
        name = Path(file_path).name
        match = file_pattern.match(name)
        if not match:
            continue

        l_values.add(int(match.group("L")))
        n_values.add(int(match.group("N")))
        rho_value = float(match.group("rho"))
        grouped.setdefault(rho_value, {"avalanches": [], "eliminated": [], "files": []})

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

        # Skip transient based on time (first column)
        time_col = data[:, 0]
        max_time = time_col[-1]
        final_times.append(max_time)
        cutoff_time = transient_time * max_time
        mask = time_col >= cutoff_time
        data_steady = data[mask, :]

        grouped[rho_value]["avalanches"].append(data_steady[:, 1])
        grouped[rho_value]["eliminated"].append(data_steady[:, 2])
        grouped[rho_value]["files"].append(file_path)

    grouped = {
        rho: {
            "avalanches": np.concatenate(values["avalanches"]),
            "eliminated": np.concatenate(values["eliminated"]),
            "files": values["files"],
        }
        for rho, values in grouped.items()
        if values["avalanches"]
    }

    if not grouped:
        raise OSError("No usable rows in matched files")

    if len(l_values) != 1 or len(n_values) != 1:
        raise ValueError(f"Expected a single L and N value, found L={sorted(l_values)} N={sorted(n_values)}")

    return grouped, next(iter(l_values)), next(iter(n_values)), int(np.ceil(max(final_times)))


def main():
    data_glob = "src/langmuirRandom/outputs/avalanche2D/*.tsv"
    file_pattern = re.compile(r"^L_(?P<L>\d+)_N_(?P<N>\d+)_rho_(?P<rho>[^_]+)_(?P<simNo>\d+)\.tsv$")
    grouped, L, N, steps = load_grouped_by_rho(data_glob, file_pattern, transient_time=0.1)
    
    total_files = sum(len(values["files"]) for values in grouped.values())
    print(f"Loaded {total_files} files from {data_glob}")
    print(f"Found rho values: {', '.join(f'{rho:.4f}' for rho in sorted(grouped.keys()))}")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    all_avalanches = np.concatenate([values["avalanches"] for values in grouped.values()])
    all_eliminated = np.concatenate([values["eliminated"] for values in grouped.values()])

    avalanche_bins = np.geomspace(1.0, max(2.0, float(np.max(all_avalanches))), 30)
    eliminated_bins = np.geomspace(1.0, max(2.0, float(np.max(all_eliminated))), 30)

    cmap = plt.get_cmap("rainbow", len(grouped))

    for i, rho_value in enumerate(sorted(grouped.keys())):
        avalanches = grouped[rho_value]["avalanches"]
        eliminated = grouped[rho_value]["eliminated"]
        n_sims = len(grouped[rho_value]["files"])
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
            label=f"rho={rho_value:.4f} ({n_sims} sims)",
            color=color,
        )

        ax2.plot(
            eliminated_bins[:-1][eliminated_density > 0],
            eliminated_density[eliminated_density > 0],
            marker="x",
            linewidth=1.2,
            label=f"rho={rho_value:.4f} ({n_sims} sims)",
            color=color,
        )

    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlabel("Avalanche size", fontsize=11)
    ax1.set_ylabel("Frequency density", fontsize=11)
    ax1.set_title("Avalanche size distribution")
    ax1.grid(alpha=0.3)
    ax1.legend(fontsize=9)

    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlabel("Total eliminated", fontsize=11)
    ax2.set_ylabel("Frequency density", fontsize=11)
    ax2.set_title("Total eliminated distribution")
    ax2.grid(alpha=0.3)
    ax2.legend(fontsize=9)

    out_path = Path(f"src/langmuirRandom/plots/avalanche/distributions_L_{L}_N_{N}_steps_{steps}.png")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()