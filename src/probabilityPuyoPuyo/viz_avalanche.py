import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from tqdm import tqdm
import glob
import re

def plot_avalanche_distributions(L, S_list, avalanches_dict, total_elim_dict, plot_dir):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    cmap = plt.get_cmap("rainbow", len(S_list))
    for i, S in enumerate(S_list):
        avalanches = avalanches_dict.get(S)
        total_elim = total_elim_dict.get(S)
        if avalanches is None or total_elim is None:
            continue
        color = cmap(i)
        bins = np.geomspace(1, avalanches.max(), 20)
        hist, edges = np.histogram(avalanches, bins)
        hist = hist/np.diff(edges)
        ax1.plot(edges[:-1][hist>0], hist[hist>0],
                 label=f"S={S:.3f}", color=color, marker='x', alpha=0.8)
        bins2 = np.geomspace(1, total_elim.max(), 20)
        hist2, edges2 = np.histogram(total_elim, bins2)
        hist2 = hist2/np.diff(edges2)
        ax2.plot(edges2[:-1][hist2>0], hist2[hist2>0],
                 label=f"S={S:.3f}", color=color, marker='x', alpha=0.8)
    ax1.set_xscale("log"); ax1.set_yscale("log")
    ax1.set_xlabel("Avalanche Size"); ax1.set_ylabel("Frequency")
    ax1.grid(); ax1.legend(ncols=2)
    ax2.set_xscale("log"); ax2.set_yscale("log")
    ax2.set_xlabel("Total Eliminated"); ax2.set_ylabel("Frequency")
    ax2.grid(); ax2.legend(ncols=2)
    plt.tight_layout()
    plt.savefig(f"{plot_dir}/avalanche/distributions_L_{L}.png", dpi=300)
    plt.show()

def plot_time_evolution(L, S_list, t_dict, mass_dict, heights_dict, slopes_dict, plot_dir, n_bins=200):
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    ax1, ax2, ax3, ax4 = axes.flatten()
    cmap = plt.get_cmap("rainbow", len(S_list))

    for i, S in enumerate(S_list):
        t = t_dict.get(S)
        mass = mass_dict.get(S)
        heights = heights_dict.get(S)
        slopes = slopes_dict.get(S)
        if t is None or mass is None or heights is None or slopes is None:
            continue
        avg_h = heights.mean(axis=1)
        max_slope = np.max(np.abs(slopes), axis=1)
        roughness = np.mean(np.abs(heights - avg_h[:,None]), axis=1)
        color = cmap(i)

        # Top: plot as usual
        ax1.plot(t, mass, label=f"S={S:.3f}", color=color)
        ax2.plot(t, avg_h, label=f"S={S:.3f}", color=color)

        # Bottom: bin and plot with errorbars
        for arr, ax in zip([max_slope, roughness], [ax3, ax4]):
            # Bin edges
            bins = np.linspace(t[0], t[-1], n_bins+1)
            bin_centers = 0.5 * (bins[:-1] + bins[1:])
            means = []
            stderrs = []
            for j in range(n_bins):
                mask = (t >= bins[j]) & (t < bins[j+1])
                vals = arr[mask]
                if len(vals) > 0:
                    means.append(np.mean(vals))
                    stderrs.append(np.std(vals) / np.sqrt(len(vals)))
                else:
                    means.append(np.nan)
                    stderrs.append(np.nan)
            means = np.array(means)
            stderrs = np.array(stderrs)
            ax.plot(bin_centers, means, label=f"S={S:.3f}", color=color, alpha=0.8)

    ax1.set_xlabel("Time");        ax1.set_ylabel("Mass");           ax1.grid(); ax1.legend(ncols=2)
    ax2.set_xlabel("Time");        ax2.set_ylabel("Avg Height");     ax2.grid(); ax2.legend(ncols=2)
    ax3.set_xlabel("Time");        ax3.set_ylabel("Max |Slope|");    ax3.grid(); ax3.legend(ncols=2)
    ax4.set_xlabel("Time");        ax4.set_ylabel(r"$R_a$");         ax4.grid(); ax4.legend(ncols=2)
    plt.tight_layout()
    plt.savefig(f"{plot_dir}/avalanche/time_evolution_L_{L}.png", dpi=300)
    plt.show()

def timeAveragedSlopes(L, S_list, slopes_dict, plot_dir):
    fig, ax = plt.subplots(figsize=(10,6))
    cmap = plt.get_cmap("rainbow", len(S_list))
    for i, S in enumerate(S_list):
        slopes = slopes_dict.get(S)
        if slopes is None:
            continue
        all_slopes = slopes.flatten()
        ax.hist(all_slopes,
                bins=np.arange(-L, L+1),
                histtype="step", density=True,
                alpha=0.6, color=cmap(i),
                label=f"S={S:.3f}")
    ax.set_xlabel("Slope"); ax.set_ylabel("Density")
    ax.set_yscale("log"); ax.grid(); ax.legend(ncols=2)
    plt.tight_layout()
    plt.savefig(f"{plot_dir}/slopeDistribution/timeAverage_L_{L}.png", dpi=300)
    plt.show()

def plotInterface(L, S_list, t_dict, heights_dict, times, plot_dir):
    fig = plt.figure(figsize=(16,12))
    gs = GridSpec(len(S_list), len(times), figure=fig)
    for r, S in enumerate(S_list):
        t = t_dict.get(S)
        heights = heights_dict.get(S)
        if t is None or heights is None:
            continue
        for c, t_query in enumerate(times):
            idx = np.argmin(np.abs(t - t_query))
            actual_t = t[idx]
            h = heights[idx]
            ax = fig.add_subplot(gs[r, c])
            x = np.arange(L)
            ax.fill_between(x, np.min(h), h, step='mid', color='blue', alpha=0.3)
            if r == len(S_list)-1:
                ax.set_xlabel("Column")
            if c == 0:
                ax.set_ylabel(f"S={S:.3f}")
            if r == 0:
                ax.set_title(f"t≈{t_query} (actual {int(actual_t)})")
            ax.grid()
    plt.tight_layout()
    plt.savefig(f"{plot_dir}/slopeDistribution/interfaceGrid_L_{L}.png", dpi=300)
    plt.show()


def entropy_from_probs(probs):
    probs = np.array(probs)
    probs = probs[probs > 0]
    return -np.sum(probs * np.log(probs))

def main():
    L = 64
    data_dir = "src/probabilityPuyoPuyo/outputs/avalanche2D"
    plot_dir = "src/probabilityPuyoPuyo/plots"

    # Preload and preprocess all data
    t_dict = {}
    mass_dict = {}
    avalanches_dict = {}
    total_elim_dict = {}
    slopes_dict = {}
    heights_dict = {}
    S_list = []
    file_info = []

    print("Starting to load data...")

    files = glob.glob(f"{data_dir}/L_{L}_P_*.tsv")
    for file in tqdm(files):
        match = re.search(r'_P_([0-9\.\-]+)\.tsv$', file)
        if not match:
            print(f"Warning: could not parse probabilities from {file}")
            continue
        prob_str = match.group(1)
        probs = np.array([float(p) for p in prob_str.split('-')])
        probs /= probs.sum()
        S = entropy_from_probs(probs)
        S_list.append(S)
        file_info.append((S, file, probs))

    # Sort by entropy
    S_list = np.array(S_list)
    sort_idx = np.argsort(S_list)
    S_list = S_list[sort_idx]
    file_info = [file_info[i] for i in sort_idx]

    for S, file, probs in tqdm(file_info):
        try:
            raw = np.loadtxt(file, delimiter="\t", skiprows=1, dtype=str)
        except OSError:
            print(f"Warning: missing {file}")
            continue
        t = raw[:,0].astype(float)
        mass = raw[:,1].astype(int)
        avalanches = raw[:,2].astype(int)
        total_elim = raw[:,3].astype(int)
        first_h = raw[:,4].astype(int)
        slopes = np.array([list(map(int, s.split(","))) for s in raw[:,5]], dtype=int)
        heights = np.empty((len(first_h), L), dtype=float)
        heights[:,0] = first_h
        for c in range(L-1):
            heights[:,c+1] = heights[:,c] + slopes[:,c]
        t_dict[S] = t
        mass_dict[S] = mass
        avalanches_dict[S] = avalanches
        total_elim_dict[S] = total_elim
        slopes_dict[S] = slopes
        heights_dict[S] = heights

    print("Finished loading data.")
    print("Plotting avalanche distributions...")
    plot_avalanche_distributions(L, S_list, avalanches_dict, total_elim_dict, plot_dir)
    print("Plotting time evolution...")
    plot_time_evolution(L, S_list, t_dict, mass_dict, heights_dict, slopes_dict, plot_dir)
    print("Plotting time-averaged slopes...")
    timeAveragedSlopes(L, S_list, slopes_dict, plot_dir)
    print("Plotting interface evolution...")
    times = [256, 1024, 4096, 8192]
    plotInterface(L, S_list, t_dict, heights_dict, times, plot_dir)


def plot_avalanche_distributions_N(L, N_list, avalanches_dict, total_elim_dict, plot_dir):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    cmap = plt.get_cmap("rainbow", len(N_list))
    for i, N in enumerate(N_list):
        avalanches = avalanches_dict.get(N)
        total_elim = total_elim_dict.get(N)
        if avalanches is None or total_elim is None:
            continue
        color = cmap(i)
        bins = np.geomspace(1, avalanches.max(), 20)
        hist, edges = np.histogram(avalanches, bins)
        hist = hist/np.diff(edges)
        ax1.plot(edges[:-1][hist>0], hist[hist>0],
                 label=f"N={N:.3f}", color=color, marker='x', alpha=0.8)
        bins2 = np.geomspace(1, total_elim.max(), 20)
        hist2, edges2 = np.histogram(total_elim, bins2)
        hist2 = hist2/np.diff(edges2)
        ax2.plot(edges2[:-1][hist2>0], hist2[hist2>0],
                 label=f"N={N:.3f}", color=color, marker='x', alpha=0.8)
    ax1.set_xscale("log"); ax1.set_yscale("log")
    ax1.set_xlabel("Avalanche Size"); ax1.set_ylabel("Frequency")
    ax1.grid(); ax1.legend(ncols=2)
    ax2.set_xscale("log"); ax2.set_yscale("log")
    ax2.set_xlabel("Total Eliminated"); ax2.set_ylabel("Frequency")
    ax2.grid(); ax2.legend(ncols=2)
    plt.tight_layout()
    plt.savefig(f"{plot_dir}/avalanche/distributions_N_L_{L}.png", dpi=300)
    plt.show()

def mainOnlyAvalanche():
    L = 64
    data_dir = "src/probabilityPuyoPuyo/outputs/avalanche2D/onlyAvalanche"
    plot_dir = "src/probabilityPuyoPuyo/plots"

    avalanches_dict = {}
    total_elim_dict = {}
    N_list = []
    file_info = []

    print("Starting to load data...")

    # Updated glob and regex for files like L_64_N_5.1.tsv
    files = glob.glob(f"{data_dir}/L_{L}_N_*.tsv")
    for file in tqdm(files):
        match = re.search(r'_N_([0-9\.]+)\.tsv$', file)
        if not match:
            print(f"Warning: could not parse N from {file}")
            continue
        N = float(match.group(1))
        N_list.append(N)
        file_info.append((N, file))

    # Sort by N
    N_list = np.array(N_list)
    sort_idx = np.argsort(N_list)
    N_list = N_list[sort_idx]
    file_info = [file_info[i] for i in sort_idx]

    for N, file in tqdm(file_info):
        try:
            raw = np.loadtxt(file, delimiter="\t", skiprows=1, dtype=str)
        except OSError:
            print(f"Warning: missing {file}")
            continue
        avalanches = raw[:,1].astype(int)
        total_elim = raw[:,2].astype(int)
        avalanches_dict[N] = avalanches
        total_elim_dict[N] = total_elim

    plot_avalanche_distributions_N(L, N_list, avalanches_dict, total_elim_dict, plot_dir)

def mainOnlyAvalancheRandomProbs():
    L = 64
    data_dir = "src/probabilityPuyoPuyo/outputs/avalanche2D/onlyAvalanche"
    plot_dir = "src/probabilityPuyoPuyo/plots"

    # Preload and preprocess all data
    avalanches_dict = {}
    total_elim_dict = {}
    S_list = []
    file_info = []

    print("Starting to load data...")

    files = glob.glob(f"{data_dir}/L_{L}_P_*.tsv")
    for file in tqdm(files):
        match = re.search(r'_P_([0-9\.\-]+)\.tsv$', file)
        if not match:
            print(f"Warning: could not parse probabilities from {file}")
            continue
        prob_str = match.group(1)
        probs = np.array([float(p) for p in prob_str.split('-')])
        probs /= probs.sum()
        S = entropy_from_probs(probs)
        S_list.append(S)
        file_info.append((S, file, probs))

    # Sort by entropy
    S_list = np.array(S_list)
    sort_idx = np.argsort(S_list)
    S_list = S_list[sort_idx]
    file_info = [file_info[i] for i in sort_idx]

    for S, file, probs in tqdm(file_info):
        try:
            raw = np.loadtxt(file, delimiter="\t", skiprows=1, dtype=str)
        except OSError:
            print(f"Warning: missing {file}")
            continue
        avalanches = raw[:,1].astype(int)
        total_elim = raw[:,2].astype(int)
        avalanches_dict[S] = avalanches
        total_elim_dict[S] = total_elim

    plot_avalanche_distributions(L, S_list, avalanches_dict, total_elim_dict, plot_dir)

if __name__ == "__main__":
    # main()
    mainOnlyAvalanche()
    # mainOnlyAvalancheRandomProbs()