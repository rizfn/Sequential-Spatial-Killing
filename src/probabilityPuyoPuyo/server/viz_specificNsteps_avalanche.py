import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import glob
import re


def plot_avalanche_distributions_N(L, N_steps, N_list, avalanches_dict, total_elim_dict, plot_dir):
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
    plt.savefig(f"{plot_dir}/distributions_L_{L}_steps_{N_steps}_2D.png", dpi=300)
    plt.show()

def mainOnlyAvalanche(L, N_steps):
    data_dir = "/nbi/nbicmplx/cell/rpw391/probabilityPuyo/avalanche/outputs/2D"
    plot_dir = "plots"

    avalanches_dict = {}
    total_elim_dict = {}
    N_list = []
    file_info = []

    print("Starting to load data...")

    # Updated glob and regex for files like L_64_N_5.1_steps_32768.tsv
    files = glob.glob(f"{data_dir}/L_{L}_N_*_steps_{N_steps}.tsv")
    for file in tqdm(files):
        match = re.search(r'_N_([0-9\.]+)_steps_(\d+)\.tsv$', file)
        if not match:
            print(f"Warning: could not parse N from {file}")
            continue
        N = float(match.group(1))
        steps = int(match.group(2))
        if steps != N_steps:
            continue
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

    plot_avalanche_distributions_N(L, N_steps, N_list, avalanches_dict, total_elim_dict, plot_dir)

if __name__ == '__main__':
    mainOnlyAvalanche(64, 32768*8)
