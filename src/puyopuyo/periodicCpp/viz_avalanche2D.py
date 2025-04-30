import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from tqdm import tqdm

def plot_avalanche_distributions(L, N_list, avalanches_dict, total_elim_dict, plot_dir):
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
                 label=f"N={N}", color=color, marker='x', alpha=0.8)
        bins2 = np.geomspace(1, total_elim.max(), 20)
        hist2, edges2 = np.histogram(total_elim, bins2)
        hist2 = hist2/np.diff(edges2)
        ax2.plot(edges2[:-1][hist2>0], hist2[hist2>0],
                 label=f"N={N}", color=color, marker='x', alpha=0.8)
    ax1.set_xscale("log"); ax1.set_yscale("log")
    ax1.set_xlabel("Avalanche Size"); ax1.set_ylabel("Frequency")
    ax1.grid(); ax1.legend(ncols=2)
    ax2.set_xscale("log"); ax2.set_yscale("log")
    ax2.set_xlabel("Total Eliminated"); ax2.set_ylabel("Frequency")
    ax2.grid(); ax2.legend(ncols=2)
    plt.tight_layout()
    plt.savefig(f"{plot_dir}/avalanche/distributions_L_{L}.png", dpi=300)
    plt.show()

def plot_time_evolution(L, N_list, t_dict, mass_dict, heights_dict, slopes_dict, plot_dir):
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    ax1, ax2, ax3, ax4 = axes.flatten()
    cmap = plt.get_cmap("rainbow", len(N_list))
    for i, N in enumerate(N_list):
        t = t_dict.get(N)
        mass = mass_dict.get(N)
        heights = heights_dict.get(N)
        slopes = slopes_dict.get(N)
        if t is None or mass is None or heights is None or slopes is None:
            continue
        avg_h = heights.mean(axis=1)
        max_slope = np.max(np.abs(slopes), axis=1)
        roughness = np.mean(np.abs(heights - avg_h[:,None]), axis=1)
        color = cmap(i)
        ax1.plot(t, mass, label=f"N={N}", color=color)
        ax2.plot(t, avg_h, label=f"N={N}", color=color)
        ax3.plot(t, max_slope, label=f"N={N}", color=color)
        ax4.plot(t, roughness, label=f"N={N}", color=color)
        # ax3.set_xscale("log");    ax3.set_yscale("log")
        # ax4.set_xscale("log");    ax4.set_yscale("log")

    ax1.set_xlabel("Time");        ax1.set_ylabel("Mass");           ax1.grid(); ax1.legend(ncols=2)
    ax2.set_xlabel("Time");        ax2.set_ylabel("Avg Height");     ax2.grid(); ax2.legend(ncols=2)
    ax3.set_xlabel("Time");        ax3.set_ylabel("Max |Slope|");    ax3.grid(); ax3.legend(ncols=2)
    ax4.set_xlabel("Time");        ax4.set_ylabel(r"$R_a$");         ax4.grid(); ax4.legend(ncols=2)
    plt.tight_layout()
    plt.savefig(f"{plot_dir}/avalanche/time_evolution_L_{L}.png", dpi=300)
    plt.show()

def timeAveragedSlopes(L, N_list, slopes_dict, plot_dir):
    fig, ax = plt.subplots(figsize=(10,6))
    cmap = plt.get_cmap("rainbow", len(N_list))
    for i, N in enumerate(N_list):
        slopes = slopes_dict.get(N)
        if slopes is None:
            continue
        all_slopes = slopes.flatten()
        ax.hist(all_slopes,
                bins=np.arange(-L, L+1),
                histtype="step", density=True,
                alpha=0.6, color=cmap(i),
                label=f"N={N}")
    ax.set_xlabel("Slope"); ax.set_ylabel("Density")
    ax.set_yscale("log"); ax.grid(); ax.legend(ncols=2)
    plt.tight_layout()
    plt.savefig(f"{plot_dir}/slopeDistribution/timeAverage_L_{L}.png", dpi=300)
    plt.show()

def plotInterface(L, N_list, t_dict, heights_dict, times, plot_dir):
    fig = plt.figure(figsize=(16,12))
    gs = GridSpec(len(N_list), len(times), figure=fig)
    for r, N in enumerate(N_list):
        t = t_dict.get(N)
        heights = heights_dict.get(N)
        if t is None or heights is None:
            continue
        for c, t_query in enumerate(times):
            idx = np.argmin(np.abs(t - t_query))
            actual_t = t[idx]
            h = heights[idx]
            ax = fig.add_subplot(gs[r, c])
            x = np.arange(L)
            ax.fill_between(x, np.min(h), h, step='mid', color='blue', alpha=0.3)
            if r == len(N_list)-1:
                ax.set_xlabel("Column")
            if c == 0:
                ax.set_ylabel(f"N={N}")
            if r == 0:
                ax.set_title(f"t≈{t_query} (actual {int(actual_t)})")
            ax.grid()
    plt.tight_layout()
    plt.savefig(f"{plot_dir}/slopeDistribution/interfaceGrid_L_{L}.png", dpi=300)
    plt.show()

def main():
    L = 64
    N_list = np.arange(2, 13)
    steps = 1024 * 8
    data_dir = "src/puyopuyo/periodicCpp/outputs"
    plot_dir = "src/puyopuyo/periodicCpp/plots"

    # Preload and preprocess all data
    t_dict = {}
    mass_dict = {}
    avalanches_dict = {}
    total_elim_dict = {}
    slopes_dict = {}
    heights_dict = {}

    print("Starting to load data...")

    for N in tqdm(N_list):
        fn = f"{data_dir}/avalanche2D/L_{L}_N_{N}_steps_{steps}.tsv"
        try:
            raw = np.loadtxt(fn, delimiter="\t", skiprows=1, dtype=str)
        except OSError:
            print(f"Warning: missing {fn}")
            continue
        t = raw[:,0].astype(float)
        mass = raw[:,1].astype(int)
        avalanches = raw[:,2].astype(int)
        total_elim = raw[:,3].astype(int)
        first_h = raw[:,4].astype(int)
        # slopes: shape (num_steps, L-1)
        slopes = np.array([list(map(int, s.split(","))) for s in raw[:,5]], dtype=int)
        # reconstruct heights: shape (num_steps, L)
        heights = np.empty((len(first_h), L), dtype=float)
        heights[:,0] = first_h
        for c in range(L-1):
            heights[:,c+1] = heights[:,c] + slopes[:,c]
        t_dict[N] = t
        mass_dict[N] = mass
        avalanches_dict[N] = avalanches
        total_elim_dict[N] = total_elim
        slopes_dict[N] = slopes
        heights_dict[N] = heights

    print("Finished loading data.")
    print("Plotting avalanche distributions...")
    plot_avalanche_distributions(L, N_list, avalanches_dict, total_elim_dict, plot_dir)
    print("Plotting time evolution...")
    plot_time_evolution(L, N_list, t_dict, mass_dict, heights_dict, slopes_dict, plot_dir)
    print("Plotting time-averaged slopes...")
    timeAveragedSlopes(L, N_list, slopes_dict, plot_dir)
    print("Plotting interface evolution...")
    times = [256, 1024, 4096, 8192]
    plotInterface(L, N_list, t_dict, heights_dict, times, plot_dir)

if __name__ == "__main__":
    main()
