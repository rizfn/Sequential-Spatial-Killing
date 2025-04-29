import numpy as np
import matplotlib.pyplot as plt

def viz_2D():
    L = 64
    N_list = np.arange(2, 12)
    steps = 1024

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 12))


    # Use a colormap (e.g., rainbow) to assign colors based on N_species
    colormap = plt.get_cmap("rainbow", len(N_list))

    for i, N in enumerate(N_list):
        try:
            # Load data
            step, mass, height, avalanches, total_eliminated = np.loadtxt(
                f"src/puyopuyo/cpp/outputs/avalanche2D/L_{L}_N_{N}_steps_{steps}.tsv",
                delimiter="\t",
                skiprows=1,
                unpack=True,
            )
            color = colormap(i / len(N_list))  # Normalize the index to [0, 1] for the colormap

            # Plot mass and height
            ax1.plot(step, mass, label=f"N={N}", color=color)
            ax2.plot(step, height, label=f"N={N}", color=color)

            # Create histograms for avalanche sizes and total eliminated puyos
            avalanche_bins = np.geomspace(1, np.max(avalanches), 20)
            avalanche_hist, avalanche_bins = np.histogram(avalanches, avalanche_bins)
            avalanche_hist = avalanche_hist / np.diff(avalanche_bins)
            eliminated_bins = np.geomspace(1, np.max(total_eliminated), 20)
            eliminated_hist, eliminated_bins = np.histogram(total_eliminated, eliminated_bins)
            eliminated_hist = eliminated_hist / np.diff(eliminated_bins)

            # Plot histograms
            ax3.plot(avalanche_bins[:-1][avalanche_hist!=0], avalanche_hist[avalanche_hist!=0], label=f"N={N}", color=color, marker='x', alpha=0.9)
            ax4.plot(eliminated_bins[:-1][eliminated_hist!=0], eliminated_hist[eliminated_hist!=0], label=f"N={N}", color=color, marker='x', alpha=0.9)

        except OSError:
            print(f"Warning: File for L={L}, N={N} not found. Skipping.")

    ylim1, ylim2 = ax1.get_ylim(), ax2.get_ylim()

    ax1.set_ylim(ylim1)
    ax2.set_ylim(ylim2)

    ax1.plot([0, steps], [L, L * (steps + 1)], color="grey", linestyle="--", label="Total mass added")
    ax2.plot([0, steps], [1, steps + 1], color="grey", linestyle="--", label="Average (dense) height")

    tau = 3
    x_vals = np.linspace(1, ax3.get_xlim()[1], 100)
    ax3.plot(x_vals, 1e5*x_vals**(-tau), label=f"$\\tau$={tau}", color="grey", linestyle="--", alpha=0.7)
    ax4.plot(x_vals, 1e5*x_vals**(-tau), label=f"$\\tau$={tau}", color="grey", linestyle="--", alpha=0.7)

    # Configure ax1 (Mass vs Time)
    ax1.set_xlabel("Time")
    ax1.set_ylabel("Mass")
    ax1.grid()
    ax1.legend(ncols=2)

    # Configure ax2 (Max Height vs Time)
    ax2.set_xlabel("Time")
    ax2.set_ylabel("Max Height")
    ax2.grid()
    ax2.legend(ncols=2)

    # Configure ax3 (Histogram of Avalanche Sizes)
    ax3.set_xlabel("Avalanche Size")
    ax3.set_ylabel("Frequency")
    ax3.set_xscale("log")
    ax3.set_yscale("log")
    ax3.grid()
    ax3.legend(ncols=2)

    # Configure ax4 (Histogram of Total Eliminated Puyos)
    ax4.set_xlabel("Number of Eliminated Puyos")
    ax4.set_ylabel("Frequency")
    ax4.set_xscale("log")
    ax4.set_yscale("log")
    ax4.grid()
    ax4.legend(ncols=2)

    # Save the plot
    plt.tight_layout()
    plt.savefig(f"src/puyopuyo/cpp/plots/avalanche/2D_L_{L}.png", dpi=300)
    plt.show()


def avg_and_max_avalanche_sizes_vs_N():
    L = 64
    N_list = np.arange(2, 12)
    steps = 1024

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # Use a colormap (e.g., rainbow) to assign colors based on N_species
    colormap = plt.get_cmap("rainbow", len(N_list))

    avg_sizes = []
    max_sizes = []

    for i, N in enumerate(N_list):
        try:
            # Load data
            step, mass, height, avalanches, total_eliminated = np.loadtxt(
                f"src/puyopuyo/cpp/outputs/avalanche2D/L_{L}_N_{N}_steps_{steps}.tsv",
                delimiter="\t",
                skiprows=1,
                unpack=True,
            )

            # Calculate average and maximum avalanche sizes
            avg_size = np.mean(avalanches)
            max_size = np.max(avalanches)

            avg_sizes.append(avg_size)
            max_sizes.append(max_size)

        except OSError:
            print(f"Warning: File for L={L}, N={N} not found. Skipping.")
            avg_sizes.append(None)
            max_sizes.append(None)

    # Filter out None values for plotting
    valid_N = [N for i, N in enumerate(N_list) if avg_sizes[i] is not None]
    valid_avg_sizes = [size for size in avg_sizes if size is not None]
    valid_max_sizes = [size for size in max_sizes if size is not None]

    # Plot average avalanche size vs. N
    ax1.plot(valid_N, valid_avg_sizes, marker="o", color="blue", label="Average Avalanche Size")
    ax1.set_xlabel("Number of Species (N)")
    ax1.set_ylabel("Average Avalanche Size")
    ax1.set_title("Average Avalanche Size vs. N")
    ax1.grid()
    ax1.legend()

    # Plot largest avalanche size vs. N
    ax2.plot(valid_N, valid_max_sizes, marker="o", color="red", label="Largest Avalanche Size")
    ax2.set_xlabel("Number of Species (N)")
    ax2.set_ylabel("Largest Avalanche Size")
    ax2.set_title("Largest Avalanche Size vs. N")
    ax2.grid()
    ax2.legend()

    # Save the plot
    plt.tight_layout()
    plt.savefig(f"src/puyopuyo/cpp/plots/avalanche/avg_and_max_sizes_vs_N_{L}.png", dpi=300)
    plt.show()


if __name__ == "__main__":
    avg_and_max_avalanche_sizes_vs_N()
    # viz_2D()