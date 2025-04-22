import numpy as np
import matplotlib.pyplot as plt

def viz_1D():
    N_list = np.arange(2, 21)
    steps = 1024

    fig, ax1 = plt.subplots(1, 1, figsize=(7, 6))

    ax1.plot([0, steps], [1, (steps+1)], color="grey", linestyle="--", label="Total mass added")

    # Use a colormap (e.g., rainbow) to assign colors based on N_species
    colormap = plt.get_cmap("rainbow", len(N_list))

    for i, N in enumerate(N_list):
        step, mass = np.loadtxt(f"src/puyopuyo/cpp/outputs/gravity1D/N_{N}_steps_{steps}.tsv", delimiter="\t", skiprows=1, unpack=True)
        color = colormap(i / len(N_list))  # Normalize the index to [0, 1] for the colormap
        ax1.plot(step, mass, label=f"N={N}", color=color)

    ax1.set_xlabel("Time")
    ax1.set_ylabel("Mass")
    ax1.grid()
    ax1.legend(ncols=2)

    plt.savefig(f"src/puyopuyo/cpp/plots/1D_gravityMassVsTime.png", dpi=300)
    plt.show()


def viz_2D():
    L = 128
    N_list = np.arange(2, 21)
    steps = 1024

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    ax1.plot([0, steps], [L, L * (steps+1)], color="grey", linestyle="--", label="Total mass added")
    ax2.plot([0, steps], [1, steps+1], color="grey", linestyle="--", label="Average (dense) height")

    # Use a colormap (e.g., rainbow) to assign colors based on N_species
    colormap = plt.get_cmap("rainbow", len(N_list))

    for i, N in enumerate(N_list):
        step, mass, height = np.loadtxt(f"src/puyopuyo/cpp/outputs/gravity2D/L_{L}_N_{N}_steps_{steps}.tsv", delimiter="\t", skiprows=1, unpack=True)
        color = colormap(i / len(N_list))  # Normalize the index to [0, 1] for the colormap
        ax1.plot(step, mass, label=f"N={N}", color=color)
        ax2.plot(step, height, label=f"N={N}", color=color)

    ax1.set_xlabel("Time")
    ax1.set_ylabel("Mass")
    ax1.grid()
    ax1.legend(ncols=2)
    ax2.set_xlabel("Time")
    ax2.set_ylabel("Max Height")
    ax2.grid()
    ax2.legend(ncols=2)

    plt.savefig(f"src/puyopuyo/cpp/plots/2D_gravityMassVsTime_L_{L}.png", dpi=300)
    plt.show()


def viz_3D():
    L = 24
    N_list = np.arange(2, 21)
    steps = 1024

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))


    # Use a colormap (e.g., rainbow) to assign colors based on N_species
    colormap = plt.get_cmap("rainbow", len(N_list))

    for i, N in enumerate(N_list):
        step, mass, height = np.loadtxt(f"src/puyopuyo/cpp/outputs/gravity3D/L_{L}_N_{N}_steps_{steps}.tsv", delimiter="\t", skiprows=1, unpack=True)
        color = colormap(i / len(N_list))  # Normalize the index to [0, 1] for the colormap
        ax1.plot(step, mass, label=f"N={N}", color=color)
        ax2.plot(step, height, label=f"N={N}", color=color)

    ax1_ylim = ax1.get_ylim()
    ax1.plot([0, steps], [L*L, L*L * (steps+1)], color="grey", linestyle="--", label="Total mass added")
    ax1.set_ylim(ax1_ylim)
    ax2_ylim = ax2.get_ylim()
    ax2.plot([0, steps], [1, steps+1], color="grey", linestyle="--", label="Average (dense) height")
    ax2.set_ylim(ax2_ylim)

    ax1.set_xlabel("Time")
    ax1.set_ylabel("Mass")
    ax1.grid()
    ax1.legend(ncols=2)
    ax2.set_xlabel("Time")
    ax2.set_ylabel("Max Height")
    ax2.grid()
    ax2.legend(ncols=2)

    plt.savefig(f"src/puyopuyo/cpp/plots/3D_gravityMassVsTime_L_{L}.png", dpi=300)
    plt.show()


if __name__ == "__main__":
    # viz_1D()
    # viz_2D()
    viz_3D()