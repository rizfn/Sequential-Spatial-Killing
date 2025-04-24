import numpy as np
import matplotlib.pyplot as plt

def main():
    L = 128
    beta = 1/3
    N_list = np.arange(2, 13)
    steps = 1024

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(21, 6))

    ax1.plot([0, steps], [L, L * (steps+1)], color="grey", linestyle="--", label="Total mass added")
    ax2.plot([0, steps], [1, steps+1], color="grey", linestyle="--", label="Average (dense) height")
    ax3.plot([1, steps], 2*np.array([1, steps])**(beta), label=f"$\\beta$={beta:.2f}", color="grey", linestyle="--", alpha=0.7)

    # Use a colormap (e.g., rainbow) to assign colors based on N_species
    colormap = plt.get_cmap("rainbow", len(N_list))

    for i, N in enumerate(N_list):
        step, mass, height, roughness = np.loadtxt(f"src/puyopuyo/cpp/outputs/roughness2D/L_{L}_N_{N}_steps_{steps}.tsv", delimiter="\t", skiprows=1, unpack=True)
        color = colormap(i / len(N_list))  # Normalize the index to [0, 1] for the colormap
        ax1.plot(step, mass, label=f"N={N}", color=color)
        ax2.plot(step, height, label=f"N={N}", color=color)
        ax3.plot(step, roughness, label=f"N={N}", color=color)

    # plot N_inf
    step, mass, height, roughness = np.loadtxt(f"src/puyopuyo/cpp/outputs/roughness2D/L_{L}_N_inf_steps_{steps}.tsv", delimiter="\t", skiprows=1, unpack=True)
    ax1.plot(step, mass, label=f"N=inf", color="black", linestyle="--")
    ax2.plot(step, height, label=f"N=inf", color="black", linestyle="--")
    ax3.plot(step, roughness, label=f"N=inf", color="black", linestyle="--")


    ax1.set_xlabel("Time")
    ax1.set_ylabel("Mass")
    ax1.grid()
    ax1.legend(ncols=2)
    ax2.set_xlabel("Time")
    ax2.set_ylabel("Max Height")
    ax2.grid()
    ax2.legend(ncols=2)
    ax3.set_xlabel("Time")
    ax3.set_ylabel("Roughness")
    ax3.set_xscale("log")
    ax3.set_yscale("log")
    ax3.grid()
    ax3.legend(ncols=2)


    plt.savefig(f"src/puyopuyo/cpp/plots/roughness/2D_L_{L}.png", dpi=300)
    plt.show()


if __name__ == "__main__":
    main()