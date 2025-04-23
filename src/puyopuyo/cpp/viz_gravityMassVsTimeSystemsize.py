import numpy as np
import matplotlib.pyplot as plt

import numpy as np
import matplotlib.pyplot as plt


def viz_2D():
    L_values = [1, 2, 3, 4, 5, 8, 16, 32, 64]
    N_list = np.arange(2, 8)
    steps = 1024

    for L in L_values:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

        # Use a colormap (e.g., rainbow) to assign colors based on N_species
        colormap = plt.get_cmap("rainbow", len(N_list))

        for i, N in enumerate(N_list):
            step, mass, height = np.loadtxt(
                f"src/puyopuyo/cpp/outputs/gravity2D_L/L_{L}_N_{N}_steps_{steps}.tsv",
                delimiter="\t",
                skiprows=1,
                unpack=True,
            )
            color = colormap(i / len(N_list))  # Normalize the index to [0, 1] for the colormap
            ax1.plot(step, mass, label=f"N={N}", color=color)
            ax2.plot(step, height, label=f"N={N}", color=color)
        ylim1 = ax1.get_ylim()
        ylim2 = ax2.get_ylim()
        ax1.plot([0, steps], [L, L * (steps + 1)], color="grey", linestyle="--", label="Total mass added")
        ax2.plot([0, steps], [1, steps + 1], color="grey", linestyle="--", label="Average (dense) height")
        ax1.set_ylim(ylim1)
        ax2.set_ylim(ylim2)


        ax1.set_xlabel("Time")
        ax1.set_ylabel("Mass")
        ax1.grid()
        ax1.legend(ncols=2)
        ax2.set_xlabel("Time")
        ax2.set_ylabel("Max Height")
        ax2.grid()
        ax2.legend(ncols=2)

        # Save the plot for the current L
        plt.savefig(f"src/puyopuyo/cpp/plots/2DSystemSize/MassVsTime_L_{L}.png", dpi=300)



if __name__ == "__main__":
    viz_2D()
