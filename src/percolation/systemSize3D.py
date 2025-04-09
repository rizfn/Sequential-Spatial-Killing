import numpy as np
import matplotlib.pyplot as plt
import cc3d
from tqdm import tqdm

def update(matrix, currentVictim):
    # Make a copy of the matrix, with 0 set to currentVictim
    copiedMatrix = matrix.copy()
    copiedMatrix[copiedMatrix == 0] = currentVictim
    # Find the connected components of the copied matrix
    labels = cc3d.connected_components(copiedMatrix, connectivity=6)  # 6-connectivity for 3D
    # Look for the cluster at the boundaries of 0. That cluster should be copied onto `matrix` (as zeros). Other clusters of 0 (due to currentVictim being set to 0) should be ignored
    boundaryCluster = labels[0, 0, 0]
    matrix[labels == boundaryCluster] = 0
    return matrix


def main():
    N_species_values = [2, 4, 6, 8, 10, 12, 40]  # Multiple values of N_species
    L_array = np.arange(10, 110, 10)
    max_steps = 10000
    N_simulations = 10  # Number of simulations for each (L, N_species)

    for N_species in tqdm(N_species_values):
        t_means = []
        t_stds = []

        for L in L_array:
            t_values = []

            for _ in range(N_simulations):
                # Generate a 3D matrix with random values
                matrix = np.random.randint(1, N_species + 1, (L, L, L))
                # Pad the matrix with zeros on all sides
                matrix = np.pad(matrix, pad_width=1, mode='constant', constant_values=0)
                victimOrder = np.arange(1, N_species + 1)

                for t in range(max_steps):
                    currentVictim = victimOrder[t % N_species]
                    matrix = update(matrix, currentVictim)
                    if np.count_nonzero(matrix) == 0:
                        t_values.append(t)
                        break

            # Compute mean and standard deviation of t for this (L, N_species)
            t_values = np.array(t_values) / N_species

            t_means.append(np.mean(t_values))
            t_stds.append(np.std(t_values))

        # Plot with error bars
        plt.errorbar(
            L_array,
            np.array(t_means),
            yerr=np.array(t_stds),
            capsize=5,
            label=f"{N_species} species"
        )

    plt.title(f"3D, {N_simulations} sims per datapoint")
    plt.xlabel("System size L")
    plt.ylabel("Number of sequential waves until eliminated (normalized)")
    plt.legend()
    plt.grid()
    plt.savefig("src/percolation/plots/systemSize/3D.png", dpi=300)
    plt.show()


if __name__ == "__main__":
    main()