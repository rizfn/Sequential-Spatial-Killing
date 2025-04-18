import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm


def update(matrix, currentVictim):
    # Iterate from the left boundary
    for i in range(len(matrix)):
        if matrix[i] == 0:
            continue
        elif matrix[i] == currentVictim:
            matrix[i] = 0
        else:
            break

    # Iterate from the right boundary
    for i in range(len(matrix) - 1, -1, -1):
        if matrix[i] == 0:
            continue
        elif matrix[i] == currentVictim:
            matrix[i] = 0
        else:
            break

    return matrix


def main():
    N_species_values = [2, 4, 6, 8, 10, 12, 40]  # Multiple values of N_species
    L_array = np.arange(100, 410, 10)
    max_steps = 10000
    N_simulations = 10  # Number of simulations for each (L, N_species)

    for N_species in tqdm(N_species_values):
        t_means = []
        t_stds = []

        for L in L_array:
            t_values = []

            for _ in range(N_simulations):
                # Generate a 1D array with random values
                matrix = np.random.randint(1, N_species + 1, (L,))
                # Pad the array with zeros on both sides
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

    plt.title(f"1D, {N_simulations} sims per datapoint")
    plt.xlabel("System size L")
    plt.ylabel("Number of sequential waves until eliminated (normalized)")
    plt.legend()
    plt.grid()
    plt.savefig("src/percolation/plots/systemSize/1D.png", dpi=300)
    plt.show()


if __name__ == "__main__":
    main()