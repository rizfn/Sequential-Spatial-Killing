import numpy as np
import matplotlib.pyplot as plt
import cc3d
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


def generate_zipf_population(N_species, L, tau):
    """
    Generate a Zipf-distributed population for N_species.
    The total population is L.
    """
    # Generate Zipf-distributed weights
    weights = np.arange(1, N_species + 1) ** (-1/tau - 1)  # Zipf distribution formula
    weights /= np.sum(weights)  # Normalize to sum to 1
    populations = (weights * L).astype(int)  # Scale to total population size

    # Ensure the total population matches exactly
    diff = L - np.sum(populations)
    for i in range(abs(diff)):
        populations[i % N_species] += np.sign(diff)

    return populations


def main():
    N_species_values = [2, 4, 6, 8, 10, 12, 40]  # Multiple values of N_species
    L_array = np.arange(10, 410, 10)
    max_steps = 10000
    N_simulations = 10  # Number of simulations for each (L, N_species)
    tau = 2.0  # species abundances distribution exponent (plot of N_species vs Abundances)

    for N_species in tqdm(N_species_values):
        t_means = []
        t_stds = []

        for L in L_array:
            t_values = []

            for _ in range(N_simulations):
                # Generate Zipf-distributed populations
                populations = generate_zipf_population(N_species, L, tau)

                # Create the array with the specified populations
                array = np.zeros(L, dtype=int)
                species = np.arange(1, N_species + 1)
                indices = np.arange(L)
                np.random.shuffle(indices)
                start = 0
                for i, pop in enumerate(populations):
                    array[indices[start:start + pop]] = species[i]
                    start += pop

                victimOrder = np.arange(1, N_species + 1)

                for t in range(max_steps):
                    currentVictim = victimOrder[t % N_species]
                    array = update(array, currentVictim)
                    if np.count_nonzero(array) == 0:
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

    plt.title(f"1D, {N_simulations} sims per datapoint, tau={tau}")
    plt.xlabel("System size L")
    plt.ylabel("Number of sequential waves until eliminated (normalized)")
    plt.legend()
    plt.grid()
    plt.savefig(f"src/percolationPowerLawPopulations/plots/systemSize/1D/tau_{tau}.png", dpi=300)
    plt.show()


if __name__ == "__main__":
    main()