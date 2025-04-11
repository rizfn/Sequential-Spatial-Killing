import numpy as np
import matplotlib.pyplot as plt
import cc3d
from tqdm import tqdm


def update(matrix, currentVictim):
    # Make a copy of the matrix, with 0 set to currentVictim
    copiedMatrix = matrix.copy()
    copiedMatrix[copiedMatrix == 0] = currentVictim

    # Find the connected components of the copied matrix
    labels = cc3d.connected_components(copiedMatrix, connectivity=4)

    # Look for the cluster at the boundaries of 0
    boundaryCluster = labels[0, 0]

    # Update the matrix
    matrix[labels == boundaryCluster] = 0

    return matrix


def main():
    N_species_values = [2, 3, 4, 8]  # Different numbers of species (subplots)
    L_array = [20, 40, 80, 160]  # Different system sizes (lines in the same plot)
    max_steps = 1000
    N_simulations = 20  # Number of simulations for each (L, N_species)

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))  # Create subplots
    axes = axes.flatten()

    for idx, N_species in enumerate(N_species_values):
        ax = axes[idx]
        for L in L_array:
            fractions_over_time = []
            max_time_steps = 0  # Track the maximum number of steps for this system size

            for _ in range(N_simulations):
                matrix = np.random.randint(1, N_species + 1, (L, L))
                matrix = np.pad(matrix, pad_width=1, mode='constant', constant_values=0)
                victimOrder = np.arange(1, N_species + 1)

                fractions = []
                for t in range(max_steps):
                    currentVictim = victimOrder[t % N_species]
                    matrix = update(matrix, currentVictim)
                    fraction = np.count_nonzero(matrix) / (L * L)  # Fraction of surviving elements
                    fractions.append(fraction)

                    if np.count_nonzero(matrix) == 0:
                        break

                max_time_steps = max(max_time_steps, len(fractions))
                fractions_over_time.append(fractions)

            # Truncate or pad all simulations to the maximum time steps
            for i in range(len(fractions_over_time)):
                fractions_over_time[i] = fractions_over_time[i][:max_time_steps]
                if len(fractions_over_time[i]) < max_time_steps:
                    fractions_over_time[i] += [0] * (max_time_steps - len(fractions_over_time[i]))

            # Compute mean and standard deviation over simulations
            fractions_over_time = np.array(fractions_over_time)
            mean_fractions = np.mean(fractions_over_time, axis=0)
            std_fractions = np.std(fractions_over_time, axis=0)

            # Plot the mean line and shaded error region
            time_steps = range(max_time_steps)
            ax.plot(time_steps, mean_fractions, label=f"L={L}")
            ax.fill_between(
                time_steps,
                mean_fractions - std_fractions,
                mean_fractions + std_fractions,
                alpha=0.2,
            )

        ax.set_title(f"{N_species} species")
        ax.set_xlabel("Time (steps)")
        ax.set_ylabel("Fraction of survivors at given time (cumulative)")
        ax.legend()
        ax.grid()

    plt.tight_layout()
    plt.savefig("src/percolation/plots/fractionKilled/2DCumulativeSurvivors.png", dpi=300)
    plt.show()


if __name__ == "__main__":
    main()