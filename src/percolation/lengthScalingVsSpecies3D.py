import numpy as np
import matplotlib.pyplot as plt
import cc3d
from scipy.optimize import curve_fit
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


def linear_model(L, a, b):
    """Linear model for fitting: y = a * L + b"""
    return a * L + b


def main():
    N_species_values = np.arange(4, 22, 2)  # Multiple values of N_species
    L_array = np.arange(10, 110, 10)
    max_steps = 10000
    N_simulations = 10  # Number of simulations for each (L, N_species)
    slopes = []  # To store the slopes for each N_species

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

        # Fit a line to the data using chi-squared minimization
        t_means = np.array(t_means)
        t_stds = np.array(t_stds)
        popt, pcov = curve_fit(
            linear_model, L_array, t_means, sigma=t_stds, absolute_sigma=True
        )
        slope, intercept = popt
        slopes.append(slope)

    # Plot slope vs N_species
    plt.plot(N_species_values, slopes, marker="o")
    plt.title("Slope vs Number of Species (3D)")
    plt.xlabel("Number of Species")
    plt.ylabel("Slope of Linear Fit")
    plt.grid()
    plt.savefig("src/percolation/plots/systemSize/3DScalingSlopeVsSpecies.png", dpi=300)
    plt.show()


if __name__ == "__main__":
    main()