import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import cc3d


def place_puyo(lattice, column, species, L):
    """
    Place a puyo in the specified column. The puyo falls to the bottom-most available cell.
    """
    for row in range(L - 1, -1, -1):  # Start from the bottom row and move up
        if lattice[row, column] == 0:  # If the cell is empty
            lattice[row, column] = species
            break


def fall(lattice, L):
    """
    Simulate gravity: all cells fall until there are no empty spaces below them.
    After each fall, check for clusters and remove them iteratively until no more clusters exist.
    """
    while True:
        moved = False
        # Simulate gravity
        for row in range(L - 1, 0, -1):  # Start from the second-to-last row and move up
            for col in range(L):
                if lattice[row, col] == 0 and lattice[row - 1, col] > 0:  # Empty cell below, occupied cell above
                    lattice[row, col] = lattice[row - 1, col]  # Move the cell down
                    lattice[row - 1, col] = 0  # Empty the original cell
                    moved = True
        remove_puyo(lattice)
        if not moved:
            break


def remove_puyo(lattice):
    """
    Remove clusters of connected puyos. A cluster is removed if it has more than one cell.
    """
    labels, N = cc3d.connected_components(lattice, connectivity=4, return_N=True)
    for label in range(1, N + 1):
        if np.sum(labels == label) > 1:  # If the cluster has more than one cell
            lattice[labels == label] = 0  # Remove the cluster


def main():
    L = 128
    N_species_list = np.arange(2, 9, 1)
    N_puyos = int(L**2 / 2)

    masses = np.zeros((len(N_species_list), N_puyos))
    heights = np.zeros((len(N_species_list), N_puyos))

    for i, N_species in enumerate(N_species_list):
        lattice = np.zeros((L, L), dtype=int)
        random_columns = np.random.randint(0, L, size=N_puyos)
        random_species = np.random.randint(1, N_species + 1, size=N_puyos)

        for step in tqdm(range(N_puyos)):

            place_puyo(lattice, random_columns[step], random_species[step], L)
            fall(lattice, L)

            masses[i, step] = np.sum(lattice > 0)
            heights[i, step] = np.max(L - np.nonzero(lattice)[0]) if np.any(lattice > 0) else 0

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    for i, N_species in enumerate(N_species_list):
        ax1.plot(np.arange(N_puyos), masses[i], label=f'{N_species} species')
        ax2.plot(np.arange(N_puyos), heights[i], label=f'{N_species} species')
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Mass')
    ax1.grid()
    ax1.legend()
    ax2.set_xlabel('Time')
    ax2.set_ylabel('Max Height')
    ax2.grid()
    ax2.legend()
    plt.savefig('src/puyopuyo/plots/gravity/2D_MassVsTime.png', dpi=300)
    plt.show()


if __name__ == "__main__":
    main()