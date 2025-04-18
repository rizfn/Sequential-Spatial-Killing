import numpy as np
import matplotlib.pyplot as plt
import cc3d
from tqdm import tqdm

def place_puyo(lattice, column, species, L):
    for row in range(L-1, 0, -1):
        if lattice[row-1, column] != 0:
            lattice[row, column] = species
            break
        if lattice[row, (column-1) % L] != 0:
            lattice[row, column] = species
            break
        if lattice[row, (column+1) % L] != 0:
            lattice[row, column] = species
            break
    else:
        lattice[0, column] = species

def remove_puyo(lattice):
    labels, N = cc3d.connected_components(lattice, connectivity=4, return_N=True)
    for label in range(1, N+1):
        if np.sum(labels == label) > 1:
            lattice[labels == label] = 0

def fall(lattice, L):
    while True:
        puyoExists = lattice > 0
        labels, N = cc3d.connected_components(puyoExists, connectivity=4, return_N=True)
        labels_at_bottom = np.unique(labels[0, :])
        floating_labels = list(set(np.arange(1, N+1)) - set(labels_at_bottom))  # Find floating clusters
        if len(floating_labels) == 0:  # No floating clusters
            break
        floating_mask = np.isin(labels, floating_labels)  # Mask for floating clusters
        floating_values = lattice[floating_mask]  # Extract the actual values of floating clusters
        lattice[floating_mask] = 0  # Remove floating clusters from their current positions
        shifted_mask = np.roll(floating_mask, shift=-1, axis=0)  # Shift the mask up by one row
        shifted_mask[-1, :] = False  # Ensure nothing wraps around to the bottom row
        lattice[shifted_mask] = floating_values  # Place the floating clusters one row down
        remove_puyo(lattice)  # Remove any clusters that are now connected after falling


def main():
    L = 128
    N_species_list = np.arange(2, 9, 1)
    N_puyos = int(L**2 / 2)

    masses = np.zeros((len(N_species_list), N_puyos))
    heights = np.zeros((len(N_species_list), N_puyos))

    for i, N_species in enumerate(N_species_list):
        lattice = np.zeros((L, L), dtype=int)
        random_columns = np.random.randint(0, L, size=N_puyos)
        random_species = np.random.randint(1, N_species+1, size=N_puyos)

        for step in tqdm(range(N_puyos)):

            place_puyo(lattice, random_columns[step], random_species[step], L)
            remove_puyo(lattice)
        
            fall(lattice, L)

            masses[i, step] = np.sum(lattice > 0)
            heights[i, step] = np.max(np.where(lattice > 0)[0]) + 1 if np.any(lattice > 0) else 0

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
    plt.savefig('src/puyopuyo/plots/ballisticDeposition/2D_MassVsTime.png', dpi=300)
    plt.show()


if __name__ == "__main__":
    main()