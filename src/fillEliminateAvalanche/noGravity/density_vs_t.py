import numpy as np
import matplotlib.pyplot as plt
import cc3d
from tqdm import tqdm

def main():
    L = 128
    N_list = [2, 3, 4, 5, 6, 7, 8]
    N_steps = 100
    densities = np.zeros((len(N_list), N_steps))

    for N in tqdm(N_list):
        lattice = np.zeros((L, L), dtype=int)

        for step in range(N_steps):
            random_species = np.random.randint(1, N+1, size=(L, L))
            mask = lattice == 0
            lattice[mask] = random_species[mask]

            labels = cc3d.connected_components(lattice, connectivity=4, periodic_boundary=True)
            for unique_label in np.unique(labels):
                if unique_label != 0:
                    cluster = (labels == unique_label)
                    if np.sum(cluster) > 1:
                        lattice[cluster] = 0
            densities[N_list.index(N), step] = np.count_nonzero(lattice) / (L * L)

        plt.plot(range(N_steps), densities[N_list.index(N)], label=f'N={N}')

    plt.xlabel('Step')
    plt.ylabel('Density')
    plt.legend()
    plt.grid()
    plt.savefig(f'src/fillEliminateAvalanche/noGravity/plots/density/L_{L}.png', dpi=300)
    plt.show()


if __name__ == "__main__":
    main()