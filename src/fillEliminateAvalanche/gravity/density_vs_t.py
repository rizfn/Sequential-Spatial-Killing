import numpy as np
import matplotlib.pyplot as plt
import cc3d
from tqdm import tqdm

def gravity(lattice):
    for col in range(lattice.shape[1]):
        col_vals = lattice[:, col]
        nonzero = col_vals[col_vals != 0]
        lattice[:, col] = 0
        if len(nonzero) > 0:
            lattice[-len(nonzero):, col] = nonzero

def merge_periodic_x_clusters(lattice, labels):
    L = lattice.shape[0]
    changed = True
    while changed:
        changed = False
        for y in range(L):
            left = (y, 0)
            right = (y, L-1)
            color = lattice[left]
            if color == 0:
                continue
            label_left = labels[left]
            label_right = labels[right]
            if label_left != 0 and label_right != 0 and color == lattice[right] and label_left != label_right:
                labels[labels == label_right] = label_left
                changed = True
    return labels

def main():
    L = 128
    N_list = [2, 3, 4, 5, 6, 7, 8]
    N_steps = 100
    densities = np.zeros((len(N_list), N_steps))

    for i, N in enumerate(tqdm(N_list)):
        lattice = np.zeros((L, L), dtype=int)
        for step in range(N_steps):
            # Fill all empty sites with random colors
            mask = lattice == 0
            lattice[mask] = np.random.randint(1, N+1, size=np.count_nonzero(mask))

            # Elimination-fall cycle
            while True:
                labels = cc3d.connected_components(lattice, connectivity=4)
                labels = merge_periodic_x_clusters(lattice, labels)
                eliminated = False
                for unique_label in np.unique(labels):
                    if unique_label != 0:
                        cluster = (labels == unique_label)
                        if np.sum(cluster) > 1:
                            lattice[cluster] = 0
                            eliminated = True
                gravity(lattice)
                if not eliminated:
                    break

            densities[i, step] = np.count_nonzero(lattice) / (L * L)

        plt.plot(range(N_steps), densities[i], label=f'N={N}')

    plt.xlabel('Step')
    plt.ylabel('Density')
    plt.legend()
    plt.grid()
    plt.savefig(f'src/fillEliminateAvalanche/gravity/plots/density/L_{L}.png', dpi=300)
    plt.show()


if __name__ == "__main__":
    main()