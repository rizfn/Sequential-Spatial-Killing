import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import cc3d

def main():
    L = 64
    N = 6
    N_steps = 100
    lattice = np.zeros((L, L), dtype=int)

    cmap = ListedColormap([(1,1,1,0)] + [plt.get_cmap('tab20', N)(i) for i in range(N)])


    for step in range(N_steps):
        random_species = np.random.randint(1, N+1, size=(L, L))
        mask = lattice == 0
        lattice[mask] = random_species[mask]
        # plt.imshow(lattice, cmap=cmap, interpolation='nearest', vmin=0, vmax=N)
        # plt.pause(1)
        # plt.cla()

        labels = cc3d.connected_components(lattice, connectivity=4, periodic_boundary=True)
        for unique_label in np.unique(labels):
            if unique_label != 0:
                cluster = (labels == unique_label)
                if np.sum(cluster) > 1:
                    lattice[cluster] = 0

        plt.imshow(lattice, cmap=cmap, interpolation='nearest', vmin=0, vmax=N)
        plt.pause(1)
        plt.cla()




if __name__ == "__main__":
    main()