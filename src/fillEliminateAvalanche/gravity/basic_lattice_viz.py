import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import cc3d

def gravity(lattice):
    # For each column, move nonzero values to the bottom
    for col in range(lattice.shape[1]):
        col_vals = lattice[:, col]
        nonzero = col_vals[col_vals != 0]
        lattice[:, col] = 0
        if len(nonzero) > 0:
            lattice[-len(nonzero):, col] = nonzero

def merge_periodic_x_clusters(lattice, labels):
    """Merge clusters that are connected across the x-boundary and have the same color."""
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
                # Merge label_right into label_left
                labels[labels == label_right] = label_left
                changed = True
    return labels

def main():
    L = 64
    N = 7
    N_steps = 100
    lattice = np.zeros((L, L), dtype=int)
    cmap = ListedColormap([(1,1,1,0)] + [plt.get_cmap('tab20', N)(i) for i in range(N)])

    for step in range(N_steps):
        # Fill all empty sites with random colors
        mask = lattice == 0
        lattice[mask] = np.random.randint(1, N+1, size=np.count_nonzero(mask))

        plt.imshow(lattice, cmap=cmap, interpolation='nearest', vmin=0, vmax=N)
        plt.pause(1)
        plt.cla()

        while True:
            # Eliminate clusters
            labels = cc3d.connected_components(lattice, connectivity=4)
            labels = merge_periodic_x_clusters(lattice, labels)
            eliminated = False
            for unique_label in np.unique(labels):
                if unique_label != 0:
                    cluster = (labels == unique_label)
                    if np.sum(cluster) > 1:
                        lattice[cluster] = 0
                        eliminated = True

            plt.imshow(lattice, cmap=cmap, interpolation='nearest', vmin=0, vmax=N)
            plt.pause(1)
            plt.cla()

            # Apply gravity
            gravity(lattice)

            # Show the current state
            plt.imshow(lattice, cmap=cmap, interpolation='nearest', vmin=0, vmax=N)
            plt.pause(1)
            plt.cla()

            # If no more eliminations, break
            if not eliminated:
                break

if __name__ == "__main__":
    main()