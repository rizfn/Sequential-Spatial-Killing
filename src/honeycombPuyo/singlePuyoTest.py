import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Rectangle, Patch
from tqdm import tqdm
import collections
import os


def get_neighbors(row, col, L, H):
    neighbors = []

    for dr in (-1, 1):
        nr = row + dr
        if 0 <= nr < H:
            neighbors.append((nr, col))

    if col % 2 == 0:
        side_rows = (row, row + 1)
    else:
        side_rows = (row - 1, row)

    for dc in (-1, 1):
        nc = (col + dc + L) % L
        for nr in side_rows:
            if 0 <= nr < H:
                neighbors.append((nr, nc))

    return neighbors


def place_puyo(lattice, moved_sites, col, species, L, H):
    for row in range(H - 1, -1, -1):
        if lattice[row, col] == 0:
            lattice[row, col] = species
            moved_sites[row, col] = True
            return True
    return False

def annihilate_puyo(lattice, moved_sites, L, H):
    visited = np.zeros((H, L), dtype=bool)
    
    for row in range(H):
        for col in range(L):
            if moved_sites[row, col] and not visited[row, col] and lattice[row, col] != 0:
                cluster = []
                q = collections.deque()
                q.append((row, col))
                visited[row, col] = True
                species = lattice[row, col]

                while q:
                    r, c = q.popleft()
                    cluster.append((r, c))
                    
                    for nr, nc in get_neighbors(r, c, L, H):
                        if not visited[nr, nc] and lattice[nr, nc] == species:
                            q.append((nr, nc))
                            visited[nr, nc] = True
                
                if len(cluster) > 1:
                    for r, c in cluster:
                        lattice[r, c] = 0
                        moved_sites[r, c] = False

def fall(lattice, moved_sites, L, H):
    new_moved = np.zeros((H, L), dtype=bool)

    for col in range(L):
        write_row = H - 1
        for row in range(H - 1, -1, -1):
            if lattice[row, col] != 0:
                if write_row != row:
                    lattice[write_row, col] = lattice[row, col]
                    lattice[row, col] = 0
                    new_moved[write_row, col] = True
                write_row -= 1

    moved_sites[:] = new_moved

def main():
    L = 6
    H = 20  # Height, matching C++ style
    N_species = 5
    steps = 20  # Number of steps, each adding L puyos
    total_drops = steps * L
    
    lattice = np.zeros((H, L), dtype=int)
    moved_sites = np.zeros((H, L), dtype=bool)
    random_columns = np.random.randint(0, L, size=(steps, L))
    random_species = np.random.randint(1, N_species + 1, size=(steps, L))
    
    colors = mcolors.ListedColormap(['white'] + list(mcolors.TABLEAU_COLORS.values()))
    
    frames_dir = os.path.join(os.path.dirname(__file__), "plots", "singlePuyo", "frames")
    os.makedirs(frames_dir, exist_ok=True)
    
    for drop in tqdm(range(total_drops)):
        step = drop // L
        i = drop % L
        col = random_columns[step, i]
        species = random_species[step, i]
        placed = place_puyo(lattice, moved_sites, col, species, L, H)
        if not placed:
            continue
        
        # Annihilation-fall cycle
        while True:
            annihilate_puyo(lattice, moved_sites, L, H)
            new_moved = np.zeros((H, L), dtype=bool)
            fall(lattice, new_moved, L, H)
            if np.array_equal(new_moved, moved_sites):
                break
            moved_sites[:] = new_moved
        
        # Save frame for each drop
        fig, ax = plt.subplots()
        for r in range(H):
            for c in range(L):
                if lattice[r, c] != 0:
                    x = c
                    y = (H - 1 - r) + (0.5 if c % 2 == 1 else 0.0)
                    color = colors(lattice[r, c])
                    square = Rectangle((x - 0.5, y - 0.5), 1.0, 1.0, facecolor=color, edgecolor='black', linewidth=0.5)
                    ax.add_patch(square)
        ax.set_aspect('equal')
        ax.autoscale()
        ax.set_xlim(-0.5, L - 0.5)
        ax.set_ylim(-0.5, H - 0.5 + 0.5)
        plt.title(f"Drop: {drop}")
        
        # Add legend for just dropped puyo
        just_dropped_color = colors(species)
        legend_patch = Patch(color=just_dropped_color, label='Just dropped')
        ax.legend(handles=[legend_patch], loc='upper right')
        
        frame_path = os.path.join(frames_dir, f"frame_{drop:04d}.png")
        plt.savefig(frame_path, dpi=150)
        plt.close(fig)

if __name__ == "__main__":
    main()
