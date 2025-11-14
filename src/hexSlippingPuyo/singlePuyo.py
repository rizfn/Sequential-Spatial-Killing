import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import RegularPolygon, Patch
from tqdm import tqdm
import collections
import os

def place_puyo(lattice, moved_sites, col, species, L, H):
    # First, select random col as before
    # Then, check col-1, col, col+1 for lowest drop height
    candidates = []
    min_height = float('inf')
    for dc in [-1, 0, 1]:
        c = (col + dc) % L
        # Find the drop row for this column
        drop_row = -1
        for row in range(H - 1, -1, -1):
            if lattice[row][c] == 0:
                drop_row = row
                break
        if drop_row != -1:
            height = H - drop_row  # Lower height means higher in lattice (lower row number)
            if height < min_height:
                min_height = height
                candidates = [c]
            elif height == min_height:
                candidates.append(c)
    
    # If no candidates (lattice full), do nothing or handle, but assume not
    if not candidates:
        return
    
    # Choose randomly among candidates
    chosen_col = np.random.choice(candidates)
    
    # Now place in chosen_col
    for row in range(H - 1, -1, -1):
        if lattice[row][chosen_col] == 0:
            lattice[row][chosen_col] = species
            moved_sites[row][chosen_col] = True
            return

def annihilate_puyo(lattice, moved_sites, L, H):
    visited = np.zeros((H, L), dtype=bool)
    
    for row in range(H):
        for col in range(L):
            if lattice[row, col] != 0 and not visited[row, col]:  # Removed moved_sites check; scan all puyos
                # Flood fill
                cluster = []
                q = collections.deque()
                q.append((row, col))
                visited[row, col] = True
                species = lattice[row, col]
                
                # Define directions based on row parity
                if row % 2 == 0:
                    directions = [(-1, 0), (-1, -1), (0, -1), (0, 1), (1, 0), (1, -1)]
                else:
                    directions = [(-1, 0), (-1, 1), (0, -1), (0, 1), (1, 0), (1, 1)]
                
                while q:
                    r, c = q.popleft()
                    cluster.append((r, c))
                    
                    for dr, dc in directions:
                        nr = (r + dr) % H
                        nc = (c + dc) % L
                        if not visited[nr, nc] and lattice[nr, nc] == species:
                            q.append((nr, nc))
                            visited[nr, nc] = True
                
                if len(cluster) > 1:
                    for r, c in cluster:
                        lattice[r, c] = 0

def fall(lattice, moved_sites, L, H):
    new_moved = np.zeros((H, L), dtype=bool)
    changed = True
    while changed:
        changed = False
        for r in range(H - 2, -1, -1):
            for c in range(L):
                if lattice[r, c] != 0:
                    # Only allow falling to the cell 2 rows below in the same column
                    br2 = r + 2
                    if br2 < H and lattice[br2, c] == 0:
                        lattice[br2, c] = lattice[r, c]
                        lattice[r, c] = 0
                        new_moved[br2, c] = True
                        changed = True
                        continue
    moved_sites[:] = new_moved

def main():
    L = 5
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
        place_puyo(lattice, moved_sites, col, species, L, H)
        
        # Annihilation-fall cycle
        while True:
            annihilate_puyo(lattice, moved_sites, L, H)
            old_moved = moved_sites.copy()
            fall(lattice, moved_sites, L, H)
            if np.array_equal(moved_sites, old_moved):
                break
        
        # Save frame for each drop
        fig, ax = plt.subplots()
        for r in range(H):
            for c in range(L):
                if lattice[r, c] != 0:
                    x = c + (r % 2) * 0.5
                    y = (H - 1 - r) * (3**0.5 / 2)  # Flip y-axis so bottom is at bottom
                    color = colors(lattice[r, c])
                    hex_patch = RegularPolygon((x, y), 6, radius=0.5, facecolor=color, edgecolor='black', linewidth=0.5)
                    ax.add_patch(hex_patch)
        ax.set_aspect('equal')
        ax.autoscale()
        ax.set_xlim(-0.5, L + 0.5)
        ax.set_ylim(-0.5, (H - 1) * (3**0.5 / 2) + 0.5)
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
