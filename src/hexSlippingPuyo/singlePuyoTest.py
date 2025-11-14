import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import RegularPolygon, Patch
import os

def main():
    L = 5
    H = 20  # Height
    lattice = np.zeros((H, L), dtype=int)
    
    # # Pick a random position
    # random_row = np.random.randint(0, H)
    # random_col = np.random.randint(0, L)
    random_row = 11
    random_col = 2
    lattice[random_row, random_col] = 1  # Red
    
    # Define directions based on row parity
    if random_row % 2 == 0:
        directions = [(-1, 0), (-1, -1), (0, -1), (0, 1), (1, 0), (1, -1)]
    else:
        directions = [(-1, 0), (-1, 1), (0, -1), (0, 1), (1, 0), (1, 1)]
    
    # Color neighbors yellow
    for dr, dc in directions:
        nr = (random_row + dr) % H
        nc = (random_col + dc) % L
        if lattice[nr, nc] == 0:
            lattice[nr, nc] = 2  # Yellow
    
    colors = mcolors.ListedColormap(['white', 'red', 'yellow'] + list(mcolors.TABLEAU_COLORS.values()))
    
    frames_dir = os.path.join(os.path.dirname(__file__), "plots", "singlePuyo", "frames")
    os.makedirs(frames_dir, exist_ok=True)
    
    # Save the plot
    fig, ax = plt.subplots()
    for r in range(H):
        for c in range(L):
            x = c + (r % 2) * 0.5
            y = (H - 1 - r) * (3**0.5 / 2)
            color = colors(lattice[r, c])
            hex_patch = RegularPolygon((x, y), 6, radius=0.5, facecolor=color, edgecolor='black', linewidth=0.5)
            ax.add_patch(hex_patch)
    ax.set_aspect('equal')
    ax.autoscale()
    ax.set_xlim(-0.5, L + 0.5)
    ax.set_ylim(-0.5, (H - 1) * (3**0.5 / 2) + 0.5)
    plt.title("Hex Neighborhood Test")
    
    frame_path = os.path.join(frames_dir, "neighborhood_test.png")
    plt.savefig(frame_path, dpi=150)
    plt.close(fig)

if __name__ == "__main__":
    main()
