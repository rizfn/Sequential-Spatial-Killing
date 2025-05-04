import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.colors import ListedColormap

def load_lattice_tsv(tsv_path, L):
    df = pd.read_csv(tsv_path, sep="\t")
    lattices = []
    for lattice_str in df['lattice']:
        flat = np.array(list(map(int, lattice_str.split(','))))
        lattices.append(flat.reshape((L, L)))
    steps = df['step'].values
    return steps, lattices

def get_custom_cmap(n_colors):
    tab20 = plt.get_cmap('tab20', n_colors)
    colors = [ (0,0,0,1) ]
    for i in range(1, n_colors):
        colors.append(tab20(i))
    return ListedColormap(colors)

def show_lattice_animation(steps, lattices, interval=0.05):
    plt.ion()
    fig, ax = plt.subplots(figsize=(6,6))
    vmax = np.max(lattices)
    cmap = get_custom_cmap(int(vmax)+1)
    im = ax.imshow(lattices[0], cmap=cmap, vmin=0, vmax=vmax)
    ax.set_title(f"Step: {steps[0]}")
    ax.axis('off')
    plt.show()

    for i in range(len(lattices)):
        ax.cla()
        im = ax.imshow(lattices[i], cmap=cmap, vmin=0, vmax=vmax)
        ax.set_title(f"Step: {steps[i]}")
        ax.axis('off')
        plt.pause(interval)
    plt.ioff()
    plt.show()

if __name__ == "__main__":
    L = 32
    N = 3
    steps = 1024
    tsv_path = f"src/edenPuyoPuyo/outputs/lattice2D/L_{L}_N_{N}_steps_{steps}.tsv"

    steps, lattices = load_lattice_tsv(tsv_path, L)
    show_lattice_animation(steps, lattices, interval=0.05)