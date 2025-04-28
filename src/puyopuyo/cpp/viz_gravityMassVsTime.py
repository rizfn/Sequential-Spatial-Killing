import numpy as np
import matplotlib.pyplot as plt

def viz_1D():
    N_list = np.arange(2, 21)
    steps = 1024

    fig, ax1 = plt.subplots(1, 1, figsize=(7, 6))

    ax1.plot([0, steps], [1, (steps+1)], color="grey", linestyle="--", label="Total mass added")

    # Use a colormap (e.g., rainbow) to assign colors based on N_species
    colormap = plt.get_cmap("rainbow", len(N_list))

    for i, N in enumerate(N_list):
        step, mass = np.loadtxt(f"src/puyopuyo/cpp/outputs/gravity1D/N_{N}_steps_{steps}.tsv", delimiter="\t", skiprows=1, unpack=True)
        color = colormap(i / len(N_list))  # Normalize the index to [0, 1] for the colormap
        ax1.plot(step, mass, label=f"N={N}", color=color)

    ax1.set_xlabel("Time")
    ax1.set_ylabel("Mass")
    ax1.grid()
    ax1.legend(ncols=2)

    plt.savefig(f"src/puyopuyo/cpp/plots/1D_gravityMassVsTime.png", dpi=300)
    plt.show()


def viz_2D():
    L = 128
    N_list = np.arange(2, 21)
    steps = 1024

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    ax1.plot([0, steps], [L, L * (steps+1)], color="grey", linestyle="--", label="Total mass added")
    ax2.plot([0, steps], [1, steps+1], color="grey", linestyle="--", label="Average (dense) height")

    # Use a colormap (e.g., rainbow) to assign colors based on N_species
    colormap = plt.get_cmap("rainbow", len(N_list))

    for i, N in enumerate(N_list):
        step, mass, height = np.loadtxt(f"src/puyopuyo/cpp/outputs/gravity2D/L_{L}_N_{N}_steps_{steps}.tsv", delimiter="\t", skiprows=1, unpack=True)
        color = colormap(i / len(N_list))  # Normalize the index to [0, 1] for the colormap
        ax1.plot(step, mass, label=f"N={N}", color=color)
        ax2.plot(step, height, label=f"N={N}", color=color)

    ax1.set_xlabel("Time")
    ax1.set_ylabel("Mass")
    ax1.grid()
    ax1.legend(ncols=2)
    ax2.set_xlabel("Time")
    ax2.set_ylabel("Max Height")
    ax2.grid()
    ax2.legend(ncols=2)

    plt.savefig(f"src/puyopuyo/cpp/plots/2D_gravityMassVsTime_L_{L}.png", dpi=300)
    plt.show()


def viz_3D():
    L = 24
    N_list = np.arange(2, 21)
    steps = 1024

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))


    # Use a colormap (e.g., rainbow) to assign colors based on N_species
    colormap = plt.get_cmap("rainbow", len(N_list))

    for i, N in enumerate(N_list):
        step, mass, height = np.loadtxt(f"src/puyopuyo/cpp/outputs/gravity3D/L_{L}_N_{N}_steps_{steps}.tsv", delimiter="\t", skiprows=1, unpack=True)
        color = colormap(i / len(N_list))  # Normalize the index to [0, 1] for the colormap
        ax1.plot(step, mass, label=f"N={N}", color=color)
        ax2.plot(step, height, label=f"N={N}", color=color)

    ax1_ylim = ax1.get_ylim()
    ax1.plot([0, steps], [L*L, L*L * (steps+1)], color="grey", linestyle="--", label="Total mass added")
    ax1.set_ylim(ax1_ylim)
    ax2_ylim = ax2.get_ylim()
    ax2.plot([0, steps], [1, steps+1], color="grey", linestyle="--", label="Average (dense) height")
    ax2.set_ylim(ax2_ylim)

    ax1.set_xlabel("Time")
    ax1.set_ylabel("Mass")
    ax1.grid()
    ax1.legend(ncols=2)
    ax2.set_xlabel("Time")
    ax2.set_ylabel("Max Height")
    ax2.grid()
    ax2.legend(ncols=2)

    plt.savefig(f"src/puyopuyo/cpp/plots/3D_gravityMassVsTime_L_{L}.png", dpi=300)
    plt.show()

def viz_4D():
    L = 8
    N_list = np.arange(2, 21)
    steps = 1024

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))


    # Use a colormap (e.g., rainbow) to assign colors based on N_species
    colormap = plt.get_cmap("rainbow", len(N_list))

    for i, N in enumerate(N_list):
        step, mass, height = np.loadtxt(f"src/puyopuyo/cpp/outputs/gravity4D/L_{L}_N_{N}_steps_{steps}.tsv", delimiter="\t", skiprows=1, unpack=True)
        color = colormap(i / len(N_list))  # Normalize the index to [0, 1] for the colormap
        ax1.plot(step, mass, label=f"N={N}", color=color)
        ax2.plot(step, height, label=f"N={N}", color=color)

    ax1_ylim = ax1.get_ylim()
    ax1.plot([0, steps], [L*L*L, L*L*L * (steps+1)], color="grey", linestyle="--", label="Total mass added")
    ax1.set_ylim(ax1_ylim)
    ax2_ylim = ax2.get_ylim()
    ax2.plot([0, steps], [1, steps+1], color="grey", linestyle="--", label="Average (dense) height")
    ax2.set_ylim(ax2_ylim)

    ax1.set_xlabel("Time")
    ax1.set_ylabel("Mass")
    ax1.grid()
    ax1.legend(ncols=2)
    ax2.set_xlabel("Time")
    ax2.set_ylabel("Max Height")
    ax2.grid()
    ax2.legend(ncols=2)

    plt.savefig(f"src/puyopuyo/cpp/plots/4D_gravityMassVsTime_L_{L}.png", dpi=300)
    plt.show()



def drift_vs_N():
    N_list = np.arange(2, 21)
    steps = 1024
    dims = [1, 2, 3, 4]
    Ls = {1: None, 2: 128, 3: 24, 4: 8}

    interval = 100
    n_bins = steps // interval

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()

    for ax, dim in zip(axes, dims):
        slopes = []
        slope_errs = []

        for N in N_list:
            # load data
            if dim == 1:
                fn = f"src/puyopuyo/cpp/outputs/gravity1D/N_{N}_steps_{steps}.tsv"
                step, mass = np.loadtxt(fn, delimiter="\t", skiprows=1, unpack=True)
            else:
                L = Ls[dim]
                fn = f"src/puyopuyo/cpp/outputs/gravity{dim}D/L_{L}_N_{N}_steps_{steps}.tsv"
                data = np.loadtxt(fn, delimiter="\t", skiprows=1, unpack=True)
                step, mass = data[0], data[1]

            # bin into intervals and compute mean mass
            t_bins = []
            m_means = []
            for b in range(n_bins):
                start = b * interval
                end = min((b+1) * interval, len(step))
                t_slice = step[start:end]
                m_slice = mass[start:end]
                if len(m_slice) == 0:
                    continue
                t_bins.append(np.mean(t_slice))
                m_means.append(np.mean(m_slice))

            t_bins = np.array(t_bins)
            m_means = np.array(m_means)

            # linear fit on the binned means
            slope, intercept = np.polyfit(t_bins, m_means, 1)

            # estimate error on slope from residuals
            m_pred = slope * t_bins + intercept
            resid = m_means - m_pred
            slope_err = np.std(resid) / np.sqrt(np.sum((t_bins - np.mean(t_bins))**2))

            slopes.append(slope)
            slope_errs.append(slope_err)

        # plot drift vs N with error bars
        ax.errorbar(N_list, slopes, yerr=slope_errs,
                    fmt='o-', capsize=4)
        ax.set_title(f"{dim}D drift vs N")
        ax.set_xlabel("Number of species N")
        ax.set_ylabel("Drift (mass/time)")
        ax.grid()

        if dim == 1:
            x_vals = np.linspace(N_list[0], N_list[-1], 100)
            ax.plot(x_vals, 1-2/x_vals, color="grey", linestyle="--", label="1-2/N")

    plt.tight_layout()
    plt.savefig(f"src/puyopuyo/cpp/plots/drift_vs_N.png", dpi=300)
    plt.show()


if __name__ == "__main__":
    # viz_1D()
    # viz_2D()
    # viz_3D()
    # viz_4D()
    drift_vs_N()