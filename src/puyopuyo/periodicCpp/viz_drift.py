import numpy as np
import matplotlib.pyplot as plt

def drift_vs_inverseN():
    N_list = np.arange(2, 21)
    steps = 1024
    dims = [1, 2, 3, 4]
    Ls = {1: None, 2: 128, 3: 24, 4: 8}
    tolerance = 0.01  # tolerance for drift to be considered nonzero

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
                normalized_mass = mass
            else:
                L = Ls[dim]
                fn = f"src/puyopuyo/periodicCpp/outputs/gravity{dim}D/L_{L}_N_{N}_steps_{steps}.tsv"
                data = np.loadtxt(fn, delimiter="\t", skiprows=1, unpack=True)
                step, mass = data[0], data[1]
                normalized_mass = mass / (L ** (dim - 1))  # as L^(d-1) puyos are added each step

            # bin into intervals and compute mean normalized mass
            t_bins = []
            m_means = []
            for b in range(n_bins):
                start = b * interval
                end = min((b+1) * interval, len(step))
                t_slice = step[start:end]
                m_slice = normalized_mass[start:end]
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

        # Plot drift vs 1/N with error bars
        inv_N = 1 / N_list
        ax.errorbar(inv_N, slopes, yerr=slope_errs, fmt='o-', capsize=4, label="All points")
        ax.set_title(f"{dim}D Drift vs 1/N")
        ax.set_xlabel("1/N")
        ax.set_ylabel("Normalized Drift (mass/column/time)")
        ax.grid()

        # Highlight fit points (drift > tolerance) in purple
        fit_mask = np.array(slopes) > tolerance
        if np.any(fit_mask):
            ax.errorbar(inv_N[fit_mask], np.array(slopes)[fit_mask], yerr=np.array(slope_errs)[fit_mask],
                        fmt='o', color='purple', capsize=4, label="Fit points")

            x_fit = inv_N[fit_mask]
            y_fit = np.array(slopes)[fit_mask]
            y_err = np.array(slope_errs)[fit_mask]

            # Weighted linear least squares (chi-square fit)
            w = 1 / (y_err ** 2)
            W = np.sum(w)
            Wx = np.sum(w * x_fit)
            Wy = np.sum(w * y_fit)
            Wxx = np.sum(w * x_fit * x_fit)
            Wxy = np.sum(w * x_fit * y_fit)
            Delta = W * Wxx - Wx ** 2

            slope = (W * Wxy - Wx * Wy) / Delta
            intercept = (Wxx * Wy - Wx * Wxy) / Delta

            fit_line = lambda x: slope * x + intercept

            x_crit = -intercept / slope
            x_fit_full = np.linspace(x_fit.min(), x_crit, 200)
            ax.plot(x_fit_full, fit_line(x_fit_full), 'r--', label=f"Fit m = {slope:.3f}")
            ax.axvline(x_crit, color='k', linestyle=':', label=f"Critical 1/N = {x_crit:.3f}")
            ax.set_title(f"{dim}D: N$_\\text{{c}}$ = {1/x_crit:.2f}")

            ax.legend()

    plt.tight_layout()
    plt.savefig(f"src/puyopuyo/periodicCpp/plots/drift/driftVsInvN.png", dpi=300)
    plt.show()

if __name__ == "__main__":
    drift_vs_inverseN()