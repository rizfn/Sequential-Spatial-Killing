import numpy as np
import matplotlib.pyplot as plt
import glob
import re
import collections

def mass_vs_time():
    steps = 1024 * 2
    L = 256

    files = []
    for file in glob.glob(f"src/probabilityPuyoPuyo/outputs/gravity2D/L_{L}_N_*_steps_{steps}.tsv"):
        N = float(file.split("_")[3])
        files.append((N, file))
    files.sort()  # Sort by N

    for N, file in files:
        step, mass, height = np.loadtxt(file, delimiter="\t", skiprows=1, unpack=True)
        normalized_mass = mass / L
        plt.plot(step, normalized_mass, label=f"N={N}")

    plt.tight_layout()
    plt.grid()
    plt.legend()
    plt.xlabel("Time")
    plt.ylabel("Normalized Mass (mass/column/time)")
    plt.savefig(f"src/probabilityPuyoPuyo/plots/drift/mass_L_{L}_steps_{steps}.png", dpi=300)
    plt.show()

def drift_vs_inverseN():
    N_list = np.array([5, 6, 6.01, 6.05, 6.1, 6.2, 6.3, 6.5, 6.75, 6.9, 7, 8, 9, 10, 11])
    steps = 1024 * 2
    dims = [2]
    Ls = {1: None, 2: 256, 3: 24, 4: 8}
    tolerance = 0.01  # tolerance for drift to be considered nonzero

    interval = 100
    n_bins = steps // interval

    fig, axes = plt.subplots(1, 2, figsize=(12, 10))
    axes = axes.flatten()

    for ax, dim in zip(axes, dims):
        slopes = []
        slope_errs = []

        for N in N_list:
            # load data
            if dim == 1:
                fn = f"src/probabilityPuyoPuyo/cpp/outputs/gravity1D/N_{N}_steps_{steps}.tsv"
                step, mass = np.loadtxt(fn, delimiter="\t", skiprows=1, unpack=True)
                normalized_mass = mass
            else:
                L = Ls[dim]
                if N == int(N):
                    fn = f"src/probabilityPuyoPuyo/outputs/gravity{dim}D/L_{L}_N_{int(N)}_steps_{steps}.tsv"
                else:
                    fn = f"src/probabilityPuyoPuyo/outputs/gravity{dim}D/L_{L}_N_{N}_steps_{steps}.tsv"
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
    plt.savefig(f"src/probabilityPuyoPuyo/plots/drift/driftVsInvN.png", dpi=300)
    plt.show()


def drift_vs_invEntropy():
    N_list = np.array([5, 6, 6.01, 6.05, 6.1, 6.2, 6.3, 6.5, 6.75, 6.9, 7, 8, 9, 10, 11])
    steps = 1024 * 2
    dims = [2]
    Ls = {1: None, 2: 256, 3: 24, 4: 8}
    tolerance = 0.01  # tolerance for drift to be considered nonzero

    interval = 100
    n_bins = steps // interval

    fig, axes = plt.subplots(1, 2, figsize=(12, 10))
    axes = axes.flatten()

    def entropy_for_N(N):
        n_int = int(np.floor(N))
        frac = N - n_int
        weights = [1.0] * n_int
        if frac > 0:
            weights.append(frac)
        weights = np.array(weights)
        probs = weights / weights.sum()
        S = -np.sum(probs * np.log(probs))
        return S

    for ax, dim in zip(axes, dims):
        slopes = []
        slope_errs = []
        entropies = []

        for N in N_list:
            # load data
            if dim == 1:
                fn = f"src/probabilityPuyoPuyo/cpp/outputs/gravity1D/N_{N}_steps_{steps}.tsv"
                step, mass = np.loadtxt(fn, delimiter="\t", skiprows=1, unpack=True)
                normalized_mass = mass
            else:
                L = Ls[dim]
                if N == int(N):
                    fn = f"src/probabilityPuyoPuyo/outputs/gravity{dim}D/L_{L}_N_{int(N)}_steps_{steps}.tsv"
                else:
                    fn = f"src/probabilityPuyoPuyo/outputs/gravity{dim}D/L_{L}_N_{N}_steps_{steps}.tsv"
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
            entropies.append(entropy_for_N(N))

        entropies = np.array(entropies)
        slopes = np.array(slopes)
        slope_errs = np.array(slope_errs)
        inv_entropies = 1 / entropies

        # Plot drift vs 1/entropy with error bars
        ax.errorbar(inv_entropies, slopes, yerr=slope_errs, fmt='o-', capsize=4, label="All points")
        ax.set_title(f"{dim}D Drift vs 1/Entropy")
        ax.set_xlabel("1 / Entropy S")
        ax.set_ylabel("Normalized Drift (mass/column/time)")
        ax.grid()

        # Highlight fit points (drift > tolerance) in purple
        fit_mask = slopes > tolerance
        if np.any(fit_mask):
            ax.errorbar(inv_entropies[fit_mask], slopes[fit_mask], yerr=slope_errs[fit_mask],
                        fmt='o', color='purple', capsize=4, label="Fit points")

            x_fit = inv_entropies[fit_mask]
            y_fit = slopes[fit_mask]
            y_err = slope_errs[fit_mask]

            # Weighted linear least squares (chi-square fit)
            w = 1 / (y_err ** 2)
            W = np.sum(w)
            Wx = np.sum(w * x_fit)
            Wy = np.sum(w * y_fit)
            Wxx = np.sum(w * x_fit * x_fit)
            Wxy = np.sum(w * x_fit * y_fit)
            Delta = W * Wxx - Wx ** 2

            slope_fit = (W * Wxy - Wx * Wy) / Delta
            intercept_fit = (Wxx * Wy - Wx * Wxy) / Delta

            fit_line = lambda x: slope_fit * x + intercept_fit

            x_crit = -intercept_fit / slope_fit
            x_fit_full = np.linspace(x_fit.min(), x_crit, 200)
            ax.plot(x_fit_full, fit_line(x_fit_full), 'r--', label=f"Fit m = {slope_fit:.3f}")
            ax.axvline(x_crit, color='k', linestyle=':', label=f"Critical 1/S = {x_crit:.3f}")
            ax.set_title(f"{dim}D: (1/S)$_\\text{{c}}$ = {x_crit:.2f}")

            ax.legend()

    plt.tight_layout()
    plt.savefig(f"src/probabilityPuyoPuyo/plots/drift/driftVsInvEntropy.png", dpi=300)
    plt.show()


def drift_vs_invEntropy_randomProbs():
    data_dir = "src/probabilityPuyoPuyo/outputs/randomProbabilities2D"
    files = glob.glob(f"{data_dir}/L_*_P_*.tsv")

    drifts = []
    drift_errs = []
    inv_entropies = []

    interval = 100
    tolerance = 0.01

    prob_count = collections.Counter()

    for file in files:
        # Extract probabilities from filename using regex
        match = re.search(r'_P_([0-9\.\-]+)\.tsv$', file)
        if not match:
            continue
        prob_str = match.group(1)
        probs = np.array([float(p) for p in prob_str.split('-')])
        probs /= probs.sum()  # Ensure normalization

        prob_count[len(probs)] += 1  # Count by number of probabilities

        # Calculate entropy
        S = -np.sum(probs * np.log(probs))
        inv_S = 1 / S

        # Load data
        try:
            step, mass = np.loadtxt(file, delimiter="\t", skiprows=1, unpack=True)
        except Exception as e:
            print(f"Could not load {file}: {e}")
            continue

        # Bin into intervals and compute mean normalized mass
        n_bins = len(step) // interval
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

        if len(t_bins) < 2:
            continue

        # Linear fit on the binned means
        slope, intercept = np.polyfit(t_bins, m_means, 1)

        # Estimate error on slope from residuals
        m_pred = slope * t_bins + intercept
        resid = m_means - m_pred
        slope_err = np.std(resid) / np.sqrt(np.sum((t_bins - np.mean(t_bins))**2))

        drifts.append(slope)
        drift_errs.append(slope_err)
        inv_entropies.append(inv_S)

    # Print the number of files for each N
    for N, count in prob_count.items():
        print(f"{count} files have {N} probabilities")


    # Convert to arrays and sort by inv_entropy for plotting
    drifts = np.array(drifts)
    drift_errs = np.array(drift_errs)
    inv_entropies = np.array(inv_entropies)
    sort_idx = np.argsort(inv_entropies)
    drifts = drifts[sort_idx]
    drift_errs = drift_errs[sort_idx]
    inv_entropies = inv_entropies[sort_idx]

    # Plot
    plt.figure(figsize=(8, 6))
    plt.errorbar(inv_entropies, drifts, yerr=drift_errs, fmt='o', capsize=4, label="All points")
    plt.xlabel("1 / Entropy S")
    plt.ylabel("Drift (mass/time)")
    plt.title("Drift vs 1/Entropy (random probabilities)")
    plt.grid()


    plt.legend()
    plt.tight_layout()
    plt.savefig("src/probabilityPuyoPuyo/plots/drift/driftVsInvEntropy_randomProbs.png", dpi=300)
    plt.show()

if __name__ == "__main__":
    # mass_vs_time()
    # drift_vs_inverseN()
    # drift_vs_invEntropy()
    drift_vs_invEntropy_randomProbs()