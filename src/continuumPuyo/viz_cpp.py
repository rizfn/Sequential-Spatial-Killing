import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def stretched_exp(x, A, decay, beta):
    return A * np.exp(-decay * np.abs(x)**beta)

def main():
    # 3. Plotting the resulting TSV files produced by the C++ engine
    # Configuration parameters to match C++ output
    L = 128
    T = 8000
    dt = 0.01
    a = 1
    c_val = 1
    lambda_abs = 10

    print("Plotting results from C++ raw data...")
    param_str = f"L{L}_T{T}_dt{dt}_a{a}_c{c_val}_lam{lambda_abs}"
    out_dir = f"src/continuumPuyo/outputs/avalanchePDE_cpp/{param_str}"

    Ns = [1, 2, 4, 8, 12, 15, 20, 40]

    fig_r, ax_r = plt.subplots(figsize=(7, 5))
    fig_p, ax_p = plt.subplots(figsize=(7, 5))
    fig_s, ax_s = plt.subplots(figsize=(7, 5))

    min_density = float('inf')
    max_density = float('-inf')

    for N in Ns:
        # 1. Load and plot Roughness
        t_vals, r_vals = np.loadtxt(f"{out_dir}/roughness_N{N}.tsv", unpack=True, skiprows=1)
        ax_r.plot(t_vals, r_vals, label=f'N={N}')
        
        # 2. Load and plot Profile
        x_vals, h_vals = np.loadtxt(f"{out_dir}/profile_N{N}.tsv", unpack=True, skiprows=1)
        ax_p.plot(x_vals, h_vals, alpha=0.8, label=f'N={N}')
        
        # 3. Load and plot Slopes
        bins, counts = np.loadtxt(f"{out_dir}/slopes_N{N}.tsv", unpack=True, skiprows=1)
        bin_width = bins[1] - bins[0] if len(bins) > 1 else 1.0
        # Convert raw histogram counts into a valid probability density distribution
        total_counts = np.sum(counts)
        density = counts / (total_counts * bin_width + 1e-12)
        
        valid_density = density[density > 0]
        if len(valid_density) > 0:
            min_density = min(min_density, np.min(valid_density))
            max_density = max(max_density, np.max(valid_density))

        mask = density > 0
        popt = None
        try:
            # Fit density function: A * exp(-decay * |x|^beta)
            popt, _ = curve_fit(stretched_exp, bins[mask], density[mask], 
                                p0=[np.max(density), 1.0, 1.0], maxfev=10000)
            A_fit, decay_fit, beta_fit = popt
            print(f"N={N:2d} | decay = {decay_fit:.4f}, exponent = {beta_fit:.4f}")
            lbl_data = fr'N={N} ($\lambda$={decay_fit:.2f}, $\beta$={beta_fit:.2f})'
        except Exception as e:
            print(f"N={N:2d} | fit failed: {e}")
            lbl_data = f'N={N}'

        p = ax_s.plot(bins, density, label=lbl_data, lw=2)
        if popt is not None:
            ax_s.plot(bins, stretched_exp(bins, *popt), '--', color=p[0].get_color(), alpha=0.7)

    # Reference lines for roughness
    t_ref = t_vals[t_vals > 0]
    mid = len(t_ref) // 2
    if len(r_vals) > mid and r_vals[mid] > 0:
        val = r_vals[mid]
        t_mid = t_ref[mid]
        ax_r.plot(t_ref, val * (t_ref/t_mid)**(1/3), 'k--', label=r'$t^{1/3}$', zorder=10)
        ax_r.plot(t_ref, val * (t_ref/t_mid)**(1/2), 'k:', label=r'$t^{1/2}$', zorder=10)

    plots_dir = f"src/continuumPuyo/plots/cpp_L{L}_T{T}_dt{dt}_a{a}_c{c_val}_lam{lambda_abs}"
    os.makedirs(plots_dir, exist_ok=True)

    # Formalize Roughness Plot
    ax_r.set_xscale('log')
    ax_r.set_yscale('log')
    ax_r.set_xlabel('Time (t)')
    ax_r.set_ylabel('Average Roughness W(t)')
    ax_r.set_title('Roughness vs Time')
    ax_r.legend()
    fig_r.tight_layout()
    fig_r.savefig(os.path.join(plots_dir, 'roughness_cpp.png'))

    # Formalize Profile Plot
    ax_p.set_xlabel('Position (x)')
    ax_p.set_ylabel('Height (h - <h>)')
    ax_p.set_title('Final Height Profile')
    ax_p.legend()
    fig_p.tight_layout()
    fig_p.savefig(os.path.join(plots_dir, 'profile_cpp.png'))

    # Formalize Slopes Plot
    ax_s.set_yscale('log')
    if min_density < float('inf'):
        ax_s.set_ylim(bottom=min_density * 0.5, top=max_density * 2.0)
    ax_s.set_xlabel('Local Slope (grad h)')
    ax_s.set_ylabel('Density')
    ax_s.set_title('Distribution of Local Slopes (Tail 10% steps)')
    ax_s.legend()
    fig_s.tight_layout()
    fig_s.savefig(os.path.join(plots_dir, 'slopes_cpp.png'))

    print(f"Done! Plots saved in {plots_dir}")

if __name__ == "__main__":
    main()