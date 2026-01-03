import os
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import glob

def load_simulation_data(file_pattern):
    """
    Load all simulation files matching the pattern.
    Returns lists of distributions and reproductivities.
    """
    files = glob.glob(file_pattern)
    
    if not files:
        raise FileNotFoundError(f"No files found matching pattern: {file_pattern}")
    
    distributions = []
    reproductivities = []
    
    for file_path in sorted(files):
        data = np.loadtxt(file_path)
        distributions.append(data[:, 0].astype(int))  # First column: counts
        reproductivities.append(data[:, 1])  # Second column: reproductivity
    
    return distributions, reproductivities, len(files)


def main():
    # Define parameters matching the simulation
    N_colors = 4000
    K_selections = 30
    N_steps = 200000
    initial_total_balls = 10000
    exponential_rate = -1
    
    # Construct file pattern
    script_dir = Path(__file__).parent
    file_pattern = str(script_dir / f"outputs/finalDistributionProbabilities/N{N_colors}_K{K_selections}_steps{N_steps}_init{initial_total_balls}_rate{exponential_rate:.2f}_sim*.txt")
    
    print(f"Loading simulation data from: {file_pattern}")
    
    # Load all simulations
    distributions, reproductivities, num_sims = load_simulation_data(file_pattern)
    
    print(f"Loaded {num_sims} simulations")
    
    # Convert to numpy arrays for easier manipulation
    distributions = np.array(distributions)  # Shape: (num_sims, N_colors)
    reproductivities = np.array(reproductivities)  # Shape: (num_sims, N_colors)
    
    # Calculate statistics across simulations
    mean_distribution = np.mean(distributions, axis=0)
    std_distribution = np.std(distributions, axis=0)
    
    # Use reproductivity from first simulation (they're all the same structure)
    reprod = reproductivities[0]
    
    # Statistics
    total_balls = np.sum(mean_distribution)
    surviving_colors = np.sum(mean_distribution > 0)
    
    print(f"Average total balls: {total_balls:.2f}")
    print(f"Average surviving colors: {surviving_colors:.2f} out of {N_colors}")
    print(f"Max average balls: {np.max(mean_distribution):.2f}")
    
    # Create figure with multiple subplots
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # Plot 1: Mean distribution by color index
    ax1 = axes[0, 0]
    colors_idx = np.arange(N_colors)
    ax1.bar(colors_idx, mean_distribution, width=1.0, edgecolor='none', alpha=0.8)
    ax1.set_xlabel('Color Index (sorted by reproductivity)', fontsize=11)
    ax1.set_ylabel('Mean Number of Balls', fontsize=11)
    ax1.set_title(f'Mean Distribution by Color Index ({num_sims} simulations)', fontsize=12)
    ax1.grid(True, alpha=0.3, axis='y')
    
    # Plot 2: Distribution vs reproductivity
    ax2 = axes[0, 1]
    # Only plot colors with mean > 0
    mask = mean_distribution > 0
    ax2.scatter(reprod[mask], mean_distribution[mask], alpha=0.6, s=20)
    ax2.errorbar(reprod[mask], mean_distribution[mask], yerr=std_distribution[mask], 
                 fmt='none', alpha=0.3, capsize=2)
    ax2.set_xlabel('Reproductivity', fontsize=11)
    ax2.set_ylabel('Mean Number of Balls', fontsize=11)
    ax2.set_title('Distribution vs Reproductivity', fontsize=12)
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Rank-ordered distribution
    ax3 = axes[1, 0]
    sorted_mean = np.sort(mean_distribution)[::-1]
    ranks = np.arange(1, N_colors + 1)
    ax3.bar(ranks, sorted_mean, width=1.0, edgecolor='none', alpha=0.8)
    ax3.set_xlabel('Rank', fontsize=11)
    ax3.set_ylabel('Mean Number of Balls', fontsize=11)
    ax3.set_title('Rank-Ordered Distribution', fontsize=12)
    ax3.grid(True, alpha=0.3, axis='y')
    ax3.set_xlim(0, N_colors)
    
    # Plot 4: Log-log plot
    ax4 = axes[1, 1]
    non_zero_sorted = sorted_mean[sorted_mean > 0]
    non_zero_ranks = np.arange(1, len(non_zero_sorted) + 1)
    ax4.loglog(non_zero_ranks, non_zero_sorted, 'o', alpha=0.6, markersize=3)
    ax4.set_xlabel('Rank (log scale)', fontsize=11)
    ax4.set_ylabel('Mean Number of Balls (log scale)', fontsize=11)
    ax4.set_title('Zipfian Plot (log-log)', fontsize=12)
    ax4.grid(True, alpha=0.3, which='both')
    
    fig.suptitle(f'Distribution with Reproductivity (N={N_colors}, K={K_selections}, rate={exponential_rate})', 
                 fontsize=14, y=0.995)
    plt.tight_layout()
    
    # Save figure
    output_dir = script_dir / "plots/finalDistributionProbabilities"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / f"distribution_N{N_colors}_K{K_selections}_rate{exponential_rate:.2f}.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nPlot saved to: {output_file}")
    

if __name__ == "__main__":
    main()
