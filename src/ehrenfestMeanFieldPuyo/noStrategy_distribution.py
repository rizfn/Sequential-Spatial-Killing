import os
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool, cpu_count

def simulate_urn(N_colors, K_selections, N_steps, initial_total_balls):
    """
    Simulate the urn evolution and return the final urn state.
    """
    # Initialize urn with equal number of balls of each color
    initial_balls_per_color = initial_total_balls // N_colors
    urn = np.full(N_colors, initial_balls_per_color, dtype=int)

    for step in range(N_steps):
        total_balls = np.sum(urn)
        
        # Can't draw K balls if we don't have enough
        if total_balls < K_selections:
            break
        
        # Create a list of all balls by expanding the urn counts
        all_balls = []
        for color in range(N_colors):
            all_balls.extend([color] * urn[color])
        
        # Randomly draw K balls without replacement
        drawn_balls = np.random.choice(all_balls, size=K_selections, replace=False)
        
        # Count occurrences of each color in the drawn balls
        unique_colors, counts = np.unique(drawn_balls, return_counts=True)
        
        # Check if all drawn balls are unique (no duplicates)
        if len(unique_colors) == K_selections:
            # All K balls are different colors - put them back and add one extra
            random_color = np.random.choice(N_colors)
            urn[random_color] += 1
        else:
            # There are duplicates - remove all balls of colors that appeared more than once
            for color, count in zip(unique_colors, counts):
                if count > 1:
                    urn[color] -= count
    
    return urn


def main():
    # Single simulation parameters
    N_colors = 300
    K_selections = 18
    N_steps = 10000
    initial_total_balls = 1000
    
    print(f"Running simulation with N={N_colors}, K={K_selections}...")
    
    # Run single simulation
    final_urn = simulate_urn(N_colors, K_selections, N_steps, initial_total_balls)
    
    # Filter out colors with zero balls
    non_zero_counts = final_urn[final_urn > 0]
    
    print(f"Final total balls: {np.sum(final_urn)}")
    print(f"Number of surviving colors: {len(non_zero_counts)} out of {N_colors}")
    print(f"Max balls of one color: {np.max(final_urn)}")
    print(f"Mean balls per surviving color: {np.mean(non_zero_counts):.2f}")
    
    # Create three subplots
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(20, 6))
    
    # Plot 1: Color index vs number of balls
    colors = np.arange(N_colors)
    ax1.bar(colors, final_urn, width=1.0, edgecolor='none', alpha=0.8)
    ax1.set_xlabel('Color Index', fontsize=12)
    ax1.set_ylabel('Number of Balls', fontsize=12)
    ax1.set_title('Color Distribution by Index', fontsize=13)
    ax1.grid(True, alpha=0.3, axis='y')
    
    # Plot 2: Rank-ordered distribution
    sorted_counts = np.sort(final_urn)[::-1]  # Sort in descending order
    ranks = np.arange(1, N_colors + 1)
    ax2.bar(ranks, sorted_counts, width=1.0, edgecolor='none', alpha=0.8)
    ax2.set_xlabel('Rank', fontsize=12)
    ax2.set_ylabel('Number of Balls', fontsize=12)
    ax2.set_title('Rank-Ordered Color Distribution', fontsize=13)
    ax2.grid(True, alpha=0.3, axis='y')
    
    # Plot 3: Log-log plot (Zipfian)
    # Filter out zeros for log-log plot
    non_zero_sorted = sorted_counts[sorted_counts > 0]
    non_zero_ranks = np.arange(1, len(non_zero_sorted) + 1)
    ax3.loglog(non_zero_ranks, non_zero_sorted, 'o', alpha=0.6, markersize=4)
    ax3.set_xlabel('Rank (log scale)', fontsize=12)
    ax3.set_ylabel('Number of Balls (log scale)', fontsize=12)
    ax3.set_title('Zipfian Plot (log-log)', fontsize=13)
    ax3.grid(True, alpha=0.3, which='both')
    
    fig.suptitle(f'Color Distribution (N={N_colors}, K={K_selections}, Steps={N_steps})', fontsize=14, y=1.02)
    plt.tight_layout()
    
    os.makedirs("src/ehrenfestMeanFieldPuyo/plots/finalDistribution", exist_ok=True)
    plt.savefig(f"src/ehrenfestMeanFieldPuyo/plots/finalDistribution/nostrategy_N{N_colors}_K{K_selections}.png", dpi=300)
    print(f"Distribution saved to src/ehrenfestMeanFieldPuyo/plots/finalDistribution/nostrategy_N{N_colors}_K{K_selections}.png")


if __name__ == "__main__":
    main()