import os
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool, cpu_count

def simulate_urn(N_colors, K_selections, N_steps, initial_total_balls):
    """
    Simulate the urn evolution and return the final total number of balls.
    """
    # Initialize urn with equal number of balls of each color
    initial_balls_per_color = initial_total_balls // N_colors
    urn = np.full(N_colors, initial_balls_per_color, dtype=int)
    
    # Track the most abundant color(s) from the last elimination
    last_eliminated_color = None

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
            # Add the color that was most abundant in the last elimination
            if last_eliminated_color is not None:
                random_color = last_eliminated_color
            else:
                # First time adding, choose randomly
                random_color = np.random.randint(0, N_colors)
            urn[random_color] += 1
        else:
            # There are duplicates - remove all balls of colors that appeared more than once
            # Track which color had the most duplicates
            max_duplicate_count = 0
            colors_with_max_duplicates = []
            
            for color, count in zip(unique_colors, counts):
                if count > 1:
                    urn[color] -= count
                    if count > max_duplicate_count:
                        max_duplicate_count = count
                        colors_with_max_duplicates = [color]
                    elif count == max_duplicate_count:
                        colors_with_max_duplicates.append(color)
            
            # Remember the most abundant eliminated color (random if tie)
            if colors_with_max_duplicates:
                last_eliminated_color = np.random.choice(colors_with_max_duplicates)
    
    return np.sum(urn)


def worker(params):
    """
    Worker function for multiprocessing.
    """
    N_colors, K_selections, N_steps, initial_total_balls = params
    final_size = simulate_urn(N_colors, K_selections, N_steps, initial_total_balls)
    return (N_colors, K_selections, final_size)


def main():
    # Parameter ranges
    N_colors_range = np.arange(2, 500, 4)  # 
    K_selections_range = np.arange(2, 20, 1)  # 2 to 8
    N_steps = 10000
    initial_total_balls = 1000
    
    # Create parameter combinations
    params_list = []
    for N_colors in N_colors_range:
        for K_selections in K_selections_range:
            if K_selections <= initial_total_balls:  # Only simulate if K <= initial total balls
                params_list.append((N_colors, K_selections, N_steps, initial_total_balls))
    
    print(f"Running {len(params_list)} simulations using {cpu_count()} cores...")
    
    # Run simulations in parallel
    with Pool() as pool:
        results = pool.map(worker, params_list)
    
    # Organize results into a 2D array for heatmap
    final_sizes = np.zeros((len(K_selections_range), len(N_colors_range)))
    
    for N_colors, K_selections, final_size in results:
        i = np.where(K_selections_range == K_selections)[0][0]
        j = np.where(N_colors_range == N_colors)[0][0]
        final_sizes[i, j] = final_size
    
    # Create heatmap
    plt.figure(figsize=(12, 10))
    im = plt.imshow(final_sizes, aspect='auto', origin='lower', 
                    extent=[N_colors_range[0], N_colors_range[-1], 
                           K_selections_range[0], K_selections_range[-1]],
                    cmap='RdBu_r', vmin=0, vmax=2*initial_total_balls)
    plt.colorbar(im, label='Final Urn Size')
    plt.title(f'Final Urn Size Heatmap (N_steps={N_steps}, Initial={initial_total_balls})', fontsize=14)
    
    ylim = plt.ylim()
    # Plot theoretical transition: K ~ 0.9*sqrt(N)
    plt.plot(N_colors_range, 0.90*np.sqrt(np.array(N_colors_range)), color='black', linestyle='--', linewidth=2, label=r'$K = 0.90\sqrt{N}$')
    plt.ylim(ylim)
    plt.tight_layout()
    
    os.makedirs("src/ehrenfestMeanFieldPuyo/plots/massVsTime", exist_ok=True)
    plt.savefig(f"src/ehrenfestMeanFieldPuyo/plots/massVsTime/prevEliminatedStrategy_steps_{N_steps}_init_{initial_total_balls}.png", dpi=300)
    print(f"Heatmap saved to src/ehrenfestMeanFieldPuyo/plots/massVsTime/prevEliminatedStrategy_steps_{N_steps}_init_{initial_total_balls}.png")


if __name__ == "__main__":
    main()