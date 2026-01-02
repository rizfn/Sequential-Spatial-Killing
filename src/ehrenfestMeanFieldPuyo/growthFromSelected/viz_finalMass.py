import os
import re
import glob
import numpy as np
import matplotlib.pyplot as plt

def main():
    # Define the directory containing the output files
    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_dir = os.path.join(script_dir, "outputs", "finalMass")
    
    # Use glob to find all output files
    file_pattern = os.path.join(output_dir, "N*_K*_steps*_init*.txt")
    files = glob.glob(file_pattern)
    
    if not files:
        print(f"No output files found in {output_dir}")
        return
    
    print(f"Found {len(files)} output files")
    
    # Parse filenames and read data
    data = []
    for file_path in files:
        filename = os.path.basename(file_path)
        
        # Extract parameters from filename using regex
        match = re.match(r'N(\d+)_K(\d+)_steps(\d+)_init(\d+)\.txt', filename)
        if match:
            N_colors = int(match.group(1))
            K_selections = int(match.group(2))
            N_steps = int(match.group(3))
            initial_total_balls = int(match.group(4))
            
            # Read final size from file
            with open(file_path, 'r') as f:
                final_size = int(f.read().strip())
            
            data.append((N_colors, K_selections, N_steps, initial_total_balls, final_size))
    
    if not data:
        print("No valid data found in files")
        return
    
    # Extract unique values for N, K, steps, and init
    N_values = sorted(set(d[0] for d in data))
    K_values = sorted(set(d[1] for d in data))
    N_steps = data[0][2]  # Assuming all files have same N_steps
    initial_total_balls = data[0][3]  # Assuming all files have same initial_total_balls
    
    print(f"Parameter space: N={len(N_values)} values, K={len(K_values)} values")
    print(f"N range: [{min(N_values)}, {max(N_values)}]")
    print(f"K range: [{min(K_values)}, {max(K_values)}]")
    
    # Create a 2D array for the heatmap
    final_sizes = np.full((len(K_values), len(N_values)), np.nan)
    
    # Fill the array with data
    for N_colors, K_selections, _, _, final_size in data:
        i = K_values.index(K_selections)
        j = N_values.index(N_colors)
        final_sizes[i, j] = final_size
    
    # Create heatmap
    plt.figure(figsize=(8, 6))
    
    # Use pcolormesh for better control with log scales
    N_array = np.array(N_values)
    K_array = np.array(K_values)
    
    # Create mesh edges for pcolormesh
    N_edges = np.concatenate([[N_array[0] * 0.9], 
                              (N_array[:-1] + N_array[1:]) / 2, 
                              [N_array[-1] * 1.1]])
    K_edges = np.concatenate([[K_array[0] - 0.5], 
                              (K_array[:-1] + K_array[1:]) / 2, 
                              [K_array[-1] + 0.5]])
    
    im = plt.pcolormesh(N_edges, K_edges, final_sizes, 
                        cmap='RdBu_r', vmin=0, vmax=2*initial_total_balls,
                        shading='flat')
    
    plt.colorbar(im, label='Final Urn Size')
    plt.xlabel('N (Number of Colors)', fontsize=12)
    plt.ylabel('K (Selection Size)', fontsize=12)
    plt.title(f'Final Urn Size Heatmap (N_steps={N_steps}, Initial={initial_total_balls})', fontsize=14)
    
    ylim = plt.ylim()

    # Analytical boundary: N = 3K(K-1)/2
    K_analytical = np.linspace(min(K_values), max(K_values), 200)
    N_analytical = 3 * K_analytical * (K_analytical - 1) / 2
    plt.plot(N_analytical, K_analytical, 'k--', linewidth=2, label=r'Analytical: $N = \frac{3}{2}K(K-1)$')
    plt.ylim(ylim)

    plt.legend(loc='upper left')
    plt.xscale('log')
    plt.yscale('log')
    plt.tight_layout()
    
    # Save plot
    plot_dir = os.path.join(script_dir, "plots")
    os.makedirs(plot_dir, exist_ok=True)
    plot_path = os.path.join(plot_dir, f"finalMass_heatmap_steps{N_steps}_init{initial_total_balls}.png")
    plt.savefig(plot_path, dpi=300)
    print(f"Heatmap saved to {plot_path}")

def main_probabilities():
    # Define the directory containing the output files
    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_dir = os.path.join(script_dir, "outputs", "finalMassProbabilities")
    
    # Use glob to find all output files
    file_pattern = os.path.join(output_dir, "N*_K*_steps*_init*_rate*.txt")
    files = glob.glob(file_pattern)
    
    if not files:
        print(f"No output files found in {output_dir}")
        return
    
    print(f"Found {len(files)} output files")
    
    # Parse filenames and read data
    data = []
    for file_path in files:
        filename = os.path.basename(file_path)
        
        # Extract parameters from filename using regex (including rate parameter)
        match = re.match(r'N(\d+)_K(\d+)_steps(\d+)_init(\d+)_rate([\d.]+)\.txt', filename)
        if match:
            N_colors = int(match.group(1))
            K_selections = int(match.group(2))
            N_steps = int(match.group(3))
            initial_total_balls = int(match.group(4))
            exponential_rate = float(match.group(5))
            
            # Read final size from file
            with open(file_path, 'r') as f:
                final_size = int(f.read().strip())
            
            data.append((N_colors, K_selections, N_steps, initial_total_balls, exponential_rate, final_size))
    
    if not data:
        print("No valid data found in files")
        return
    
    # Group data by rate parameter
    rate_values = sorted(set(d[4] for d in data))
    
    for rate in rate_values:
        # Filter data for this rate
        rate_data = [(d[0], d[1], d[2], d[3], d[5]) for d in data if d[4] == rate]
        
        # Extract unique values for N, K, steps, and init
        N_values = sorted(set(d[0] for d in rate_data))
        K_values = sorted(set(d[1] for d in rate_data))
        N_steps = rate_data[0][2]  # Assuming all files have same N_steps
        initial_total_balls = rate_data[0][3]  # Assuming all files have same initial_total_balls
        
        print(f"\nRate {rate}: Parameter space: N={len(N_values)} values, K={len(K_values)} values")
        print(f"N range: [{min(N_values)}, {max(N_values)}]")
        print(f"K range: [{min(K_values)}, {max(K_values)}]")
        
        # Create a 2D array for the heatmap
        final_sizes = np.full((len(K_values), len(N_values)), np.nan)
        
        # Fill the array with data
        for N_colors, K_selections, _, _, final_size in rate_data:
            i = K_values.index(K_selections)
            j = N_values.index(N_colors)
            final_sizes[i, j] = final_size
        
        # Create heatmap
        plt.figure(figsize=(8, 6))
        
        # Use pcolormesh for better control with log scales
        N_array = np.array(N_values)
        K_array = np.array(K_values)
        
        # Create mesh edges for pcolormesh
        N_edges = np.concatenate([[N_array[0] * 0.9], 
                                  (N_array[:-1] + N_array[1:]) / 2, 
                                  [N_array[-1] * 1.1]])
        K_edges = np.concatenate([[K_array[0] - 0.5], 
                                  (K_array[:-1] + K_array[1:]) / 2, 
                                  [K_array[-1] + 0.5]])
        
        im = plt.pcolormesh(N_edges, K_edges, final_sizes, 
                            cmap='RdBu_r', vmin=0, vmax=2*initial_total_balls,
                            shading='flat')
        
        plt.colorbar(im, label='Final Urn Size')
        plt.xlabel('N (Number of Colors)', fontsize=12)
        plt.ylabel('K (Selection Size)', fontsize=12)
        plt.title(f'Final Urn Size (Exponential Rate={rate}, Steps={N_steps}, Init={initial_total_balls})', fontsize=13)
        
        ylim = plt.ylim()
        
        # Analytical boundary: N = 3K(K-1)/2
        K_analytical = np.linspace(min(K_values), max(K_values), 200)
        N_analytical = 3 * K_analytical * (K_analytical - 1) / 2
        plt.plot(N_analytical, K_analytical, 'k--', linewidth=2, label=r'Analytical: $N = \frac{3}{2}K(K-1)$')
        plt.ylim(ylim)
        
        plt.legend(loc='upper left')
        plt.xscale('log')
        plt.yscale('log')
        plt.tight_layout()
        
        # Save plot
        plot_dir = os.path.join(script_dir, "plots")
        os.makedirs(plot_dir, exist_ok=True)
        plot_path = os.path.join(plot_dir, f"finalMassProbabilities_rate{rate}_steps{N_steps}_init{initial_total_balls}.png")
        plt.savefig(plot_path, dpi=300)
        print(f"Heatmap saved to {plot_path}")
        plt.close()

    
if __name__ == "__main__":
    # main()
    main_probabilities()