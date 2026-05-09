import matplotlib.pyplot as plt
import numpy as np
from glob import glob
import os
import re

def main():
    L = 64
    N = 128
    STEPS = 1000

    files = glob(f"src/langmuirRandom/outputs/finalMass/L_{L}_N_{N}/rho_*_steps_{STEPS}_sim_*.txt")

    # Extract rho, simno, and final mass from each file
    densities = []
    final_masses = []
    
    rho_pattern = r"rho_([\d.]+)"
    simno_pattern = r"sim_(\d+)"
    
    for filepath in sorted(files):
        # Extract rho from filename
        rho_match = re.search(rho_pattern, filepath)
        if not rho_match:
            continue
        rho = float(rho_match.group(1))
        
        # Extract simno from filename
        simno_match = re.search(simno_pattern, filepath)
        if not simno_match:
            continue
        simno = int(simno_match.group(1))
        
        # Read final mass from file
        try:
            with open(filepath, 'r') as f:
                final_mass = int(f.read().strip())
        except (IOError, ValueError):
            continue
        
        densities.append(rho)
        final_masses.append(final_mass)

    if len(final_masses) == 0:
        print("No valid data found.")
        return
    
    # Create scatterplot
    plt.figure(figsize=(10, 6))
    plt.scatter(densities, np.array(final_masses) / (STEPS * L), alpha=0.6, s=50)
    plt.xlabel(r'Density ($\rho$)', fontsize=12)
    plt.ylabel('Mass added per step per site', fontsize=12)
    plt.title(f'Drift Velocity vs Density (L={L}, N={N}, STEPS={STEPS})', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    output_path = f'src/langmuirRandom/plots/finalMass_vs_density/L{L}_N{N}_steps{STEPS}.png'
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, dpi=300)

if __name__ == "__main__":
    main()

