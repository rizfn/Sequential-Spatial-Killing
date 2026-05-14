import os
import glob
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def get_roughness(slope_str, h0):
    slopes = np.fromstring(slope_str, sep=',')
    # Reconstruct heights: h_{i+1} = h_i + slope_i
    # We only need up to the second to last slope to get all L heights
    heights = np.concatenate(([0], np.cumsum(slopes[:-1]))) + h0
    return np.std(heights)

def load_data(steps):
    data_frames = []
    pattern = re.compile(r"L_(\d+)_N_(\d+)_steps_\d+_sim_(\d+)\.tsv")
    
    # Locate output directory
    output_dir = "src/paperDraft/roughness/outputs/slopeDist"
    files = glob.glob(f"{output_dir}/L_*_N_*_steps_{steps}_sim_*.tsv")
    
    if not files:
        raise FileNotFoundError(f"No files found in {output_dir} for steps={steps}")

    for f in files:
        match = pattern.search(f)
        if not match:
            continue
            
        L, N, simno = int(match.group(1)), int(match.group(2)), int(match.group(3))
        
        df = pd.read_csv(f, sep='\t')
        
        # Ensure we have the expected columns
        if 'slope_distribution' in df.columns and 'first_col_height' in df.columns:
            df['roughness'] = df.apply(lambda r: get_roughness(r['slope_distribution'], r['first_col_height']), axis=1)
        elif 'roughness' not in df.columns:
            print(f"Warning: Unexpected columns in {f}")
            continue
            
        df['L'] = L
        df['N'] = N
        df['sim_no'] = simno
        
        data_frames.append(df[['L', 'N', 'sim_no', 'step', 'roughness']])

    return pd.concat(data_frames, ignore_index=True)

def plot_roughness(df, fixed_col, fixed_val, vary_col, output_path, exponent):
    plt.figure(figsize=(10, 6))
    
    subset = df[df[fixed_col] == fixed_val]
    # Aggregate mean and standard error of the mean (sem) across sim_no
    grouped = subset.groupby([vary_col, 'step'])['roughness'].agg(['mean', 'sem']).reset_index()
    
    colormap = plt.get_cmap("viridis", len(grouped[vary_col].unique()))
    
    for i, val in enumerate(sorted(grouped[vary_col].unique())):
        curr = grouped[grouped[vary_col] == val]
        plt.plot(curr['step'], curr['mean'], label=f"{vary_col}={val}", color=colormap(i))
        plt.fill_between(curr['step'], 
                         curr['mean'] - curr['sem'], 
                         curr['mean'] + curr['sem'], 
                         color=colormap(i), alpha=0.3)

    x = np.linspace(1, df['step'].max(), 100)
    plt.plot(x, np.power(x, exponent), 'k--', label=f'KPZ Scaling ($t^{{ {exponent:.2f} }}$)')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Time')
    plt.ylabel('Roughness')
    plt.title(f'Roughness vs Time ({fixed_col}={fixed_val})')
    plt.grid(True, which="both", ls="--", alpha=0.4)
    plt.legend()
    
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()

def main():
    steps = 8000
    df = load_data(steps)
    
    # Get median values
    sorted_L = sorted(df['L'].unique())
    sorted_N = sorted(df['N'].unique())
    # median_L = sorted_L[len(sorted_L) // 2]
    # median_N = sorted_N[len(sorted_N) // 2]
    median_L = 128
    median_N = 10

    print(f"Discovered L values: {sorted_L}")
    print(f"Discovered N values: {sorted_N}")
    print(f"Using Median L = {median_L}")
    print(f"Using Median N = {median_N}")
    
    plot_dir = "src/paperDraft/roughness/plots/roughnessVsTime"
    os.makedirs(plot_dir, exist_ok=True)
    
    # plot_1 = f"{plot_dir}/fix_L_{median_L}_vary_N_steps_{steps}.png"
    # plot_roughness(df, 'L', median_L, 'N', plot_1, 1/3)
    # print(f"Saved: {plot_1}")
    
    plot_2 = f"{plot_dir}/fix_N_{median_N}_vary_L_steps_{steps}.png"
    plot_roughness(df, 'N', median_N, 'L', plot_2, 1/3)
    print(f"Saved: {plot_2}")

if __name__ == "__main__":
    main()
