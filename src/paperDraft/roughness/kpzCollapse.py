import os
import glob
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def get_roughness(slope_str, h0):
    slopes = np.fromstring(slope_str, sep=',')
    # Reconstruct heights: h_{i+1} = h_i + slope_i
    heights = np.concatenate(([0], np.cumsum(slopes[:-1]))) + h0
    return np.std(heights)

def load_collapse_data(steps, fixed_N):
    data_frames = []
    pattern = re.compile(f"L_(\\d+)_N_{fixed_N}_steps_{steps}_sim_(\\d+)\\.tsv")
    
    # Locate output directory
    output_dir = "src/paperDraft/roughness/outputs/slopeDist"
    files = glob.glob(f"{output_dir}/L_*_N_{fixed_N}_steps_{steps}_sim_*.tsv")
    
    if not files:
        raise FileNotFoundError(f"No files found in {output_dir} for N={fixed_N}, steps={steps}")

    for f in files:
        match = pattern.search(f)
        if not match:
            continue
            
        L, simno = int(match.group(1)), int(match.group(2))
        
        df = pd.read_csv(f, sep='\t')
        
        # Ensure we have the expected columns
        if 'slope_distribution' in df.columns and 'first_col_height' in df.columns:
            df['roughness'] = df.apply(lambda r: get_roughness(r['slope_distribution'], r['first_col_height']), axis=1)
        elif 'roughness' not in df.columns:
            print(f"Warning: Unexpected columns in {f}")
            continue
            
        df['L'] = L
        df['sim_no'] = simno
        
        data_frames.append(df[['L', 'sim_no', 'step', 'roughness']])
    
    return pd.concat(data_frames, ignore_index=True)

def plot_kpz_collapse(df, steps, fixed_N, output_path):
    plt.figure(figsize=(10, 6))
    
    # Aggregate mean across sim_no for each L and step
    grouped = df.groupby(['L', 'step'])['roughness'].mean().reset_index()
    
    L_values = sorted(grouped['L'].unique())
    colormap = plt.get_cmap("viridis", len(L_values))
    
    # According to 1D KPZ scaling: W(L, t) ~ L^alpha * f(t / L^z)
    alpha = 0.5
    z = 1.5
    
    for i, L in enumerate(L_values):
        curr = grouped[grouped['L'] == L].copy()
        
        # Avoid division by zero at t=0
        mask = curr['step'] > 0
        curr = curr[mask]
        
        t_scaled = curr['step'] / (L ** z)
        w_scaled = curr['roughness'] / (L ** alpha)
        
        plt.plot(t_scaled, w_scaled, label=f"L={L}", color=colormap(i))

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'Scaled Time $t / L^{3/2}$')
    plt.ylabel(r'Scaled Roughness $W / L^{1/2}$')
    plt.title(f'KPZ Data Collapse (N={fixed_N})')
    plt.grid(True, which="both", ls="--", alpha=0.4)
    plt.legend()
    
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()

def main():
    steps = 8000
    fixed_N = 10
    
    print(f"Loading data for N={fixed_N}, steps={steps}...")
    df = load_collapse_data(steps, fixed_N)
    
    plot_dir = "src/paperDraft/roughness/plots/kpzCollapse"
    os.makedirs(plot_dir, exist_ok=True)
    
    plot_path = f"{plot_dir}/kpz_collapse_N_{fixed_N}_steps_{steps}.png"
    plot_kpz_collapse(df, steps, fixed_N, plot_path)
    print(f"Saved KPZ collapse plot: {plot_path}")

if __name__ == "__main__":
    main()
