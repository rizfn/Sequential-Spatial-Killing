import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

def slopeDistVsTime3D():
    L = 128
    N = 10  # Choose a single N value for visualization
    N_steps = 8000

    # Load slope data
    files = glob.glob(f"src/paperDraft/roughness/outputs/slopeDist/L_{L}_N_{N}_steps_{N_steps}_sim_*.tsv")
    if not files:
        print(f"Skipping slopeDistVsTime3D, no files found for N={N}.")
        return

    slope_data = np.loadtxt(
        files[0], # Just use the first sim for the 3D plot
        delimiter="\t",
        skiprows=1,
        dtype=str
    )

    # Extract steps and slopes
    steps = slope_data[:, 0].astype(float).astype(int)
    slopes = slope_data[:, 2]  # The third column contains the comma-separated slopes

    # Prepare data for 3D plotting
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Define bins for the histogram
    bins = np.arange(-L, L + 1)

    # Use a colormap to assign colors based on the step index
    colormap = plt.get_cmap("viridis", len(steps) // 10)

    for i in range(0, len(steps), 10):
        # Combine slopes for the current group of 10 steps
        combined_slopes = []
        for j in range(i, min(i + 10, len(steps))):
            slope_values = list(map(int, slopes[j].split(",")))
            combined_slopes.extend(slope_values)

        # Calculate the averaged histogram
        hist, bin_edges = np.histogram(combined_slopes, bins=bins, density=True)

        # Use the bin centers for plotting
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

        # Plot the averaged histogram as a line in 3D space
        x = [steps[i]] * len(bin_centers)  # Use the first step of the group as the x-coordinate
        y = bin_centers                    # Slope values (bins)
        z = hist                           # Histogram density

        ax.plot(x, y, z, color=colormap(i // 10), alpha=0.8)

    # Set labels and title
    ax.set_xlabel("Time (Steps)")
    ax.set_ylabel("Slope")
    ax.set_zlabel("Density")
    ax.set_title(f"Averaged Slope Distribution Over Time (L={L}, N={N})")

    # Adjust the view angle to look head-on at the distribution
    ax.view_init(elev=0, azim=0)

    # Save and show the plot
    import os
    os.makedirs("src/paperDraft/roughness/plots/slopeDistribution", exist_ok=True)
    plt.savefig(f"src/paperDraft/roughness/plots/slopeDistribution/2D_L_{L}_N_{N}.png", dpi=300)


import glob

def timeAveragedSlopes():
    L = 128
    N_list = [10]
    N_steps = 8000

    fig, ax = plt.subplots(figsize=(10, 6))

    # Use a colormap (e.g., rainbow) to assign colors based on N_species
    colormap = plt.get_cmap("rainbow", len(N_list))

    for i, N in enumerate(N_list):
        all_slopes = []
        files = glob.glob(f"src/paperDraft/roughness/outputs/slopeDist/L_{L}_N_{N}_steps_{N_steps}_sim_*.tsv")

        for file in files:
            # Load slope data
            slope_data = np.loadtxt(
                file,
                delimiter="\t",
                skiprows=1,
                dtype=str
            )

            if len(slope_data) == 0:
                continue

            slopes = slope_data[:, 2]

            for slope_row in slopes:
                all_slopes.extend(map(int, slope_row.split(",")))

        if not all_slopes:
            print(f"No slope data found for N={N}")
            continue

        # Plot the histogram of slopes
        color = colormap(i / len(N_list))  # Normalize the index to [0, 1] for the colormap
        ax.hist(
            all_slopes,
            bins=np.arange(-L, L + 1),
            alpha=0.5,
            label=f"N={N}",
            color=color,
            density=True,
            histtype="step"
        )

    ax.set_xlabel("Slope")
    ax.set_ylabel("Probability Density")
    ax.set_yscale("log")
    ax.set_title("Slope Distribution")
    ax.grid()
    ax.legend(ncols=2)

    import os
    os.makedirs(f"src/paperDraft/roughness/plots/slopeDistribution", exist_ok=True)
    plt.savefig(f"src/paperDraft/roughness/plots/slopeDistribution/2D_timeAverage_L_{L}.png", dpi=300)




def plotInterface():
    L = 128
    N_list = [10]  # Number of species
    time_snapshots = [128, 256, 512, 1024, 2048, 4096, 8192]  # Time steps to plot
    N_steps = 8000

    # Create a grid of subplots
    fig = plt.figure(figsize=(16, 12))
    gs = GridSpec(len(N_list), len(time_snapshots), figure=fig)

    for i, N in enumerate(N_list):
        files = glob.glob(f"src/paperDraft/roughness/outputs/slopeDist/L_{L}_N_{N}_steps_{N_steps}_sim_*.tsv")
        if not files:
            print(f"Skipping plotInterface N={N}, no files found.")
            continue
        
        # Load slope data for the current number of species from the first available file
        slope_data = np.loadtxt(
            files[0],
            delimiter="\t",
            skiprows=1,
            dtype=str
        )

        # Extract steps, first column heights, and slopes
        steps = slope_data[:, 0].astype(float).astype(int)
        first_column_heights = slope_data[:, 1].astype(float).astype(int)
        slopes = slope_data[:, 2]  # The third column contains the comma-separated slopes

        for col, time in enumerate(time_snapshots):
            # Find the index of the given time step
            if time not in steps:
                continue
            time_index = np.where(steps == time)[0][0]

            # Extract the height of the first column and the slopes
            first_height = first_column_heights[time_index]
            slope_values = list(map(int, slopes[time_index].split(",")))

            # Calculate the absolute heights
            heights = [first_height]
            for slope in slope_values:
                heights.append(heights[-1] + slope)

            # Plot the interface as a step plot
            ax = fig.add_subplot(gs[i, col])
            x = range(L)
            ax.fill_between(x, np.min(heights[:L]), heights[:L], step='mid', color="blue", alpha=0.3)

            # Set labels and titles
            if i == len(N_list) - 1:
                ax.set_xlabel("Column Index")
            if col == 0:
                ax.set_ylabel(f"N={N}")
            if i == 0:
                ax.set_title(f"Time={time}")

            ax.grid()
            ax.legend()

    # Adjust layout and save the plot
    plt.tight_layout()
    import os
    os.makedirs("src/paperDraft/roughness/plots/slopeDistribution", exist_ok=True)
    plt.savefig(f"src/paperDraft/roughness/plots/slopeDistribution/interfaceGrid_L_{L}.png", dpi=300)


if __name__ == "__main__":
    # slopeDistVsTime3D()
    # timeAveragedSlopes()
    plotInterface()