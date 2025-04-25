import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

def slopeDistVsTime3D():
    L = 128
    N = 6  # Choose a single N value for visualization
    N_steps = 1024

    # Load slope data
    slope_data = np.loadtxt(
        f"src/puyopuyo/cpp/outputs/slopeDistribution2D/L_{L}_N_{N}_steps_{N_steps}.tsv",
        delimiter="\t",
        skiprows=1,
        dtype=str
    )

    # Extract steps and slopes
    steps = slope_data[:, 0].astype(int)
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
    plt.savefig(f"src/puyopuyo/cpp/plots/slopeDistribution/2D_L_{L}_N_{N}.png", dpi=300)
    plt.show()


def timeAveragedSlopes():
    L = 128
    N_list = np.arange(2, 13)
    N_steps = 1024

    fig, ax = plt.subplots(figsize=(10, 6))

    # Use a colormap (e.g., rainbow) to assign colors based on N_species
    colormap = plt.get_cmap("rainbow", len(N_list))

    for i, N in enumerate(N_list):
        # Load slope data
        slope_data = np.loadtxt(
            f"src/puyopuyo/cpp/outputs/slopeDistribution2D/L_{L}_N_{N}_steps_{N_steps}.tsv",
            delimiter="\t",
            skiprows=1,
            dtype=str
        )

        # Extract steps and slopes
        steps = slope_data[:, 0].astype(int)
        slopes = slope_data[:, 2]  # The third column contains the comma-separated slopes

        # Flatten the slope data into a single array for distribution
        all_slopes = []
        for slope_row in slopes:
            all_slopes.extend(map(int, slope_row.split(",")))

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

    plt.savefig(f"src/puyopuyo/cpp/plots/slopeDistribution/2D_timeAverage_L_{L}.png", dpi=300)
    plt.show()




def plotInterface():
    L = 128
    N_list = [4, 6, 10]  # Number of species
    time_snapshots = [128, 256, 512, 1024]  # Time steps to plot
    N_steps = 1024

    # Create a grid of subplots
    fig = plt.figure(figsize=(16, 12))
    gs = GridSpec(len(N_list), len(time_snapshots), figure=fig)

    for row, N in enumerate(N_list):
        # Load slope data for the current number of species
        slope_data = np.loadtxt(
            f"src/puyopuyo/cpp/outputs/slopeDistribution2D/L_{L}_N_{N}_steps_{N_steps}.tsv",
            delimiter="\t",
            skiprows=1,
            dtype=str
        )

        # Extract steps, first column heights, and slopes
        steps = slope_data[:, 0].astype(int)
        first_column_heights = slope_data[:, 1].astype(int)
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
            ax = fig.add_subplot(gs[row, col])
            x = range(L)
            ax.fill_between(x, np.min(heights[:L]), heights[:L], step='mid', color="blue", alpha=0.3)

            # Set labels and titles
            if row == len(N_list) - 1:
                ax.set_xlabel("Column Index")
            if col == 0:
                ax.set_ylabel(f"N={N}")
            if row == 0:
                ax.set_title(f"Time={time}")

            ax.grid()
            ax.legend()

    # Adjust layout and save the plot
    plt.tight_layout()
    plt.savefig(f"src/puyopuyo/cpp/plots/slopeDistribution/interfaceGrid_L_{L}.png", dpi=300)
    plt.show()


if __name__ == "__main__":
    # slopeDistVsTime3D()
    # timeAveragedSlopes()
    plotInterface()