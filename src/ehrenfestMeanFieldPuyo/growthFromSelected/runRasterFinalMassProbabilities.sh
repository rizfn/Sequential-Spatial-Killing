#!/bin/bash

# Define constants
N_STEPS=20000
INITIAL_TOTAL_BALLS=10000
EXPONENTIAL_RATE=0.01

# Define parameter ranges
N_COLORS_START=2
N_COLORS_END=4000
N_COLORS_BINS=100 

K_SELECTIONS_START=2
K_SELECTIONS_END=50
K_SELECTIONS_BINS=100  # Number of points in log space

# Set to "log" for logarithmic spacing, "linear" for linear spacing
N_COLORS_SCALE="log"
K_SELECTIONS_SCALE="log"

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXECUTABLE="$SCRIPT_DIR/finalMassProbabilities"

# Check if executable exists
if [ ! -f "$EXECUTABLE" ]; then
    echo "Error: Executable not found at $EXECUTABLE"
    echo "Please compile first with: g++ -std=c++17 -O3 -march=native finalMass.cpp -o finalMass"
    exit 1
fi

# Get total number of processors and leave 2 free
TOTAL_PROCS=$(nproc)
MAX_JOBS=$((TOTAL_PROCS - 2))
if [ $MAX_JOBS -lt 1 ]; then
    MAX_JOBS=1
fi

echo "Running simulations with $MAX_JOBS parallel jobs (leaving 2 cores free)..."
echo "N scale: $N_COLORS_SCALE, K scale: $K_SELECTIONS_SCALE"

# Function to generate log-spaced integers
generate_logspace() {
    local start=$1
    local end=$2
    local num_bins=$3
    python3 -c "import numpy as np; print(' '.join(map(str, np.unique(np.logspace(np.log10($start), np.log10($end), $num_bins, dtype=int)))))"
}

# Function to generate linear-spaced integers
generate_linspace() {
    local start=$1
    local end=$2
    local step=$3
    seq $start $step $end
}

# Generate N_colors values based on scale
if [ "$N_COLORS_SCALE" = "log" ]; then
    N_COLORS_VALUES=($(generate_logspace $N_COLORS_START $N_COLORS_END $N_COLORS_BINS))
else
    N_COLORS_VALUES=($(generate_linspace $N_COLORS_START $N_COLORS_END $N_COLORS_STEP))
fi

# Generate K_selections values based on scale
if [ "$K_SELECTIONS_SCALE" = "log" ]; then
    K_SELECTIONS_VALUES=($(generate_logspace $K_SELECTIONS_START $K_SELECTIONS_END $K_SELECTIONS_BINS))
else
    K_SELECTIONS_VALUES=($(generate_linspace $K_SELECTIONS_START $K_SELECTIONS_END $K_SELECTIONS_STEP))
fi

echo "N values: ${#N_COLORS_VALUES[@]} points from ${N_COLORS_VALUES[0]} to ${N_COLORS_VALUES[-1]}"
echo "K values: ${#K_SELECTIONS_VALUES[@]} points from ${K_SELECTIONS_VALUES[0]} to ${K_SELECTIONS_VALUES[-1]}"

# Create temporary file for parameter combinations
PARAM_FILE=$(mktemp)

# Generate all parameter combinations
for N_colors in "${N_COLORS_VALUES[@]}"; do
    for K_selections in "${K_SELECTIONS_VALUES[@]}"; do
        if [ $K_selections -le $INITIAL_TOTAL_BALLS ]; then
            echo "$N_colors $K_selections $N_STEPS $INITIAL_TOTAL_BALLS $EXPONENTIAL_RATE" >> "$PARAM_FILE"
        fi
    done
done

# Count total simulations
TOTAL_SIMS=$(wc -l < "$PARAM_FILE")
echo "Total simulations to run: $TOTAL_SIMS"

# Run simulations in parallel using GNU parallel if available, otherwise use xargs
if command -v parallel &> /dev/null; then
    echo "Using GNU parallel..."
    cat "$PARAM_FILE" | parallel --progress -j $MAX_JOBS --colsep ' ' "$EXECUTABLE {1} {2} {3} {4} {5}"
else
    echo "Using xargs (install GNU parallel for better progress reporting)..."
    
    # Create a simple progress tracker
    COMPLETED=0
    
    # Function to update progress
    update_progress() {
        COMPLETED=$((COMPLETED + 1))
        PERCENT=$((COMPLETED * 100 / TOTAL_SIMS))
        printf "\rProgress: %d/%d (%d%%) " $COMPLETED $TOTAL_SIMS $PERCENT
    }
    
    export -f update_progress
    export COMPLETED TOTAL_SIMS
    
    # Run with xargs and track completion
    cat "$PARAM_FILE" | while read -r params; do
        # Wait if we're at max jobs
        while [ $(jobs -r | wc -l) -ge $MAX_JOBS ]; do
            sleep 0.1
        done
        
        # Run simulation in background
        (
            $EXECUTABLE $params > /dev/null 2>&1
            update_progress
        ) &
    done
    
    # Wait for all background jobs to complete
    wait
    echo ""
fi

# Clean up
rm "$PARAM_FILE"

echo ""
echo "All simulations completed!"
echo "Results saved to: $SCRIPT_DIR/outputs/finalMass/"
