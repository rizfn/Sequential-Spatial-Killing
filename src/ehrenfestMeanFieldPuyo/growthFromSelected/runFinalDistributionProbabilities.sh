#!/bin/bash

# Define fixed parameters
N_COLORS=4000
K_SELECTIONS=30
N_STEPS=200000
INITIAL_TOTAL_BALLS=10000
EXPONENTIAL_RATE=-1

# Number of simulations to run
NUM_SIMS=10

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXECUTABLE="$SCRIPT_DIR/finalDistributionProbabilities"

# Check if executable exists
if [ ! -f "$EXECUTABLE" ]; then
    echo "Error: Executable not found at $EXECUTABLE"
    echo "Please compile first with: g++ -std=c++17 -O3 -march=native finalDistributionProbabilities.cpp -o finalDistributionProbabilities"
    exit 1
fi

# Get total number of processors and leave 2 free
TOTAL_PROCS=$(nproc)
MAX_JOBS=$((TOTAL_PROCS - 2))
if [ $MAX_JOBS -lt 1 ]; then
    MAX_JOBS=1
fi

echo "Running $NUM_SIMS simulations with $MAX_JOBS parallel jobs (leaving 2 cores free)..."
echo "Parameters: N=$N_COLORS, K=$K_SELECTIONS, steps=$N_STEPS, initial=$INITIAL_TOTAL_BALLS, rate=$EXPONENTIAL_RATE"
echo "Total simulations to run: $NUM_SIMS"

# Run simulations in parallel using GNU parallel if available, otherwise use xargs
if command -v parallel &> /dev/null; then
    echo "Using GNU parallel..."
    seq 0 $((NUM_SIMS-1)) | parallel --progress -j $MAX_JOBS \
        "$EXECUTABLE $N_COLORS $K_SELECTIONS $N_STEPS $INITIAL_TOTAL_BALLS $EXPONENTIAL_RATE {}"
else
    echo "Using xargs (install GNU parallel for better progress reporting)..."
    
    # Create a simple progress tracker
    COMPLETED=0
    
    # Function to update progress
    update_progress() {
        COMPLETED=$((COMPLETED + 1))
        PERCENT=$((COMPLETED * 100 / NUM_SIMS))
        printf "\rProgress: %d/%d (%d%%) " $COMPLETED $NUM_SIMS $PERCENT
    }
    
    export -f update_progress
    export COMPLETED NUM_SIMS
    
    # Run with xargs and track completion
    for ((sim=0; sim<$NUM_SIMS; sim++)); do
        # Wait if we're at max jobs
        while [ $(jobs -r | wc -l) -ge $MAX_JOBS ]; do
            sleep 0.1
        done
        
        # Run simulation in background
        (
            $EXECUTABLE $N_COLORS $K_SELECTIONS $N_STEPS $INITIAL_TOTAL_BALLS $EXPONENTIAL_RATE $sim > /dev/null 2>&1
            update_progress
        ) &
    done
    
    # Wait for all background jobs to complete
    wait
    echo ""
fi

echo ""
echo "All simulations completed!"
echo "Results saved to: $SCRIPT_DIR/outputs/finalDistributionProbabilities/"
