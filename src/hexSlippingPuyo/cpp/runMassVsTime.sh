#!/bin/bash

# Fixed parameters
L=64
STEPS=1024

# Array of N_SPECIES values
N_SPECIES_VALUES=(2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25)

# Get number of CPUs and set max jobs (leave 2 cores free)
MAX_JOBS=$(( $(nproc) - 2 ))
if [ $MAX_JOBS -lt 1 ]; then
    MAX_JOBS=1
fi

running=0

# Run in parallel with limited concurrency
for N in "${N_SPECIES_VALUES[@]}"; do
    ./massVsTime2D $L $N $STEPS &
    ((running++))
    
    if [ $running -ge $MAX_JOBS ]; then
        wait -n  # Wait for one job to finish
        ((running--))
    fi
done

# Wait for all remaining jobs to finish
wait

echo "All runs completed."
