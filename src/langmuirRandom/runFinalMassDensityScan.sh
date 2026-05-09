#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXE="$SCRIPT_DIR/finalMass"

L=64
N=128
STEPS=1000
N_SIMS=20
STARTING_SIM_NO=${1:-0}
DENSITIES=(0.00 0.01 0.02 0.03 0.04 0.05 
		   0.06 0.07 0.08 0.09 0.10
		   0.11 0.12 0.13 0.14 0.15 
		   0.16 0.17 0.18 0.19 0.20
		   0.21 0.22 0.23 0.24 0.25
		   0.26 0.27 0.28 0.29 0.30
		   0.31 0.32 0.33 0.34 0.35
		   0.36 0.37 0.38 0.39 0.40
		   0.41 0.42 0.43 0.44 0.45
		   0.46 0.47 0.48 0.49 0.50)

TOTAL_CPUS="$(getconf _NPROCESSORS_ONLN)"
MAX_JOBS=$(( TOTAL_CPUS > 2 ? TOTAL_CPUS - 2 : 1 ))

if [[ ! -x "$EXE" ]]; then
	echo "Executable not found or not executable: $EXE" >&2
	exit 1
fi

echo "Running ${#DENSITIES[@]} densities x $N_SIMS sims (starting from sim $STARTING_SIM_NO) with max $MAX_JOBS concurrent jobs (keeping 2 CPUs free)..."

job_count=0
for density in "${DENSITIES[@]}"; do
	for sim in $(seq 0 $(( N_SIMS - 1 ))); do
		while (( job_count >= MAX_JOBS )); do
			wait -n
			job_count=$(( job_count - 1 ))
		done

		actual_sim=$(( STARTING_SIM_NO + sim ))
		"$EXE" "$L" "$N" "$STEPS" "$density" "$actual_sim" &
		job_count=$(( job_count + 1 ))
	done
done

wait

echo "All simulations completed."