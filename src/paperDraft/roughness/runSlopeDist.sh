#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXE="$SCRIPT_DIR/slopeDist"

STEPS=8000
N_SIMS=15
STARTING_SIM_NO=${1:-0}

L_VALUES=(16 32 64 128)
N_VALUES=(10)

TOTAL_CPUS="$(getconf _NPROCESSORS_ONLN)"
CPUS_FREE=2
MAX_JOBS=$(( TOTAL_CPUS > CPUS_FREE ? TOTAL_CPUS - CPUS_FREE : 1 ))
TOTAL_SIMS=$(( ${#L_VALUES[@]} * ${#N_VALUES[@]} * N_SIMS ))

if [[ ! -x "$EXE" ]]; then
	echo "Executable not found or not executable: $EXE" >&2
	exit 1
fi

echo "Running $TOTAL_SIMS sims with max $MAX_JOBS concurrent jobs (keeping $CPUS_FREE CPUs free)"

echo "Sweep: ${#L_VALUES[@]} lattice sizes x ${#N_VALUES[@]} species counts x $N_SIMS sims each"

job_count=0
for L in "${L_VALUES[@]}"; do
	for N in "${N_VALUES[@]}"; do
		for sim in $(seq 0 $(( N_SIMS - 1 ))); do
			while (( job_count >= MAX_JOBS )); do
				wait -n
				job_count=$(( job_count - 1 ))
			done

			actual_sim=$(( STARTING_SIM_NO + sim ))
			"$EXE" "$L" "$N" "$STEPS" "$actual_sim" &
			job_count=$(( job_count + 1 ))
		 done
	done
done

wait

echo "All simulations completed."
