#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXE="$SCRIPT_DIR/massVsTime2D"

L=128
STEPS=4096
N_SIMS=4
N_VALUES=(8 8.2 8.5 9 9.2 9.5 10)

TOTAL_CPUS="$(getconf _NPROCESSORS_ONLN)"
MAX_JOBS=$(( TOTAL_CPUS > 2 ? TOTAL_CPUS - 2 : 1 ))

if [[ ! -x "$EXE" ]]; then
	echo "Executable not found or not executable: $EXE" >&2
	exit 1
fi

echo "Running ${#N_VALUES[@]} N values x $N_SIMS sims with max $MAX_JOBS concurrent jobs (keeping 2 CPUs free)..."

job_count=0
for N in "${N_VALUES[@]}"; do
	for sim in $(seq 0 $(( N_SIMS - 1 ))); do
		while (( job_count >= MAX_JOBS )); do
			wait -n
			job_count=$(( job_count - 1 ))
		done

		"$EXE" "$L" "$N" "$STEPS" "$sim" &
		job_count=$(( job_count + 1 ))
	done
done

wait

echo "All simulations completed."
