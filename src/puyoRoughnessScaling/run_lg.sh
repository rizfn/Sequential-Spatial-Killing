#!/usr/bin/env bash
cd "$(dirname "$0")"
: > lgw.log
{ for s in $(seq 0 7); do echo "2048 12 120000 $s 74000 48"; done
  for s in $(seq 0 7); do echo "2048 15 120000 $s 85000 48"; done
  for s in $(seq 0 7); do echo "2048 20 120000 $s 96000 48"; done
  for s in $(seq 0 7); do echo "2048 30 120000 $s 107000 48"; done
  for s in $(seq 0 7); do echo "2048 40 120000 $s 112500 48"; done; } | \
  xargs -P 4 -n6 sh -c './slopeDistFast "$1" "$2" "$3" "$4" "$5" "$6" >/dev/null 2>>lgw.log' _
echo "LG_RUNS_DONE ceiling=$(wc -l < lgw.log)"
