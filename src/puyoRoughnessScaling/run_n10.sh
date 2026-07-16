#!/usr/bin/env bash
cd "$(dirname "$0")"
: > n10w.log
# L=256 (158MB/sim): extend to 32 sims
for s in $(seq 16 31); do echo "256 10 150000 $s 77000 60"; done | \
  xargs -P 10 -n6 sh -c './slopeDistFast "$1" "$2" "$3" "$4" "$5" "$6" >/dev/null 2>>n10w.log' _
echo "L=256 done"
# L=512 (520MB/sim): extend to 24 sims
for s in $(seq 12 23); do echo "512 10 250000 $s 127000 100"; done | \
  xargs -P 6 -n6 sh -c './slopeDistFast "$1" "$2" "$3" "$4" "$5" "$6" >/dev/null 2>>n10w.log' _
echo "L=512 done"
# L=1024 (2.03GB/sim): 16 sims, low concurrency for RAM
for s in $(seq 0 15); do echo "1024 10 500000 $s 253000 200"; done | \
  xargs -P 4 -n6 sh -c './slopeDistFast "$1" "$2" "$3" "$4" "$5" "$6" >/dev/null 2>>n10w.log' _
echo "L=1024 done"
echo "ALL N=10 RUNS COMPLETE; ceiling warnings: $(wc -l < n10w.log)"
