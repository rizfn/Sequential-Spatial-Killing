#!/usr/bin/env bash
cd "$(dirname "$0")"
: > n10bw.log
# L=1024, N=10: 2.4x more steps so it actually saturates. H=605k -> ~4.96GB/sim, so -P2.
for s in $(seq 0 11); do echo "1024 10 1200000 $s 605000 480"; done | \
  xargs -P 2 -n6 sh -c './slopeDistFast "$1" "$2" "$3" "$4" "$5" "$6" >/dev/null 2>>n10bw.log' _
echo "N10_L1024_LONG_DONE ceiling=$(wc -l < n10bw.log)"
