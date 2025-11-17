#!/usr/bin/env bash
set -euo pipefail

STATS_DIR="${STATS_DIR:-perf_data/intervals/_summaries}"
CFG_LIST="${CFG_LIST:-fp01 fp02 fp03 fp04}"

extract_series() {
  local cfg="$1"    # e.g., fp02
  local key="$2"    # stalls_per_store | miss_rate | cycles_per_store | instructions_per_store | store_forward_blocks_per_store
  grep -h '^ratio_of_sums:' "${STATS_DIR}"/perf_"${cfg}"_rep*.stats | \
  awk -v K="$key" '{
    for (i=1;i<=NF;i++) if ($i ~ (K "=")) { split($i,a,"="); if (a[2] != "") print a[2] }
  }'
}

summ_stats() {
  awk '
  { n++; s+=$1; ss+=$1*$1 }
  END {
    if (n==0) { print "n=0"; exit 1 }
    m=s/n; sd=(n>1)?sqrt(ss/n - m*m):0
    printf("n=%d mean=%.6f sd=%.6f\n", n, m, sd)
  }'
}

for cfg in $CFG_LIST; do
  echo "=== ${cfg} ==="
  for metric in stalls_per_store miss_rate cycles_per_store instructions_per_store store_forward_blocks_per_store; do
    printf "%-25s" "$metric"
    extract_series "$cfg" "$metric" | summ_stats || echo "n=0"
  done
done
