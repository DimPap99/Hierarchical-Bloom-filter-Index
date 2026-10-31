#!/usr/bin/env bash
set -euo pipefail

FP_LIST="${FP_LIST:-0.2}"
REPS="${REPS:-5}"
SAMPLE_MS="${SAMPLE_MS:-250}"
CORE="${CORE:-3}"
OUTDIR="${OUTDIR:-perf_data/intervals}"

JVM_FLAGS='-XX:+UnlockDiagnosticVMOptions -XX:+PreserveFramePointer -XX:+DebugNonSafepoints'
CP_FLAGS='-cp target/Hierarchical-Bloom-filter-Index-1.0-SNAPSHOT.jar'
MAIN_CLASS='HBIDatasetBenchmark'
APP_COMMON='--runs 1 --run-suffix false --skip-queries true'

EVENTS='cycles,instructions,cpu_core/mem_inst_retired.all_stores/,cpu_core/resource_stalls.sb/,cpu_core/ld_blocks.store_forward/,cpu_core/cache-references/,cpu_core/cache-misses/'

mkdir -p "$OUTDIR"
echo "Output directory: $OUTDIR"
echo "Sampling period: ${SAMPLE_MS} ms. Core: $CORE"
echo "False positive rates: $FP_LIST. Repetitions: $REPS"
echo "Events: $EVENTS"

for FP in $FP_LIST; do
  echo "Pre-warming Java for fp=$FP"
  taskset -c "$CORE" java $JVM_FLAGS $CP_FLAGS $MAIN_CLASS $APP_COMMON --fp "$FP" >/dev/null 2>&1 || true

  for r in $(seq 1 "$REPS"); do
    out="$OUTDIR/perf_fp${FP/./}_rep${r}.csv"
    echo "Running fp=$FP rep=$r → $out"
    sudo perf stat -I "$SAMPLE_MS" -x, --no-big-num \
      -e "$EVENTS" \
      -- taskset -c "$CORE" \
         java $JVM_FLAGS $CP_FLAGS $MAIN_CLASS $APP_COMMON --fp "$FP" \
      2> "$out"
  done
done

echo "Done. CSV files are under $OUTDIR."
