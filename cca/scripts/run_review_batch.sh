#!/usr/bin/env bash
# Review re-runs (2026-06-11): trajectory at window=15 and the upgraded KCCA
# (30 shuffles, +-8 lags), each at BOTH bin widths (25 & 50 ms) and BOTH FS
# conditions, to distinct output files. Long-running; launch detached with nohup.
set -u
cd /Users/theoamvr/Desktop/Experiments/TomLearning/cca || exit 1
export PYTHONPATH=src
LOG=results/review_batch.log
: > "$LOG"
echo "=== review batch START $(date) ===" | tee -a "$LOG"

run() {
  echo ">>> $(date +%H:%M:%S)  $*" | tee -a "$LOG"
  python "$@" >> "$LOG" 2>&1
  echo "<<< exit $? $(date +%H:%M:%S)" | tee -a "$LOG"
}

# --- Trajectory, window=15 step=5, both bins, both FS ---
for bin in 25 50; do
  run scripts/run_trajectory.py --bin-ms "$bin" --window 15 --step 5 --min-window 15 --out "trajectory_w15_bin${bin}"
  run scripts/run_trajectory.py --bin-ms "$bin" --window 15 --step 5 --min-window 15 --include-fs --out "trajectory_w15_bin${bin}"
done

# --- KCCA upgraded (30 shuffles, +-8 lags), both bins, both FS ---
for bin in 25 50; do
  run scripts/run_kcca.py --bin-ms "$bin" --out "kcca_up_bin${bin}"
  run scripts/run_kcca.py --bin-ms "$bin" --include-fs --out "kcca_up_bin${bin}_fsincl"
done

echo "=== review batch DONE $(date) ===" | tee -a "$LOG"
