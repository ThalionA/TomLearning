#!/usr/bin/env bash
# Re-run run_lag_curves.py UNCAPPED (all running bins, whole session) for both FS
# conditions, 3 processes per condition split by animal (balanced on pair count),
# BLAS capped at 3 threads per process (one 10-thread process spent more time in
# system than in user — thread spin on 30x30 problems). Parts are merged back into
# results/lag_curves_bin10{,_fsincl}.csv by merge_lag_curve_parts.py at the end.
#
# Usage: nohup bash scripts/run_lag_curves_uncapped_batch.sh > results/lag_curves_uncapped.log 2>&1 &
set -uo pipefail
cd "$(dirname "$0")/.."
export PYTHONPATH=src
export OMP_NUM_THREADS=3 OPENBLAS_NUM_THREADS=3 MKL_NUM_THREADS=3 VECLIB_MAXIMUM_THREADS=3
G1="71,52,61,28"            # 8+6+6+4 = 24 pairs (FS-excl)
G2="66,70,73,34,36,63,98,41" # 6+6+6+2+1+1+1+2 = 25
G3="75,77,100,68"           # 6+6+6+4 = 22
COMMON="--bin-ms 10 --smooth-ms 2.5 --max-lag 25 --n-shuffles 200"
echo "START $(date)"
pids=()
for fs in excl incl; do
  flag=""; suf=""
  if [ "$fs" = incl ]; then flag="--include-fs"; suf="_fsincl"; fi
  i=0
  for grp in "$G1" "$G2" "$G3"; do
    i=$((i+1))
    python scripts/run_lag_curves.py $COMMON $flag --animals "$grp" \
        --out "lag_curves_bin10${suf}_part${i}" > "results/lag_curves_uncapped_${fs}_part${i}.log" 2>&1 &
    pids+=($!)
    echo "  launched fs=$fs part$i animals=$grp pid=$!"
  done
done
fail=0
for p in "${pids[@]}"; do wait "$p" || fail=1; done
echo "PARTS DONE $(date) fail=$fail"
python scripts/merge_lag_curve_parts.py || fail=1
echo "END $(date) fail=$fail"
