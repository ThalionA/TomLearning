#!/usr/bin/env bash
# Regenerate ALL landmark-arm figures from the current pkl/csv results.
#
# For each of the 44 landmark configs (25/50 ms x res/sig x k-rule) it runs the
# two per-config figure scripts, then runs the sweep summary once at the end.
# Stale figures on disk are overwritten. Progress is logged to stdout; per-config
# script output goes to results/.sandbox_scratch/regen_<tag>.log.
#
# Usage:  PYTHONPATH=src bash scripts/regen_all_landmark_figs.sh
set -uo pipefail

cd "$(dirname "$0")/.." || exit 1
export PYTHONPATH=src
LOGDIR=results/.sandbox_scratch
mkdir -p "$LOGDIR"

CONFIGS=()
while IFS= read -r line; do
  CONFIGS+=("$line")
done < <(ls results/landmark*_*.done \
  | sed 's#.*/##; s/\.done$//' \
  | grep -E '^landmark(25|50)_(res|sig)_' | sort)

n=${#CONFIGS[@]}
echo "[regen] $n configs to process"
fail=0
i=0
for tag in "${CONFIGS[@]}"; do
  i=$((i+1))
  log="$LOGDIR/regen_${tag}.log"
  printf '[regen] (%2d/%2d) %s ... ' "$i" "$n" "$tag"
  if python3 scripts/plot_landmark.py --tag "$tag" >"$log" 2>&1 \
     && python3 scripts/learning_changes.py --tag "$tag" >>"$log" 2>&1; then
    echo "ok"
  else
    echo "FAIL (see $log)"
    fail=$((fail+1))
  fi
done

echo "[regen] sweep summary across all done configs ..."
if python3 scripts/summarise_landmark_sweep.py >"$LOGDIR/regen_sweep.log" 2>&1; then
  echo "[regen] sweep ok"
else
  echo "[regen] sweep FAIL (see $LOGDIR/regen_sweep.log)"
  fail=$((fail+1))
fi

echo "[regen] done; failures=$fail"
exit $fail
