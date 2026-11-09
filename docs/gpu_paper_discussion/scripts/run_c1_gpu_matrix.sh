#!/usr/bin/env bash
# C1 GPU matrix: require per-dp relaxed solids, then smoke → benchmark ×3.
# Usage:
#   bash run_c1_gpu_matrix.sh              # uses RESOLUTION=0.0012 store
#   SKIP_SMOKE=1 bash run_c1_gpu_matrix.sh
#   RESOLUTION=0.0012 bash run_c1_gpu_matrix.sh
set -eo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT=${GPU_PAPER_ROOT:-$(cd "$SCRIPT_DIR/../../.." && pwd)}
SCRIPTS=$ROOT/docs/gpu_paper_discussion/scripts
export CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-1}
export RESOLUTION=${RESOLUTION:-0.0012}
export GPU_PAPER_ROOT=$ROOT
# EINSIMO_ROOT is consumed by run_c1.sh; keep caller's value if set.

DP_TAG=$(python3 - <<PY
r=float("$RESOLUTION")
s=f"{r:.4f}".rstrip("0").rstrip(".")
print(s.replace(".","p"))
PY
)
STORE=$ROOT/docs/gpu_paper_discussion/results/C1_gearbox/reload/dp${DP_TAG}

cd "$ROOT"
mkdir -p "$ROOT/docs/gpu_paper_discussion/logs"

{
  echo "===== C1 GPU MATRIX START $(date -Is) CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES RESOLUTION=$RESOLUTION ====="

  for f in Gear1_rld.xml Gear2_rld.xml Gear3_rld.xml Wall_rld.xml; do
    if [[ ! -f "$STORE/$f" ]]; then
      echo "ERROR: missing $STORE/$f"
      echo "First run: RESOLUTION=$RESOLUTION bash $SCRIPTS/run_c1.sh gpu_fp32 relax"
      exit 2
    fi
  done
  echo "Using relaxed solids from $STORE"

  if [[ "${SKIP_SMOKE:-0}" != "1" ]]; then
    echo "===== smoke r01 ====="
    RUN_SUFFIX=r01 bash "$SCRIPTS/run_c1.sh" gpu_fp32 smoke 0.01
  else
    echo "===== skip smoke ====="
  fi

  for suf in r01 r02 r03; do
    echo "===== benchmark $suf End_Time=1.05 ====="
    RUN_SUFFIX=$suf bash "$SCRIPTS/run_c1.sh" gpu_fp32 benchmark 1.05
  done

  echo "===== C1 GPU MATRIX DONE $(date -Is) ====="
} 2>&1 | tee -a "$ROOT/docs/gpu_paper_discussion/logs/gp_C1_gpu_matrix.log"
