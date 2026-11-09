#!/usr/bin/env bash
# B2+B4 L3 GPU group scan: G in {32,64,128,256} ×3 (r01..r03).
# Does not touch src/. Requires rebuilt GPU binaries with gpu_paper_sycl_group.h.
# Usage: bash run_group_scan_b2_b4.sh
# Optional: ONLY=B2|B4  SCAN_GROUPS="32 64 128 256"  REPS="r01 r02 r03"
# NOTE: use SCAN_GROUPS (not GROUPS) — GROUPS is a common inherited env and caused G=1007.
set -eo pipefail
ROOT=/home/yongchuan/SPHinXsys-gpu-paper
R=$ROOT/docs/gpu_paper_discussion/scripts
ONLY=${ONLY:-ALL}
# Prefer SCAN_GROUPS; ignore bare GROUPS to avoid accidental pollution.
SCAN_GROUPS=${SCAN_GROUPS:-"32 64 128 256"}
REPS=${REPS:-"r01 r02 r03"}

# Clear stale group from parent shell so only the loop value applies.
unset GPU_PAPER_GROUP

run_one() {
  local case=$1 group=$2 rep=$3
  echo "==== ${case} L3 gpu_fp32 G=${group} ${rep} ===="
  case "$case" in
    B2)
      env -u GPU_PAPER_GROUP CUDA_VISIBLE_DEVICES=1 GPU_PAPER_GROUP=$group RUN_SUFFIX=$rep \
        bash "$R/run_b2.sh" L3 gpu_fp32 benchmark
      ;;
    B4)
      env -u GPU_PAPER_GROUP CUDA_VISIBLE_DEVICES=1 GPU_PAPER_GROUP=$group RUN_SUFFIX=$rep \
        bash "$R/run_b4.sh" L3 gpu_fp32 benchmark
      ;;
    *) echo "bad case $case"; exit 2 ;;
  esac
}

cases=()
[[ "$ONLY" == "ALL" || "$ONLY" == "B2" ]] && cases+=(B2)
[[ "$ONLY" == "ALL" || "$ONLY" == "B4" ]] && cases+=(B4)

echo "SCAN_GROUPS=$SCAN_GROUPS REPS=$REPS ONLY=$ONLY"
echo "WARNING: do not export GPU_PAPER_GROUP=1007 in the parent shell."

for c in "${cases[@]}"; do
  for g in $SCAN_GROUPS; do
    for rep in $REPS; do
      run_one "$c" "$g" "$rep"
    done
  done
done

echo "DONE group scan ONLY=$ONLY SCAN_GROUPS=$SCAN_GROUPS REPS=$REPS"
