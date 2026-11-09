#!/usr/bin/env bash
# Usage: run_b4.sh <L1|L2|L3|L4|L5> <gpu_fp32|cpu_fp32|cpu_fp64> <benchmark|validate|showcase|smoke> [END_TIME]
# PH/dp: L1=6 (source), L2=12, L3=24, L4=48, L5=96；默认 End_Time=1.0；smoke 默认 0.1
# Geometry/material: NeoHookean passive cantilever (classic test_3d_passive_cantilever_neohookean)
# Note: classic physical damping omitted (no CK); gravity via GravityForceCK
# GPU: CUDA_VISIBLE_DEVICES=1
set -eo pipefail
LEVEL=${1:?level L1..L5}
DEVICE=${2:?device gpu_fp32|cpu_fp32|cpu_fp64}
MODE=${3:?mode benchmark|validate|showcase|smoke}
ROOT=/home/yongchuan/SPHinXsys-gpu-paper
DISC=$ROOT/docs/gpu_paper_discussion
case "$LEVEL" in
  L1) DH=6 ;;
  L2) DH=12 ;;
  L3) DH=24 ;;
  L4) DH=48 ;;
  L5) DH=96 ;;
  *) echo "bad LEVEL=$LEVEL"; exit 2 ;;
esac
case "$DEVICE" in
  gpu_fp32)
    BIN=$ROOT/build/tests/extra_source_and_tests/3d_examples/test_3d_gpu_paper_b4_passive_cantilever/bin/test_3d_gpu_paper_b4_passive_cantilever
    NEED_ONEAPI=1
    ;;
  cpu_fp32)
    BIN=$ROOT/build_cpu_fp32/tests/extra_source_and_tests/3d_examples/test_3d_gpu_paper_b4_passive_cantilever/bin/test_3d_gpu_paper_b4_passive_cantilever
    NEED_ONEAPI=0
    ;;
  cpu_fp64)
    BIN=$ROOT/build_cpu_fp64/tests/extra_source_and_tests/3d_examples/test_3d_gpu_paper_b4_passive_cantilever/bin/test_3d_gpu_paper_b4_passive_cantilever
    NEED_ONEAPI=0
    ;;
  *) echo "bad DEVICE=$DEVICE"; exit 2 ;;
esac

VTP=0
case "$MODE" in
  benchmark)
    END_TIME=${4:-1.0}
    SCHED=$DISC/results/B4_passive_cantilever/schedules/B4_benchmark_T1_bound.json
    ;;
  validate)
    END_TIME=${4:-1.0}
    SCHED=$DISC/results/B4_passive_cantilever/schedules/B4_validate_T1_bound.json
    ;;
  showcase)
    END_TIME=${4:-1.0}
    SCHED=$DISC/results/B4_passive_cantilever/schedules/B4_validate_T1_bound.json
    if [[ "$DEVICE" == "gpu_fp32" ]]; then VTP=1; fi
    ;;
  smoke)
    END_TIME=${4:-0.1}
    SCHED=$DISC/results/B4_passive_cantilever/schedules/B4_smoke_T0p1_bound.json
    ;;
  *) echo "bad MODE=$MODE"; exit 2 ;;
esac

RUN_SUFFIX=${RUN_SUFFIX:-r01}
GROUP_TAG=""
if [[ -n "${GPU_PAPER_GROUP:-}" ]]; then
  GROUP_TAG=_g${GPU_PAPER_GROUP}
fi
RUN_ID=B4_${LEVEL}_${DEVICE}_${MODE}${GROUP_TAG}_${RUN_SUFFIX}
RUN=$DISC/results/B4_passive_cantilever/runs/$RUN_ID
LOG=$DISC/logs/gp_B4_${LEVEL}_${DEVICE}_${MODE}${GROUP_TAG}_${RUN_SUFFIX}.log
mkdir -p "$RUN" "$(dirname "$LOG")"

if [[ ! -x "$BIN" ]]; then
  echo "missing binary: $BIN" | tee "$LOG"
  echo "EXIT:127" | tee -a "$LOG"
  exit 127
fi

if [[ "$NEED_ONEAPI" == "1" ]]; then
  set +u
  # shellcheck disable=SC1091
  source /opt/intel/oneapi/setvars.sh --force >/dev/null 2>&1 || true
  set -eo pipefail
  export CUDA_VISIBLE_DEVICES=1
fi

export GPU_PAPER_RUN_ID=$RUN_ID
export GPU_PAPER_MODE=$MODE
export GPU_PAPER_RUN_DIR=$RUN
export GPU_PAPER_END_TIME=$END_TIME
export GPU_PAPER_DH_OVER_DP=$DH
export GPU_PAPER_ENABLE_VTP=$VTP
export GPU_PAPER_SCHEDULE_JSON=$SCHED
# optional: GPU_PAPER_GROUP=32|64|128|256 (SYCL work-group); leave unset for default

{
  echo "START:$(date -Is)"
  echo "RUN_ID=$RUN_ID LEVEL=$LEVEL DEVICE=$DEVICE MODE=$MODE PH_over_dp=$DH VTP=$VTP End_Time=$END_TIME GPU_PAPER_GROUP=${GPU_PAPER_GROUP:-default}"
  echo "BIN=$BIN"
  echo "NOTE: NeoHookean cantilever; physical damping omitted (no CK); gravity on"
  ls -la "$BIN"
} | tee "$LOG"

cd "$(dirname "$BIN")"
set +e
if command -v stdbuf >/dev/null 2>&1; then
  stdbuf -oL -eL "$BIN" 2>&1 | tee -a "$LOG"
  rc=${PIPESTATUS[0]}
else
  "$BIN" 2>&1 | tee -a "$LOG"
  rc=${PIPESTATUS[0]}
fi
set -e
echo "EXIT:${rc}" | tee -a "$LOG"
exit "$rc"
