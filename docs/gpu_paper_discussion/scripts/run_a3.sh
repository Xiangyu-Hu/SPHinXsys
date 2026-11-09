#!/usr/bin/env bash
# Usage: run_a3.sh <L1|L2|L3> <gpu_fp32|cpu_fp32|cpu_fp64> <smoke|validate> [END_TIME]
# A3 D4 windows-frame: resolution scale L1=1,L2=2,L3=4；validate 默认 End_Time=0.02（伪时间）
# CK path: DiffusionRelaxationCK + Robin；GPU: CUDA_VISIBLE_DEVICES=1
set -eo pipefail
LEVEL=${1:?level L1/L2/L3}
DEVICE=${2:?device gpu_fp32|cpu_fp32|cpu_fp64}
MODE=${3:?mode smoke|validate}
ROOT=/home/yongchuan/SPHinXsys-gpu-paper
DISC=$ROOT/docs/gpu_paper_discussion
case "$LEVEL" in
  L1) DH=1 ;;
  L2) DH=2 ;;
  L3) DH=4 ;;
  *) echo "bad LEVEL=$LEVEL"; exit 2 ;;
esac
case "$DEVICE" in
  gpu_fp32)
    BIN=$ROOT/build/tests/extra_source_and_tests/2d_examples/test_2d_gpu_paper_a3_windows_frame_D4/bin/test_2d_gpu_paper_a3_windows_frame_D4
    NEED_ONEAPI=1
    ;;
  cpu_fp32)
    BIN=$ROOT/build_cpu_fp32/tests/extra_source_and_tests/2d_examples/test_2d_gpu_paper_a3_windows_frame_D4/bin/test_2d_gpu_paper_a3_windows_frame_D4
    NEED_ONEAPI=0
    ;;
  cpu_fp64)
    BIN=$ROOT/build_cpu_fp64/tests/extra_source_and_tests/2d_examples/test_2d_gpu_paper_a3_windows_frame_D4/bin/test_2d_gpu_paper_a3_windows_frame_D4
    NEED_ONEAPI=0
    ;;
  *) echo "bad DEVICE=$DEVICE"; exit 2 ;;
esac

VTP=0
case "$MODE" in
  smoke)
    END_TIME=${4:-0.002}
    SCHED=$DISC/results/A3_windows_frame_D4/schedules/A3_validate_T0p02_bound.json
    ;;
  validate)
    END_TIME=${4:-0.02}
    SCHED=$DISC/results/A3_windows_frame_D4/schedules/A3_validate_T0p02_bound.json
    if [[ "$DEVICE" == "gpu_fp32" && "$LEVEL" == "L2" ]]; then VTP=1; fi
    ;;
  *) echo "bad MODE=$MODE (use smoke|validate)"; exit 2 ;;
esac

RUN_SUFFIX=${RUN_SUFFIX:-r01}
RUN_ID=A3_D4_${LEVEL}_${DEVICE}_${MODE}_${RUN_SUFFIX}
RUN=$DISC/results/A3_windows_frame_D4/runs/$RUN_ID
LOG=$DISC/logs/gp_A3_D4_${LEVEL}_${DEVICE}_${MODE}_${RUN_SUFFIX}.log
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

{
  echo "START:$(date -Is)"
  echo "RUN_ID=$RUN_ID LEVEL=$LEVEL DEVICE=$DEVICE MODE=$MODE res_scale=$DH VTP=$VTP End_Time=$END_TIME(pseudo)"
  echo "BIN=$BIN"
  echo "RUN_DIR=$RUN"
  echo "NOTE: A3 CK path DiffusionRelaxationCK+Robin; backend sycl_ck_gpu or cpu_ck"
  echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-unset}"
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
