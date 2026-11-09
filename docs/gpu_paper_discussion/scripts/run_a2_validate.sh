#!/usr/bin/env bash
# Usage: run_a2_validate.sh <L1|L2|L3> <gpu_fp32|cpu_fp32|cpu_fp64>
# stdout+log via tee — tmux attach 可见 N=/Time=
set -eo pipefail
LEVEL=${1:?level L1/L2/L3}
DEVICE=${2:?device gpu_fp32|cpu_fp32|cpu_fp64}
ROOT=/home/yongchuan/SPHinXsys-gpu-paper
DISC=$ROOT/docs/gpu_paper_discussion
case "$LEVEL" in
  L1) DH=10 ;;
  L2) DH=20 ;;
  L3) DH=40 ;;
  *) echo "bad LEVEL=$LEVEL"; exit 2 ;;
esac
case "$DEVICE" in
  gpu_fp32)
    BIN=$ROOT/build/tests/extra_source_and_tests/2d_examples/test_2d_gpu_paper_a2_oscillating_beam/bin/test_2d_gpu_paper_a2_oscillating_beam
    if [[ "$LEVEL" == "L2" ]]; then VTP=1; else VTP=0; fi
    NEED_ONEAPI=1
    ;;
  cpu_fp32)
    BIN=$ROOT/build_cpu_fp32/tests/extra_source_and_tests/2d_examples/test_2d_gpu_paper_a2_oscillating_beam/bin/test_2d_gpu_paper_a2_oscillating_beam
    VTP=0
    NEED_ONEAPI=0
    ;;
  cpu_fp64)
    BIN=$ROOT/build_cpu_fp64/tests/extra_source_and_tests/2d_examples/test_2d_gpu_paper_a2_oscillating_beam/bin/test_2d_gpu_paper_a2_oscillating_beam
    VTP=0
    NEED_ONEAPI=0
    ;;
  *) echo "bad DEVICE=$DEVICE"; exit 2 ;;
esac

RUN_SUFFIX=${RUN_SUFFIX:-r01}
RUN_ID=A2_${LEVEL}_${DEVICE}_validate_${RUN_SUFFIX}
RUN=$DISC/results/A2_oscillating_beam/runs/$RUN_ID
LOG=$DISC/logs/gp_A2_${LEVEL}_${DEVICE}_${RUN_SUFFIX}.log
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
export GPU_PAPER_MODE=validate
export GPU_PAPER_RUN_DIR=$RUN
export GPU_PAPER_END_TIME=2.0
export GPU_PAPER_DH_OVER_DP=$DH
export GPU_PAPER_ENABLE_VTP=$VTP
export GPU_PAPER_SCHEDULE_JSON=$DISC/results/A2_oscillating_beam/schedules/A2_validate_T2_bound.json

{
  echo "START:$(date -Is)"
  echo "RUN_ID=$RUN_ID LEVEL=$LEVEL DEVICE=$DEVICE PH/dp=$DH VTP=$VTP"
  echo "BIN=$BIN"
  echo "RUN_DIR=$RUN"
  echo "SCHEDULE=$GPU_PAPER_SCHEDULE_JSON"
  echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-}"
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
