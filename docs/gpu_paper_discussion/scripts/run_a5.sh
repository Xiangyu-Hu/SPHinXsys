#!/usr/bin/env bash
# Usage: run_a5.sh <L1|L2|L3> <gpu_fp32|cpu_fp32|cpu_fp64> <relax|smoke|validate> [END_TIME]
# D/dp: L1=10,L2=20,L3=40
# Workflow: 先 relax 写出 bin/reload/Reload.xml，再 smoke/validate（默认 RELOAD=1）
# GPU: CUDA_VISIBLE_DEVICES=1
#
# Shared Reload (CPU/GPU alignment):
#   GPU_PAPER_PUBLISH_SHARED_RELOAD=1  — after successful relax, copy Reload.xml to shared_reload/<LEVEL>/
#   GPU_PAPER_SHARED_RELOAD=1          — before smoke/validate, install shared Reload.xml into this BINDIR
#   GPU_PAPER_SHARED_RELOAD_SRC=/path  — optional explicit source xml (overrides shared_reload/<LEVEL>)
set -eo pipefail
LEVEL=${1:?level L1/L2/L3}
DEVICE=${2:?device gpu_fp32|cpu_fp32|cpu_fp64}
MODE=${3:?mode relax|smoke|validate}
ROOT=/home/yongchuan/SPHinXsys-gpu-paper
DISC=$ROOT/docs/gpu_paper_discussion
SHARED_RELOAD_DIR=$DISC/results/A5_fsi2/shared_reload/$LEVEL
SHARED_RELOAD_XML=${GPU_PAPER_SHARED_RELOAD_SRC:-$SHARED_RELOAD_DIR/Reload.xml}
case "$LEVEL" in
  L1) DH=10 ;;
  L2) DH=20 ;;
  L3) DH=40 ;;
  *) echo "bad LEVEL=$LEVEL"; exit 2 ;;
esac
case "$DEVICE" in
  gpu_fp32)
    BINDIR=$ROOT/build/tests/extra_source_and_tests/2d_examples/test_2d_gpu_paper_a5_fsi2/bin
    NEED_ONEAPI=1
    ;;
  cpu_fp32)
    BINDIR=$ROOT/build_cpu_fp32/tests/extra_source_and_tests/2d_examples/test_2d_gpu_paper_a5_fsi2/bin
    NEED_ONEAPI=0
    ;;
  cpu_fp64)
    BINDIR=$ROOT/build_cpu_fp64/tests/extra_source_and_tests/2d_examples/test_2d_gpu_paper_a5_fsi2/bin
    NEED_ONEAPI=0
    ;;
  *) echo "bad DEVICE=$DEVICE"; exit 2 ;;
esac
BIN=$BINDIR/test_2d_gpu_paper_a5_fsi2
RELOAD_XML=$BINDIR/reload/Reload.xml

VTP=0
RELAX=0
RELOAD=${GPU_PAPER_RELOAD:-1}
USE_SHARED=${GPU_PAPER_SHARED_RELOAD:-0}
PUBLISH_SHARED=${GPU_PAPER_PUBLISH_SHARED_RELOAD:-0}
# ALIGN=1: shared Reload + no particle sort + fixed acoustic dt (CPU/GPU tip alignment)
ALIGN=${GPU_PAPER_ALIGN:-0}
case "$MODE" in
  relax)
    RELAX=1
    RELOAD=0
    END_TIME=${4:-0}
    SCHED=$DISC/results/A5_fsi2/schedules/A5_smoke_T2_bound.json
    ;;
  smoke)
    END_TIME=${4:-2}
    SCHED=$DISC/results/A5_fsi2/schedules/A5_smoke_T2_bound.json
    ;;
  validate)
    END_TIME=${4:-200}
    SCHED=$DISC/results/A5_fsi2/schedules/A5_validate_T200_bound.json
    if [[ "$DEVICE" == "gpu_fp32" && "$LEVEL" == "L2" ]]; then VTP=1; fi
    ;;
  *) echo "bad MODE=$MODE (use relax|smoke|validate)"; exit 2 ;;
esac
if [[ "$ALIGN" == "1" && "$MODE" != "relax" ]]; then
  USE_SHARED=1
  export GPU_PAPER_SHARED_RELOAD=1
  export GPU_PAPER_DISABLE_PARTICLE_SORT=${GPU_PAPER_DISABLE_PARTICLE_SORT:-1}
  # L1: adaptive acoustic is ~O(1e-3); 0.003 keeps similar cost while syncing time grids
  export GPU_PAPER_FIXED_ACOUSTIC_DT=${GPU_PAPER_FIXED_ACOUSTIC_DT:-0.003}
  if [[ "$DEVICE" == cpu_fp32 || "$DEVICE" == cpu_fp64 ]]; then
    export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
  fi
fi

RUN_SUFFIX=${RUN_SUFFIX:-r01}
RUN_ID=A5_${LEVEL}_${DEVICE}_${MODE}_${RUN_SUFFIX}
RUN=$DISC/results/A5_fsi2/runs/$RUN_ID
LOG=$DISC/logs/gp_A5_${LEVEL}_${DEVICE}_${MODE}_${RUN_SUFFIX}.log
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

install_shared_reload() {
  if [[ ! -f "$SHARED_RELOAD_XML" ]]; then
    echo "ERROR: shared Reload missing: $SHARED_RELOAD_XML" | tee -a "$LOG"
    echo "HINT: first run gpu_fp32 relax with GPU_PAPER_PUBLISH_SHARED_RELOAD=1" | tee -a "$LOG"
    echo "EXIT:2" | tee -a "$LOG"
    exit 2
  fi
  mkdir -p "$BINDIR/reload"
  cp -f "$SHARED_RELOAD_XML" "$RELOAD_XML"
  echo "INSTALLED_SHARED_RELOAD src=$SHARED_RELOAD_XML -> $RELOAD_XML sha256=$(sha256sum "$RELOAD_XML" | awk '{print $1}')" | tee -a "$LOG"
}

# Align mode: force identical InsertedBody IC from shared Reload.xml
if [[ "$MODE" != "relax" && "$RELOAD" == "1" && "$USE_SHARED" == "1" ]]; then
  : >"$LOG"
  install_shared_reload
fi

# smoke/validate with reload: auto-run relax for this resolution if Reload.xml missing
if [[ "$MODE" != "relax" && "$RELOAD" == "1" && ! -f "$RELOAD_XML" ]]; then
  echo "NOTE: $RELOAD_XML missing → auto-running relax for D/dp=$DH first" | tee -a "$LOG"
  GPU_PAPER_RELAX=1 GPU_PAPER_RELOAD=0 GPU_PAPER_DH_OVER_DP=$DH \
    bash "$0" "$LEVEL" "$DEVICE" relax 2>&1 | tee -a "$LOG"
  if [[ ! -f "$RELOAD_XML" ]]; then
    echo "ERROR: relax did not produce $RELOAD_XML" | tee -a "$LOG"
    echo "EXIT:2" | tee -a "$LOG"
    exit 2
  fi
fi

export GPU_PAPER_RUN_ID=$RUN_ID
export GPU_PAPER_MODE=$MODE
export GPU_PAPER_RUN_DIR=$RUN
export GPU_PAPER_END_TIME=$END_TIME
export GPU_PAPER_DH_OVER_DP=$DH
export GPU_PAPER_ENABLE_VTP=$VTP
export GPU_PAPER_SCHEDULE_JSON=$SCHED
export GPU_PAPER_RELAX=$RELAX
export GPU_PAPER_RELOAD=$RELOAD

{
  echo "START:$(date -Is)"
  echo "RUN_ID=$RUN_ID LEVEL=$LEVEL DEVICE=$DEVICE MODE=$MODE D_over_dp=$DH VTP=$VTP End_Time=$END_TIME RELAX=$RELAX RELOAD=$RELOAD USE_SHARED=$USE_SHARED PUBLISH_SHARED=$PUBLISH_SHARED ALIGN=$ALIGN"
  echo "DISABLE_SORT=${GPU_PAPER_DISABLE_PARTICLE_SORT:-0} FIXED_ACOUSTIC_DT=${GPU_PAPER_FIXED_ACOUSTIC_DT:-adaptive} OMP_NUM_THREADS=${OMP_NUM_THREADS:-unset}"
  echo "BIN=$BIN"
  echo "RELOAD_XML=$RELOAD_XML exists=$([[ -f $RELOAD_XML ]] && echo yes || echo no)"
  if [[ -f "$RELOAD_XML" ]]; then
    echo "RELOAD_SHA256=$(sha256sum "$RELOAD_XML" | awk '{print $1}')"
  fi
  echo "NOTE: SYCL/CK FSI2; insert body relax→reload; emitter/disposer BC"
  echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-unset}"
  ls -la "$BIN"
} | tee -a "$LOG"

cd "$BINDIR"
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
if [[ "$MODE" == "relax" && "$rc" -eq 0 ]]; then
  ls -la "$BINDIR/reload" 2>&1 | tee -a "$LOG" || true
  if [[ "$PUBLISH_SHARED" == "1" && -f "$RELOAD_XML" ]]; then
    mkdir -p "$SHARED_RELOAD_DIR"
    cp -f "$RELOAD_XML" "$SHARED_RELOAD_XML"
    echo "PUBLISHED_SHARED_RELOAD $RELOAD_XML -> $SHARED_RELOAD_XML sha256=$(sha256sum "$SHARED_RELOAD_XML" | awk '{print $1}')" | tee -a "$LOG"
  fi
fi
exit "$rc"
