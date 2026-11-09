#!/usr/bin/env bash
# Usage: run_b1.sh <L1|L2|L3|L4|L5> <gpu_fp32|cpu_fp32|cpu_fp64> <benchmark|validate|showcase|smoke> [RELAX_ITERS]
# L_char/dp（B1 复杂几何从 L1=80 起；旧 20/40 兔子不成形，已废弃）:
#   L1=80, L2=160, L3=320, L4=640, L5=1280
# 默认 iterations: smoke=100, else=1000
# Geometry: frozen StanfordBunny_fixed.stl (sha256 a1b1badb...)
# GPU: CUDA_VISIBLE_DEVICES=1
set -eo pipefail
LEVEL=${1:?level L1..L5}
DEVICE=${2:?device gpu_fp32|cpu_fp32|cpu_fp64}
MODE=${3:?mode benchmark|validate|showcase|smoke}
ROOT=/home/yongchuan/SPHinXsys-gpu-paper
DISC=$ROOT/docs/gpu_paper_discussion
case "$LEVEL" in
  L1) DH=80 ;;
  L2) DH=160 ;;
  L3) DH=320 ;;
  L4) DH=640 ;;
  L5) DH=1280 ;;
  *) echo "bad LEVEL=$LEVEL"; exit 2 ;;
esac
case "$DEVICE" in
  gpu_fp32)
    BINDIR=$ROOT/build/tests/extra_source_and_tests/3d_examples/test_3d_gpu_paper_b1_particle_relaxation/bin
    NEED_ONEAPI=1
    ;;
  cpu_fp32)
    BINDIR=$ROOT/build_cpu_fp32/tests/extra_source_and_tests/3d_examples/test_3d_gpu_paper_b1_particle_relaxation/bin
    NEED_ONEAPI=0
    ;;
  cpu_fp64)
    BINDIR=$ROOT/build_cpu_fp64/tests/extra_source_and_tests/3d_examples/test_3d_gpu_paper_b1_particle_relaxation/bin
    NEED_ONEAPI=0
    ;;
  *) echo "bad DEVICE=$DEVICE"; exit 2 ;;
esac
BIN=$BINDIR/test_3d_gpu_paper_b1_particle_relaxation

VTP=0
case "$MODE" in
  benchmark)
    RELAX_ITERS=${4:-1000}
    SCHED=$DISC/results/B1_stanford_bunny/schedules/B1_benchmark_I1000_bound.json
    ;;
  validate)
    RELAX_ITERS=${4:-1000}
    SCHED=$DISC/results/B1_stanford_bunny/schedules/B1_validate_I1000_bound.json
    ;;
  showcase)
    RELAX_ITERS=${4:-1000}
    SCHED=$DISC/results/B1_stanford_bunny/schedules/B1_validate_I1000_bound.json
    if [[ "$DEVICE" == "gpu_fp32" ]]; then VTP=1; fi
    ;;
  smoke)
    RELAX_ITERS=${4:-100}
    SCHED=$DISC/results/B1_stanford_bunny/schedules/B1_smoke_I100_bound.json
    ;;
  *) echo "bad MODE=$MODE"; exit 2 ;;
esac

RUN_SUFFIX=${RUN_SUFFIX:-r01}
RUN_ID=B1_${LEVEL}_${DEVICE}_${MODE}_${RUN_SUFFIX}
RUN=$DISC/results/B1_stanford_bunny/runs/$RUN_ID
LOG=$DISC/logs/gp_B1_${LEVEL}_${DEVICE}_${MODE}_${RUN_SUFFIX}.log
mkdir -p "$RUN" "$(dirname "$LOG")"

if [[ ! -x "$BIN" ]]; then
  echo "missing binary: $BIN" | tee "$LOG"
  echo "EXIT:127" | tee -a "$LOG"
  exit 127
fi

# Ensure STL is present in bin/input (cmake copy; re-copy if missing)
if [[ ! -f "$BINDIR/input/StanfordBunny_fixed.stl" ]]; then
  mkdir -p "$BINDIR/input"
  cp -f "$ROOT/tests/extra_source_and_tests/3d_examples/test_3d_gpu_paper_b1_particle_relaxation/data/StanfordBunny_fixed.stl" \
    "$BINDIR/input/StanfordBunny_fixed.stl"
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
export GPU_PAPER_DH_OVER_DP=$DH
export GPU_PAPER_RELAX_ITERS=$RELAX_ITERS
export GPU_PAPER_ENABLE_VTP=$VTP
export GPU_PAPER_SCHEDULE_JSON=$SCHED

{
  echo "START:$(date -Is)"
  echo "RUN_ID=$RUN_ID LEVEL=$LEVEL DEVICE=$DEVICE MODE=$MODE L_over_dp=$DH VTP=$VTP RELAX_ITERS=$RELAX_ITERS"
  echo "BIN=$BIN"
  echo "STL=$BINDIR/input/StanfordBunny_fixed.stl sha256=$(sha256sum "$BINDIR/input/StanfordBunny_fixed.stl" | awk '{print $1}')"
  echo "NOTE: B1 Stanford bunny unmodified STL solid body-fitted + staged timing (levelset/particle_gen/relax); coordinate=iteration"
  ls -la "$BIN"
} | tee "$LOG"

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
exit "$rc"
