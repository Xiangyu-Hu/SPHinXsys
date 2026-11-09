#!/usr/bin/env bash
# C1 industrial gearbox (EinSIMO reducer) — run in-place, archive into gpu-paper.
# Usage: run_c1.sh <gpu_fp32|cpu_fp32> <relax|smoke|benchmark|showcase> [END_TIME]
#
# Particles:
#   Solids (Gear/Wall) MUST be relaxed per resolution, then reloaded.
#   Oil is Lattice-generated at runtime (case design).
#   Relax wall time is a separate process and is NOT in solver timer_total.
#   Per-dp store: results/C1_gearbox/reload/dp<tag>/  (e.g. dp0p0012)
#
# Env: RESOLUTION=0.0012 CASE_ID=8 RUN_SUFFIX=r01 ENABLE_CONTACT_FACTOR=0
set -eo pipefail

DEVICE=${1:?device gpu_fp32|cpu_fp32}
MODE=${2:?mode relax|smoke|benchmark|showcase}
END_TIME=${3:-}

# Portable roots (override on other machines):
#   export GPU_PAPER_ROOT=...  EINSIMO_ROOT=...
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT=${GPU_PAPER_ROOT:-$(cd "$SCRIPT_DIR/../../.." && pwd)}
EINSIMO=${EINSIMO_ROOT:-/home/yongchuan/SPHinXsys-EinSIMO-gearbox-paper}
DISC=$ROOT/docs/gpu_paper_discussion
RESOLUTION=${RESOLUTION:-0.0012}
CASE_ID=${CASE_ID:-8}
RUN_SUFFIX=${RUN_SUFFIX:-r01}
ENABLE_CONTACT_FACTOR=${ENABLE_CONTACT_FACTOR:-0}

# Tag like 0.0012 -> 0p0012 for paths
DP_TAG=$(python3 - <<PY
r=float("$RESOLUTION")
s=f"{r:.4f}".rstrip("0").rstrip(".")
print(s.replace(".","p"))
PY
)
RELOAD_STORE=$DISC/results/C1_gearbox/reload/dp${DP_TAG}

# NOTE: oneAPI setvars.sh exports BIN_DIR=bin64 — never use that name for our paths.
case "$DEVICE" in
  gpu_fp32)
    C1_BIN_DIR=$EINSIMO/build/tests/tests_extra/test_3d_gearbox_real_sycl/bin
    C1_BIN=$C1_BIN_DIR/test_3d_gearbox_real_sycl
    NEED_ONEAPI=1
    BACKEND_TAG=gpu_fp32
    ;;
  cpu_fp32)
    C1_BIN_DIR=$EINSIMO/build/tests/tests_extra/test_3d_gearbox_real_cpu/bin
    C1_BIN=$C1_BIN_DIR/test_3d_gearbox_real_cpu
    NEED_ONEAPI=1
    BACKEND_TAG=cpu_fp32
    ;;
  *)
    echo "bad DEVICE=$DEVICE (gpu_fp32|cpu_fp32)"
    exit 2
    ;;
esac

C1_BIN_RELOAD=$C1_BIN_DIR/reload

stage_reload_from_store() {
  mkdir -p "$C1_BIN_RELOAD"
  local missing=0
  for f in Gear1_rld.xml Gear2_rld.xml Gear3_rld.xml Wall_rld.xml; do
    if [[ ! -f "$RELOAD_STORE/$f" ]]; then
      echo "MISSING $RELOAD_STORE/$f — run relax for RESOLUTION=$RESOLUTION first"
      missing=1
    fi
  done
  if [[ "$missing" == "1" ]]; then
    return 1
  fi
  cp -f "$RELOAD_STORE"/Gear1_rld.xml "$RELOAD_STORE"/Gear2_rld.xml \
        "$RELOAD_STORE"/Gear3_rld.xml "$RELOAD_STORE"/Wall_rld.xml \
        "$C1_BIN_RELOAD/"
  echo "Staged reload dp=$RESOLUTION from $RELOAD_STORE -> $C1_BIN_RELOAD"
}

archive_reload_to_store() {
  mkdir -p "$RELOAD_STORE"
  local missing=0
  for f in Gear1_rld.xml Gear2_rld.xml Gear3_rld.xml Wall_rld.xml; do
    if [[ ! -f "$C1_BIN_RELOAD/$f" ]]; then
      echo "MISSING after relax: $C1_BIN_RELOAD/$f"
      missing=1
    fi
  done
  if [[ "$missing" == "1" ]]; then
    return 1
  fi
  cp -f "$C1_BIN_RELOAD"/Gear1_rld.xml "$C1_BIN_RELOAD"/Gear2_rld.xml \
        "$C1_BIN_RELOAD"/Gear3_rld.xml "$C1_BIN_RELOAD"/Wall_rld.xml \
        "$RELOAD_STORE/"
  cat > "$RELOAD_STORE/README.txt" <<EOF
dp=$RESOLUTION relaxed solid reload (Gear1/2/3/Wall)
Produced by: run_c1.sh $DEVICE relax
Oil is Lattice at runtime (not stored here).
Relax time is excluded from GPU-paper solver metrics.
Archived: $(date -Is)
EOF
  echo "Archived reload -> $RELOAD_STORE"
}

case "$MODE" in
  relax)
    # Separate job: lattice solids -> 1000 relax steps -> write reload; then exit.
    END_TIME=${END_TIME:-0}
    ENABLE_BODY=0
    ENABLE_INIT_VTP=0
    ENABLE_REGULAR_VTP=0
    ENABLE_CYCLE_SNAP=0
    DO_RELAX=1
    DO_RELOAD=0
    ;;
  smoke)
    END_TIME=${END_TIME:-0.01}
    ENABLE_BODY=0
    ENABLE_INIT_VTP=0
    ENABLE_REGULAR_VTP=0
    ENABLE_CYCLE_SNAP=0
    DO_RELAX=0
    DO_RELOAD=1
    ;;
  benchmark)
    END_TIME=${END_TIME:-1.05}
    ENABLE_BODY=1
    ENABLE_INIT_VTP=0
    ENABLE_REGULAR_VTP=0
    ENABLE_CYCLE_SNAP=1
    DO_RELAX=0
    DO_RELOAD=1
    ;;
  showcase)
    END_TIME=${END_TIME:-1.05}
    ENABLE_BODY=1
    ENABLE_INIT_VTP=1
    ENABLE_REGULAR_VTP=0
    ENABLE_CYCLE_SNAP=1
    DO_RELAX=0
    DO_RELOAD=1
    ;;
  *)
    echo "bad MODE=$MODE (relax|smoke|benchmark|showcase)"
    exit 2
    ;;
esac

RUN_ID=C1_case${CASE_ID}_${BACKEND_TAG}_dp${DP_TAG}_${MODE}_${RUN_SUFFIX}
RUN=$DISC/results/C1_gearbox/runs/$RUN_ID
LOG=$DISC/logs/gp_${RUN_ID}.log
mkdir -p "$RUN" "$(dirname "$LOG")" "$DISC/results/C1_gearbox/schedules" "$RELOAD_STORE"

if [[ "$NEED_ONEAPI" == "1" ]]; then
  set +u
  # shellcheck disable=SC1091
  source /opt/intel/oneapi/setvars.sh --force >/dev/null 2>&1 || true
  set -eo pipefail
  # setvars overwrites BIN_DIR=bin64; keep our paths on C1_* names only.
fi

if [[ ! -x "$C1_BIN" ]]; then
  echo "missing binary: $C1_BIN" | tee "$LOG"
  echo "EXIT:127" | tee -a "$LOG"
  exit 127
fi

if [[ "$DEVICE" == "gpu_fp32" ]]; then
  export CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-1}
fi

export GPU_PAPER_RUN_ID=$RUN_ID
export GPU_PAPER_MODE=$MODE
export GPU_PAPER_RUN_DIR=$RUN
export GPU_PAPER_END_TIME=$END_TIME
export GPU_PAPER_RESOLUTION=$RESOLUTION

if [[ "$DO_RELOAD" == "1" ]]; then
  stage_reload_from_store | tee -a "$LOG"
fi

DISABLE_VTP=1
DISABLE_RESTART=1
if [[ "$ENABLE_BODY" == "1" && ( "$ENABLE_CYCLE_SNAP" == "1" || "$ENABLE_INIT_VTP" == "1" || "$ENABLE_REGULAR_VTP" == "1" ) ]]; then
  DISABLE_VTP=0
fi

ARGS=(
  --state_recording=false
  --case_id="$CASE_ID"
  --enable_case_preset=1
  --resolution_override="$RESOLUTION"
  --case_run_relaxation="$DO_RELAX"
  --case_reload_particles="$DO_RELOAD"
  --enable_contact_factor="$ENABLE_CONTACT_FACTOR"
  --enable_thermal_model=1
  --enable_body_state_recorder="$ENABLE_BODY"
  --enable_initial_vtp_output="$ENABLE_INIT_VTP"
  --enable_regular_vtp_output="$ENABLE_REGULAR_VTP"
  --disable_vtp_output="$DISABLE_VTP"
  --disable_restart_output="$DISABLE_RESTART"
  --snapshot_times_s=
  --enable_cycle_snapshots="$ENABLE_CYCLE_SNAP"
  --cycle_snapshot_window_rev=1
  --cycle_snapshot_interval_rev=0.25
  --run_tag="$RUN_ID"
)

if [[ "$DO_RELAX" == "1" ]]; then
  ARGS+=(--relax=true --reload=false)
else
  ARGS+=(--relax=false --reload=true --end_time="$END_TIME")
fi

{
  echo "START:$(date -Is)"
  echo "RUN_ID=$RUN_ID DEVICE=$DEVICE MODE=$MODE END_TIME=$END_TIME CASE_ID=$CASE_ID RESOLUTION=$RESOLUTION"
  echo "DO_RELAX=$DO_RELAX DO_RELOAD=$DO_RELOAD enable_contact_factor=$ENABLE_CONTACT_FACTOR"
  echo "RELOAD_STORE=$RELOAD_STORE"
  echo "BIN=$C1_BIN"
  echo "C1_BIN_DIR=$C1_BIN_DIR"
  echo "RUN_DIR=$RUN"
  echo "ARGS=${ARGS[*]}"
  ls -la "$C1_BIN"
  if [[ "$DO_RELOAD" == "1" ]]; then
    ls -la "$C1_BIN_RELOAD"/Gear*_rld.xml "$C1_BIN_RELOAD"/Wall_rld.xml
  fi
} | tee "$LOG"

cd "$C1_BIN_DIR"
set +e
if command -v stdbuf >/dev/null 2>&1; then
  stdbuf -oL -eL "$C1_BIN" "${ARGS[@]}" 2>&1 | tee -a "$LOG"
  rc=${PIPESTATUS[0]}
else
  "$C1_BIN" "${ARGS[@]}" 2>&1 | tee -a "$LOG"
  rc=$?
fi
set -e

if [[ "$DO_RELAX" == "1" && "$rc" == "0" ]]; then
  archive_reload_to_store | tee -a "$LOG"
fi

# Fallback archive of physics outputs
if [[ "$DO_RELAX" != "1" && -d "$C1_BIN_DIR" ]]; then
  out_candidate=$(ls -dt "$C1_BIN_DIR"/output_"${RUN_ID}"* 2>/dev/null | head -1 || true)
  if [[ -n "$out_candidate" && -d "$out_candidate" ]]; then
    for f in run_manifest.json runtime_breakdown.csv runtime_breakdown_progress.csv \
             reducer_thermal_history.csv run_summary.csv stage_timing.csv exit_status.json \
             selected_snapshot_times.txt SingularVariable_ActionOnGear*.dat \
             Observer*_Temperature.dat Oil_ReducedQuantity.dat *.vtp *.vtu; do
      # shellcheck disable=SC2086
      for src in "$out_candidate"/$f; do
        [[ -f "$src" ]] && cp -f "$src" "$RUN/" || true
      done
    done
  fi
fi

echo "EXIT:${rc}" | tee -a "$LOG"
exit "$rc"
