#!/usr/bin/env bash
# Phase5 French stir matrix: 9 cases × 60 s × state every 5 s.
# Fig8a T IC · EM pure · Q←Fig10-correct (PLACEHOLDER seed until true Fig.10 digitization).
# NOT literature PASS. See PLAN_PHASE5_FR_STIR_MATRIX_60s_20261126.md
#
# Env overrides:
#   END_PT=60 STATE_DT=5 WALL_H=2.5 BLEND=0.8 EVERY=1
#   CASES="SW_rpm10 EM_P60"   # subset
#   SKIP_EXISTING=1           # skip case dirs that already have exit_code.txt==0
#   ACOUSTIC_CFL_RPM5=0.35    # soft fix for 5 rpm startup
set -euo pipefail

# Prefer env SRC; else infer repo root from this script; else 4090 default path.
_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
_REPO_FROM_SCRIPT="$(cd "$_SCRIPT_DIR/../../../../.." && pwd)"
SRC=${SRC:-${_REPO_FROM_SCRIPT:-/home/ubuntu/SPHinXsys}}
BUILD=${BUILD:-$SRC/build}
OUT_ROOT=${OUT_ROOT:-$SRC/docs/ophelie/reference/doc247_longrun_four_cases/overnight_6h/phase5_fr_stir_matrix_60s}
FR_STIR=$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_stirring_em/bin
RELOAD_FR_STIR=$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_stirring_glass_relax/bin/reload
HO=$SRC/docs/ophelie/reference/doc247_longrun_four_cases/handoff_4090_fr_stir_phase5
SEED_T=$HO/seeds/SEED_STIR_Fig8a_axisym_ASSUMED_T_rz.csv
SEED_Q=$HO/seeds/SEED_STIR_Fig10_Qshape_ASSUMED_rz.csv
POST=$SRC/docs/ophelie/reference/doc247_longrun_four_cases/02_fr_stir/post_phase5_matrix_compare.py

END_PT=${END_PT:-60}
STATE_DT=${STATE_DT:-5}
WALL_H=${WALL_H:-2.5}
BLEND=${BLEND:-0.8}
EVERY=${EVERY:-1}
SKIP_EXISTING=${SKIP_EXISTING:-0}
ACOUSTIC_CFL_RPM5=${ACOUSTIC_CFL_RPM5:-0.35}
CASES=${CASES:-"SW_rpm05 SW_rpm10 SW_rpm20 EM_P56 EM_P60 EM_P70 ST_d3 ST_d5 ST_d8"}

export ONEAPI_DEVICE_SELECTOR=${ONEAPI_DEVICE_SELECTOR:-cuda:gpu}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-8}
export TBB_NUM_THREADS=${TBB_NUM_THREADS:-8}

BIN=$FR_STIR/test_3d_ophelie_french_stirring_em
[[ -x "$BIN" ]] || { echo "FATAL: missing $BIN — compile first"; exit 1; }
[[ -f "$SEED_T" && -f "$SEED_Q" ]] || { echo "FATAL: missing seeds"; exit 1; }
[[ -f "$RELOAD_FR_STIR/Reload.xml" ]] || { echo "FATAL: missing STIR reload"; exit 1; }

mkdir -p "$OUT_ROOT"
echo "[phase5] start $(date -Is) END=$END_PT STATE_DT=$STATE_DT BLEND=$BLEND CASES=$CASES OUT=$OUT_ROOT" \
  | tee -a "$OUT_ROOT/run_matrix.log"

# case_id → rpm P_abs skull_delta
case_params() {
  local id=$1
  case "$id" in
    SW_rpm05) echo 5 60000 0.005 ;;
    SW_rpm10) echo 10 60000 0.005 ;;
    SW_rpm20) echo 20 60000 0.005 ;;
    EM_P56)   echo 10 56000 0.005 ;;
    EM_P60)   echo 10 60000 0.005 ;;
    EM_P70)   echo 10 70000 0.005 ;;
    ST_d3)    echo 10 60000 0.003 ;;
    ST_d5)    echo 10 60000 0.005 ;;
    ST_d8)    echo 10 60000 0.008 ;;
    *) echo "FATAL: unknown case_id=$id" >&2; return 1 ;;
  esac
}

run_one() {
  local id=$1
  read -r rpm pabs skull < <(case_params "$id")
  local OUT=$OUT_ROOT/$id
  mkdir -p "$OUT/out" "$OUT/trz_snaps"
  if [[ "$SKIP_EXISTING" == "1" && -f "$OUT/exit_code.txt" && "$(cat "$OUT/exit_code.txt")" == "0" ]]; then
    echo "[phase5] skip $id (exit_code=0)" | tee -a "$OUT_ROOT/run_matrix.log"
    return 0
  fi

  local extra=()
  local rpm_i=${rpm%.*}
  if [[ "$rpm_i" -le 5 ]]; then
    extra+=(--acoustic-cfl="$ACOUSTIC_CFL_RPM5")
  fi

  echo "[phase5] $id rpm=$rpm P=$pabs skull=$skull start $(date -Is) extra=${extra[*]:-}" \
    | tee -a "$OUT_ROOT/run_matrix.log"

  cd "$FR_STIR"
  rm -rf output && mkdir -p output
  set +e
  ./test_3d_ophelie_french_stirring_em \
    --reload-dir="$RELOAD_FR_STIR" \
    --init-T-rz-csv="$SEED_T" \
    --init-Q-rz-csv="$SEED_Q" \
    --fig10-q-correct \
    --fig10-q-correct-blend="$BLEND" \
    --fig10-q-correct-every="$EVERY" \
    --rotation-rpm="$rpm" \
    --em-control=fixed-power --target-absorbed-power="$pabs" \
    --thermo-em-coupling=periodic --em-update-interval=0.5 \
    --aind=off \
    --skull-resistance --skull-delta="$skull" \
    --free-crust-resistance --free-crust-delta="$skull" --free-crust-k=1 \
    --end-time="$END_PT" --max-wall-hours="$WALL_H" \
    --state-record-interval="$STATE_DT" \
    --screen-every=20 --csv-every=10 \
    "${extra[@]}" \
    2>&1 | tee "$OUT/run.log"
  local rc=${PIPESTATUS[0]}
  set -e

  echo "[phase5] $id exit=$rc $(date -Is)" | tee -a "$OUT_ROOT/run_matrix.log"
  cp -f output/*.csv "$OUT/" 2>/dev/null || true
  cp -f output/french_T_rz_t*.csv "$OUT/trz_snaps/" 2>/dev/null || true
  cp -a output/*.vtp "$OUT/out/" 2>/dev/null || true
  echo "$rc" >"$OUT/exit_code.txt"
  echo "$rpm" >"$OUT/rpm.txt"
  echo "$pabs" >"$OUT/p_abs.txt"
  echo "$skull" >"$OUT/skull_delta.txt"
  printf '%s\n' "case_id=$id rpm=$rpm P_abs=$pabs skull_delta=$skull blend=$BLEND every=$EVERY end=$END_PT" \
    >"$OUT/case_meta.txt"
  return 0
}

for id in $CASES; do
  run_one "$id" || true
done

echo "[phase5] all cases done $(date -Is)" | tee -a "$OUT_ROOT/run_matrix.log"
if [[ -f "$POST" ]]; then
  python3 "$POST" --root "$OUT_ROOT" --seed-t "$SEED_T" --seed-q "$SEED_Q" \
    --times 20,40,60 || true
fi
echo "[phase5] finished $(date -Is)" | tee -a "$OUT_ROOT/run_matrix.log"
