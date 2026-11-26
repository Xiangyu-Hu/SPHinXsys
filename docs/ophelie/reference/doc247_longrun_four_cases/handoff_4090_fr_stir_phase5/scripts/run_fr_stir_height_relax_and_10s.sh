#!/usr/bin/env bash
# French stir — impeller immersion (z) sweep: full relax (glass+rotor+wall) per height, then 10 s EM.
#
# Baseline paddle z≈[0.04,0.21] m in glass z∈[0,0.21].
#   HZ_dn02: Δz=-0.02  (deeper / 更深入)
#   HZ_z00 : Δz= 0     (baseline; default REUSE stock reload unless FORCE_RELAX_BASELINE=1)
#   HZ_up02: Δz=+0.02  (raised / 抬高；桨顶略出自由面，ASSUMED)
#
# Env:
#   CASES="HZ_dn02 HZ_z00 HZ_up02"
#   SKIP_RELAX=1          # only run EM if reload_* already exist
#   SKIP_EXISTING=1       # skip EM case if exit_code.txt==0
#   FORCE_RELAX_BASELINE=1
#   END_PT=10 STATE_DT=5 WALL_H=1
#   ONEAPI_DEVICE_SELECTOR=cuda:1
set -euo pipefail

_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
_REPO_FROM_SCRIPT="$(cd "$_SCRIPT_DIR/../../../../.." && pwd)"
SRC=${SRC:-${_REPO_FROM_SCRIPT:-/home/ubuntu/SPHinXsys}}
BUILD=${BUILD:-$SRC/build}
OUT_ROOT=${OUT_ROOT:-$SRC/docs/ophelie/reference/doc247_longrun_four_cases/overnight_6h/phase5_fr_stir_height10s}

RELAX_BIN_DIR=$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_stirring_glass_relax/bin
RELAX_BIN=$RELAX_BIN_DIR/test_3d_ophelie_french_stirring_glass_relax
STOCK_RELOAD=$RELAX_BIN_DIR/reload
FR_STIR=$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_stirring_em/bin
EM_BIN=$FR_STIR/test_3d_ophelie_french_stirring_em

HO=$SRC/docs/ophelie/reference/doc247_longrun_four_cases/handoff_4090_fr_stir_phase5
SEED_T=$HO/seeds/SEED_STIR_Fig8a_axisym_ASSUMED_T_rz.csv
SEED_Q=$HO/seeds/SEED_STIR_Fig10_Qshape_ASSUMED_rz.csv

END_PT=${END_PT:-10}
STATE_DT=${STATE_DT:-5}
WALL_H=${WALL_H:-1}
BLEND=${BLEND:-0.8}
EVERY=${EVERY:-1}
RPM=${RPM:-10}
PABS=${PABS:-60000}
SKULL=${SKULL:-0.005}
SKIP_RELAX=${SKIP_RELAX:-0}
SKIP_EXISTING=${SKIP_EXISTING:-0}
FORCE_RELAX_BASELINE=${FORCE_RELAX_BASELINE:-0}
CASES=${CASES:-"HZ_dn02 HZ_z00 HZ_up02"}

export ONEAPI_DEVICE_SELECTOR=${ONEAPI_DEVICE_SELECTOR:-cuda:gpu}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-8}
export TBB_NUM_THREADS=${TBB_NUM_THREADS:-8}

[[ -x "$RELAX_BIN" ]] || { echo "FATAL: missing $RELAX_BIN — compile test_3d_ophelie_french_stirring_glass_relax"; exit 1; }
[[ -x "$EM_BIN" ]] || { echo "FATAL: missing $EM_BIN — compile test_3d_ophelie_french_stirring_em"; exit 1; }
[[ -f "$SEED_T" && -f "$SEED_Q" ]] || { echo "FATAL: missing seeds"; exit 1; }

mkdir -p "$OUT_ROOT"
echo "[height10s] start $(date -Is) CASES=$CASES END=$END_PT OUT=$OUT_ROOT" | tee -a "$OUT_ROOT/run_matrix.log"

# case_id → dz
case_dz() {
  case "$1" in
    HZ_dn02) echo -0.02 ;;
    HZ_z00)  echo 0.0 ;;
    HZ_up02) echo 0.02 ;;
    *) echo "FATAL: unknown case $1" >&2; return 1 ;;
  esac
}

reload_dir_for() {
  echo "$RELAX_BIN_DIR/reload_$1"
}

relax_one() {
  local id=$1
  local dz
  dz=$(case_dz "$id")
  local rdir
  rdir=$(reload_dir_for "$id")

  if [[ "$SKIP_RELAX" == "1" && -f "$rdir/Reload.xml" ]]; then
    echo "[height10s] skip relax $id (found $rdir/Reload.xml)" | tee -a "$OUT_ROOT/run_matrix.log"
    return 0
  fi

  # Baseline: reuse stock reload unless forced (same Δz=0 geometry).
  if [[ "$id" == "HZ_z00" && "$FORCE_RELAX_BASELINE" != "1" && -f "$STOCK_RELOAD/Reload.xml" ]]; then
    mkdir -p "$rdir"
    # copy stock into named dir (do not move — keep stock intact)
    rm -rf "$rdir"
    cp -a "$STOCK_RELOAD" "$rdir"
    echo "[height10s] HZ_z00 reuse stock reload → $rdir" | tee -a "$OUT_ROOT/run_matrix.log"
    return 0
  fi

  echo "[height10s] relax $id dz=$dz start $(date -Is)" | tee -a "$OUT_ROOT/run_matrix.log"
  cd "$RELAX_BIN_DIR"
  rm -rf reload output
  mkdir -p output
  # Full bodies: glass cavity must be re-cut for translated paddle.
  ./test_3d_ophelie_french_stirring_glass_relax \
    --bodies=all \
    --rotor-translation="0,0,$dz" \
    --no-auto-center-rotor \
    2>&1 | tee "$OUT_ROOT/relax_${id}.log"
  [[ -f reload/Reload.xml ]] || { echo "FATAL: relax $id produced no Reload.xml"; return 1; }
  rm -rf "$rdir"
  mv reload "$rdir"
  mkdir -p reload  # leave empty placeholder
  echo "[height10s] relax $id done → $rdir $(date -Is)" | tee -a "$OUT_ROOT/run_matrix.log"
}

run_em_one() {
  local id=$1
  local dz
  dz=$(case_dz "$id")
  local rdir
  rdir=$(reload_dir_for "$id")
  local OUT=$OUT_ROOT/$id
  mkdir -p "$OUT/out" "$OUT/trz_snaps"

  if [[ "$SKIP_EXISTING" == "1" && -f "$OUT/exit_code.txt" && "$(cat "$OUT/exit_code.txt")" == "0" ]]; then
    echo "[height10s] skip EM $id (exit_code=0)" | tee -a "$OUT_ROOT/run_matrix.log"
    return 0
  fi
  [[ -f "$rdir/Reload.xml" ]] || { echo "FATAL: missing $rdir/Reload.xml — relax first"; return 1; }

  echo "[height10s] EM $id rpm=$RPM dz=$dz start $(date -Is)" | tee -a "$OUT_ROOT/run_matrix.log"
  cd "$FR_STIR"
  rm -rf output && mkdir -p output
  set +e
  ./test_3d_ophelie_french_stirring_em \
    --reload-dir="$rdir" \
    --rotor-translation="0,0,$dz" \
    --no-auto-center-rotor \
    --init-T-rz-csv="$SEED_T" \
    --init-Q-rz-csv="$SEED_Q" \
    --fig10-q-correct \
    --fig10-q-correct-blend="$BLEND" \
    --fig10-q-correct-every="$EVERY" \
    --rotation-rpm="$RPM" \
    --em-control=fixed-power --target-absorbed-power="$PABS" \
    --thermo-em-coupling=periodic --em-update-interval=0.5 \
    --aind=off \
    --skull-resistance --skull-delta="$SKULL" \
    --free-crust-resistance --free-crust-delta="$SKULL" --free-crust-k=1 \
    --end-time="$END_PT" --max-wall-hours="$WALL_H" \
    --state-record-interval="$STATE_DT" \
    --screen-every=20 --csv-every=10 \
    2>&1 | tee "$OUT/run.log"
  local rc=${PIPESTATUS[0]}
  set -e

  echo "[height10s] EM $id exit=$rc $(date -Is)" | tee -a "$OUT_ROOT/run_matrix.log"
  cp -f output/*.csv "$OUT/" 2>/dev/null || true
  cp -f output/french_T_rz_t*.csv "$OUT/trz_snaps/" 2>/dev/null || true
  cp -a output/*.vtp "$OUT/out/" 2>/dev/null || true
  echo "$rc" >"$OUT/exit_code.txt"
  echo "$RPM" >"$OUT/rpm.txt"
  echo "$dz" >"$OUT/rotor_dz.txt"
  echo "$PABS" >"$OUT/p_abs.txt"
  printf '%s\n' "case_id=$id rpm=$RPM rotor_dz=$dz P_abs=$PABS skull=$SKULL blend=$BLEND end=$END_PT reload=$rdir" \
    >"$OUT/case_meta.txt"
}

for id in $CASES; do
  relax_one "$id"
  run_em_one "$id"
done

echo "[height10s] finished $(date -Is)" | tee -a "$OUT_ROOT/run_matrix.log"
