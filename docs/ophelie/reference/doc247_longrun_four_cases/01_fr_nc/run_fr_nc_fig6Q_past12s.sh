#!/usr/bin/env bash
# Trial: push French NC Fig6Q past the ~12 s blowup seen with Fig8 seed.
#
# Modes (SEED_MODE):
#   fig8      — same Fig8 hot-melt IC as short test (default); adds softer dt/CFL
#   coldwall  — Doc246 C0-stable coldwall skin IC (control: can hydro pass 12 s?)
#
# Soft knobs (override via env):
#   DT_SCALE=0.5  ACOUSTIC_CFL=0.35  END_NC=20  WALL_NC_S=1800
set -euo pipefail

SRC=/home/ubuntu/SPHinXsys
BUILD=${BUILD:-$SRC/build}
OUT_ROOT=${OUT_ROOT:-$SRC/docs/ophelie/reference/doc247_longrun_four_cases/overnight_6h}
FROZEN=$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_natural_convection_frozen_q/bin
RELOAD_FR_NC=$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_natural_glass_relax/bin/reload
SEED_Q=$SRC/docs/ophelie/reference/standard_data/doc247_seeds/SEED_NC_Fig6_Qshape_ASSUMED_rz.csv
TABLES=$SRC/docs/ophelie/reference/doc239_runs/solver_tables

SEED_MODE=${SEED_MODE:-fig8}
case "$SEED_MODE" in
  fig8)
    SEED_T=$SRC/docs/ophelie/reference/standard_data/doc247_seeds/SEED_NC_Fig8_hotmelt_101_ASSUMED_T_rz.csv
    TAG=Fig8T_softDT
    ;;
  coldwall)
    SEED_T=$SRC/docs/ophelie/reference/doc243_runs/03_thermal_initial_guess/SEED_coldwall_skin_ASSUMED_T_rz.csv
    TAG=coldwall_softDT
    ;;
  *)
    echo "SEED_MODE=fig8|coldwall"; exit 2
    ;;
esac

WALL_NC_S=${WALL_NC_S:-1800}
END_NC=${END_NC:-20}
P_ABS=${P_ABS:-50000}
DT_SCALE=${DT_SCALE:-0.5}
ACOUSTIC_CFL=${ACOUSTIC_CFL:-0.35}
ADVECTION_CFL=${ADVECTION_CFL:-0.2}

export ONEAPI_DEVICE_SELECTOR=${ONEAPI_DEVICE_SELECTOR:-cuda:gpu}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-8}
export TBB_NUM_THREADS=${TBB_NUM_THREADS:-8}

BIN=$FROZEN/test_3d_ophelie_french_natural_convection_frozen_q
[[ -x "$BIN" ]] || { echo "FATAL: missing $BIN"; exit 1; }
[[ -f "$SEED_T" ]] || { echo "FATAL: missing $SEED_T"; exit 1; }
[[ -f "$SEED_Q" ]] || { echo "FATAL: missing $SEED_Q"; exit 1; }
[[ -f "$RELOAD_FR_NC/Reload.xml" ]] || { echo "FATAL: missing NC reload"; exit 1; }

OUT=$OUT_ROOT/01_FR_NC_Fig6Q_past12_${TAG}_end${END_NC}
mkdir -p "$OUT/out" "$OUT/trz_snaps" "$OUT_ROOT"
cd "$FROZEN"
rm -rf output output_backup && mkdir -p output

echo "[doc247][NC-past12] start $(date -Is) mode=$SEED_MODE end=$END_NC dt_scale=$DT_SCALE acfl=$ACOUSTIC_CFL OUT=$OUT" \
  | tee -a "$OUT_ROOT/overnight.log"

set +e
timeout --signal=TERM --kill-after=120 "${WALL_NC_S}s" \
  ./test_3d_ophelie_french_natural_convection_frozen_q \
    --reload-dir="$RELOAD_FR_NC" \
    --dp=0.015 \
    --doc223-paper-flat-slip-top \
    --prescribed-heat-power="$P_ABS" \
    --init-Q-rz-csv="$SEED_Q" \
    --enable-thermal-diffusion \
    --temperature-dependent-properties \
    --doc226-local-particle-viscosity \
    --viscosity-mode=implicit_split \
    --material-dataset=JAC_THESIS_RAW_BORROWED \
    --material-tables-dir="$TABLES" \
    --material-evidence-dir="$OUT" \
    --init-T-rz-csv="$SEED_T" \
    --skull-resistance \
    --skull-delta-side=0.005 --skull-delta-bottom=0.005 \
    --skull-k=1 \
    --skull-h-cool-side=300 --skull-h-cool-bottom=35 \
    --free-crust-resistance --free-crust-delta=0.005 --free-crust-k=1 \
    --doc224-u-blowup-gate=50 \
    --dt-scale="$DT_SCALE" \
    --acoustic-cfl="$ACOUSTIC_CFL" \
    --advection-cfl="$ADVECTION_CFL" \
    --end-time="$END_NC" --hydro-settle-time=0.05 \
    --state-record-interval=1.0 \
    --meridional-average-start=0.5 \
    --boundary-face-t-mode=local_dt_dz \
    --summary-json="$OUT/summary.json" \
    2>&1 | tee "$OUT/run.log"
rc=${PIPESTATUS[0]}
set -e

echo "[doc247][NC-past12] exit=$rc $(date -Is) mode=$SEED_MODE" | tee -a "$OUT_ROOT/overnight.log"
cp -f output/french_energy_budget.csv \
      output/french_meridional_timeseries.csv \
      output/french_spatial_stats.csv \
      output/french_T_rz.csv \
      output/french_q_rz.csv \
      "$OUT/" 2>/dev/null || true
cp -f output/french_T_rz_t*.csv "$OUT/trz_snaps/" 2>/dev/null || true
cp -a output/*.vtp "$OUT/out/" 2>/dev/null || true
echo "$rc" >"$OUT/exit_code.txt"

# quick verdict line
if rg -q 'numerical-abort' "$OUT/run.log"; then
  rg -n 'numerical-abort|U_max_above' "$OUT/run.log" | tail -3 | tee -a "$OUT_ROOT/overnight.log"
elif rg -q 'completed_requested_time=1' "$OUT/run.log"; then
  echo "[doc247][NC-past12] OK completed end-time=$END_NC" | tee -a "$OUT_ROOT/overnight.log"
fi
exit "$rc"
