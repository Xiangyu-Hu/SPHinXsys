#!/usr/bin/env bash
# French NC: live EM → blend toward Fig6 Q shape → next EM (ASSUMED engineering mix).
# NOT pure prescribed Fig6. Formal default: blend=0.8, every=1 (correct after each EM).
#
# Env:
#   END_NC=60 WALL_NC_S=7200 BLEND=0.8 EVERY=1 P_ABS=50000
#   DT_SCALE=0.5 ACOUSTIC_CFL=0.35  (optional soft step; helped past ~12 s on prescribed path)
#   EVERY=2 → correct only on EM updates 2,4,6...
set -euo pipefail

SRC=/home/ubuntu/SPHinXsys
BUILD=${BUILD:-$SRC/build}
OUT_ROOT=${OUT_ROOT:-$SRC/docs/ophelie/reference/doc247_longrun_four_cases/overnight_6h}
FROZEN=$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_natural_convection_frozen_q/bin
RELOAD_FR_NC=$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_natural_glass_relax/bin/reload
SEED_T=$SRC/docs/ophelie/reference/standard_data/doc247_seeds/SEED_NC_Fig8_hotmelt_101_ASSUMED_T_rz.csv
SEED_Q=$SRC/docs/ophelie/reference/standard_data/doc247_seeds/SEED_NC_Fig6_Qshape_ASSUMED_rz.csv
TABLES=$SRC/docs/ophelie/reference/doc239_runs/solver_tables

WALL_NC_S=${WALL_NC_S:-7200}
END_NC=${END_NC:-60}
P_ABS=${P_ABS:-50000}
BLEND=${BLEND:-0.8}
EVERY=${EVERY:-1}
DT_SCALE=${DT_SCALE:-}
ACOUSTIC_CFL=${ACOUSTIC_CFL:-}
ADVECTION_CFL=${ADVECTION_CFL:-}

export ONEAPI_DEVICE_SELECTOR=${ONEAPI_DEVICE_SELECTOR:-cuda:gpu}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-8}
export TBB_NUM_THREADS=${TBB_NUM_THREADS:-8}

BIN=$FROZEN/test_3d_ophelie_french_natural_convection_frozen_q
[[ -x "$BIN" ]] || { echo "FATAL: missing $BIN — recompile frozen_q"; exit 1; }
[[ -f "$SEED_T" && -f "$SEED_Q" ]] || { echo "FATAL: missing seeds"; exit 1; }
[[ -f "$RELOAD_FR_NC/Reload.xml" ]] || { echo "FATAL: missing NC reload"; exit 1; }

SOFT_TAG=""
EXTRA_DT=()
if [[ -n "${DT_SCALE}" ]]; then
  EXTRA_DT+=(--dt-scale="$DT_SCALE")
  SOFT_TAG="${SOFT_TAG}_dts${DT_SCALE}"
fi
if [[ -n "${ACOUSTIC_CFL}" ]]; then
  EXTRA_DT+=(--acoustic-cfl="$ACOUSTIC_CFL")
  SOFT_TAG="${SOFT_TAG}_acfl${ACOUSTIC_CFL}"
fi
if [[ -n "${ADVECTION_CFL}" ]]; then
  EXTRA_DT+=(--advection-cfl="$ADVECTION_CFL")
fi

OUT=$OUT_ROOT/01_FR_NC_EM_Fig6correct_b${BLEND}_e${EVERY}_end${END_NC}${SOFT_TAG}
mkdir -p "$OUT/out" "$OUT/trz_snaps" "$OUT_ROOT"
cd "$FROZEN"
rm -rf output output_backup && mkdir -p output

echo "[doc247][EM-Fig6correct] start $(date -Is) blend=$BLEND every=$EVERY end=$END_NC dt_scale=${DT_SCALE:-1} OUT=$OUT" \
  | tee -a "$OUT_ROOT/overnight.log"

set +e
timeout --signal=TERM --kill-after=120 "${WALL_NC_S}s" \
  ./test_3d_ophelie_french_natural_convection_frozen_q \
    --reload-dir="$RELOAD_FR_NC" \
    --dp=0.015 \
    --doc223-paper-flat-slip-top \
    --em-control=fixed-power --target-absorbed-power="$P_ABS" \
    --thermo-em-coupling=periodic --em-update-interval=0.5 --aind=off \
    --fig6-q-correct \
    --fig6-q-correct-blend="$BLEND" \
    --fig6-q-correct-every="$EVERY" \
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
    "${EXTRA_DT[@]}" \
    --end-time="$END_NC" --hydro-settle-time=0.05 \
    --state-record-interval=5.0 \
    --meridional-average-start=0.5 \
    --boundary-face-t-mode=local_dt_dz \
    --summary-json="$OUT/summary.json" \
    2>&1 | tee "$OUT/run.log"
rc=${PIPESTATUS[0]}
set -e

echo "[doc247][EM-Fig6correct] exit=$rc $(date -Is)" | tee -a "$OUT_ROOT/overnight.log"
cp -f output/french_energy_budget.csv \
      output/french_thermo_em_coupling.csv \
      output/french_spatial_stats.csv \
      output/french_T_rz.csv \
      output/french_q_rz.csv \
      "$OUT/" 2>/dev/null || true
cp -f output/french_T_rz_t*.csv "$OUT/trz_snaps/" 2>/dev/null || true
cp -a output/*.vtp "$OUT/out/" 2>/dev/null || true
echo "$rc" >"$OUT/exit_code.txt"
rg -n 'Fig6Q-correct|thermo-em\] update|numerical-abort|completed_requested' "$OUT/run.log" | head -40 || true
exit "$rc"
