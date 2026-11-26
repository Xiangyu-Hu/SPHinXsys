#!/usr/bin/env bash
# French NC — Fig8 hot-melt IC + C0 + real EM (Doc247 locked).
# Defaults: ~3h wall / end-time=600. Override: WALL_NC_S=... END_NC=...
set -euo pipefail

SRC=/home/ubuntu/SPHinXsys
BUILD=${BUILD:-$SRC/build}
OUT_ROOT=${OUT_ROOT:-$SRC/docs/ophelie/reference/doc247_longrun_four_cases/overnight_6h}
FROZEN=$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_natural_convection_frozen_q/bin
RELOAD_FR_NC=$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_natural_glass_relax/bin/reload
SEED_NC=$SRC/docs/ophelie/reference/standard_data/doc247_seeds/SEED_NC_Fig8_hotmelt_101_ASSUMED_T_rz.csv
TABLES=$SRC/docs/ophelie/reference/doc239_runs/solver_tables

WALL_NC_S=${WALL_NC_S:-10800}
END_NC=${END_NC:-600}

export ONEAPI_DEVICE_SELECTOR=${ONEAPI_DEVICE_SELECTOR:-cuda:gpu}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-8}
export TBB_NUM_THREADS=${TBB_NUM_THREADS:-8}

BIN=$FROZEN/test_3d_ophelie_french_natural_convection_frozen_q
[[ -x "$BIN" ]] || { echo "FATAL: missing $BIN — compile first"; exit 1; }
[[ -f "$SEED_NC" ]] || { echo "FATAL: missing $SEED_NC"; exit 1; }
[[ -f "$RELOAD_FR_NC/Reload.xml" ]] || { echo "FATAL: missing NC reload"; exit 1; }

OUT=$OUT_ROOT/01_FR_NC_Fig8seed_EM
mkdir -p "$OUT/out" "$OUT/trz_snaps" "$OUT_ROOT"
cd "$FROZEN"
rm -rf output output_backup && mkdir -p output

echo "[doc247][NC] start $(date -Is) end-time=$END_NC wall=${WALL_NC_S}s OUT=$OUT" \
  | tee -a "$OUT_ROOT/overnight.log"

set +e
timeout --signal=TERM --kill-after=120 "${WALL_NC_S}s" \
  ./test_3d_ophelie_french_natural_convection_frozen_q \
    --reload-dir="$RELOAD_FR_NC" \
    --dp=0.015 \
    --doc223-paper-flat-slip-top \
    --em-control=fixed-power --target-absorbed-power=50000 \
    --thermo-em-coupling=periodic --em-update-interval=0.5 --aind=off \
    --enable-thermal-diffusion \
    --temperature-dependent-properties \
    --doc226-local-particle-viscosity \
    --viscosity-mode=implicit_split \
    --material-dataset=JAC_THESIS_RAW_BORROWED \
    --material-tables-dir="$TABLES" \
    --material-evidence-dir="$OUT" \
    --init-T-rz-csv="$SEED_NC" \
    --skull-resistance \
    --skull-delta-side=0.005 --skull-delta-bottom=0.005 \
    --skull-k=1 \
    --skull-h-cool-side=300 --skull-h-cool-bottom=35 \
    --free-crust-resistance --free-crust-delta=0.005 --free-crust-k=1 \
    --doc224-u-blowup-gate=50 \
    --end-time="$END_NC" --hydro-settle-time=0.05 \
    --state-record-interval=5.0 \
    --meridional-average-start=0.5 \
    --boundary-face-t-mode=local_dt_dz \
    --summary-json="$OUT/summary.json" \
    2>&1 | tee "$OUT/run.log"
rc=${PIPESTATUS[0]}
set -e

echo "[doc247][NC] exit=$rc $(date -Is)" | tee -a "$OUT_ROOT/overnight.log"
cp -f output/french_energy_budget.csv \
      output/french_thermo_em_coupling.csv \
      output/french_meridional_timeseries.csv \
      output/french_meridional_ur_uz_mean.csv \
      output/french_spatial_stats.csv \
      output/french_T_rz.csv \
      "$OUT/" 2>/dev/null || true
cp -f output/french_T_rz_t*.csv "$OUT/trz_snaps/" 2>/dev/null || true
cp -a output/*.vtp "$OUT/out/" 2>/dev/null || true
echo "$rc" >"$OUT/exit_code.txt"
exit "$rc"
