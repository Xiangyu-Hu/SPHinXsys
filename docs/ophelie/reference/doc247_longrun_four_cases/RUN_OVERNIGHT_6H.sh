#!/usr/bin/env bash
# Doc247 overnight: French NC (~3h wall) then French stirring (~3h wall).
# Compile yourself first (see header). Dense VTP + timed french_T_rz_tXXXX.csv snapshots.
set -euo pipefail

########################################
# 0) 你先编译（本脚本不编译）
########################################
# cmake --build /home/ubuntu/SPHinXsys/build --target \
#   test_3d_ophelie_french_natural_convection_frozen_q \
#   test_3d_ophelie_french_stirring_em \
#   -j"$(nproc)"

SRC=/home/ubuntu/SPHinXsys
BUILD=${BUILD:-$SRC/build}
OUT_ROOT=$SRC/docs/ophelie/reference/doc247_longrun_four_cases/overnight_6h
FROZEN=$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_natural_convection_frozen_q/bin
FR_STIR=$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_stirring_em/bin
RELOAD_FR_NC=$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_natural_glass_relax/bin/reload
RELOAD_FR_STIR=$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_stirring_glass_relax/bin/reload
SEED_NC=$SRC/docs/ophelie/reference/standard_data/doc247_seeds/SEED_NC_Fig8_hotmelt_101_ASSUMED_T_rz.csv
SEED_STIR=$SRC/docs/ophelie/reference/standard_data/doc247_seeds/SEED_STIR_Fig8a_axisym_ASSUMED_T_rz.csv
TABLES=$SRC/docs/ophelie/reference/doc239_runs/solver_tables

# Wall budgets (seconds). 3h each; NC also has high end-time capped by timeout.
WALL_NC_S=${WALL_NC_S:-10800}
WALL_STIR_H=${WALL_STIR_H:-3}
# Physical end-time ceilings (wall budget usually hits first for stir).
END_NC=${END_NC:-600}
END_STIR=${END_STIR:-450}

export ONEAPI_DEVICE_SELECTOR=${ONEAPI_DEVICE_SELECTOR:-cuda:gpu}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-8}
export TBB_NUM_THREADS=${TBB_NUM_THREADS:-8}

mkdir -p "$OUT_ROOT"
echo "[doc247] overnight start $(date -Is) OUT=$OUT_ROOT" | tee -a "$OUT_ROOT/overnight.log"

need() { [[ -x "$1" ]] || { echo "FATAL: missing binary $1 — compile first"; exit 1; }; }
need "$FROZEN/test_3d_ophelie_french_natural_convection_frozen_q"
need "$FR_STIR/test_3d_ophelie_french_stirring_em"
[[ -f "$SEED_NC" ]] || { echo "FATAL: missing $SEED_NC"; exit 1; }
[[ -f "$SEED_STIR" ]] || { echo "FATAL: missing $SEED_STIR"; exit 1; }
[[ -f "$RELOAD_FR_NC/Reload.xml" ]] || { echo "FATAL: missing NC reload"; exit 1; }
[[ -f "$RELOAD_FR_STIR/Reload.xml" ]] || { echo "FATAL: missing STIR reload"; exit 1; }

########################################
# 1) French NC — Fig8 hot-melt IC + C0 + real EM
#    Estimate: ~0.05–0.12 s_phys / s_wall → 3h ≈ 500–1300 s_phys
#    Cap end-time=600; timeout ${WALL_NC_S}s
########################################
run_nc() {
  local OUT=$OUT_ROOT/01_FR_NC_Fig8seed_EM
  mkdir -p "$OUT/out" "$OUT/trz_snaps"
  cd "$FROZEN"
  rm -rf output output_backup && mkdir -p output
  echo "[doc247][NC] start $(date -Is) end-time=$END_NC wall=${WALL_NC_S}s" | tee -a "$OUT_ROOT/overnight.log"

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
  local rc=${PIPESTATUS[0]}
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
}

########################################
# 2) French stirring — Fig8a axisym IC + 60kW/10rpm + skull/crust
#    Estimate: ~0.039 s_phys / s_wall → 3h ≈ 420 s_phys
#    end-time=450; max-wall-hours=3
########################################
run_stir() {
  local OUT=$OUT_ROOT/02_FR_STIR_Fig8a_EM
  mkdir -p "$OUT/out" "$OUT/trz_snaps"
  cd "$FR_STIR"
  rm -rf output && mkdir -p output
  echo "[doc247][STIR] start $(date -Is) end-time=$END_STIR max-wall-hours=$WALL_STIR_H" \
    | tee -a "$OUT_ROOT/overnight.log"

  set +e
  ./test_3d_ophelie_french_stirring_em \
    --reload-dir="$RELOAD_FR_STIR" \
    --init-T-rz-csv="$SEED_STIR" \
    --rotation-rpm=10 \
    --thermo-em-coupling=periodic --em-update-interval=0.5 \
    --em-control=fixed-current --aind=off \
    --skull-resistance --skull-delta=0.005 \
    --free-crust-resistance --free-crust-delta=0.005 --free-crust-k=1 \
    --end-time="$END_STIR" --max-wall-hours="$WALL_STIR_H" \
    --state-record-interval=5 \
    --screen-every=10 --csv-every=10 \
    2>&1 | tee "$OUT/run.log"
  local rc=${PIPESTATUS[0]}
  set -e
  echo "[doc247][STIR] exit=$rc $(date -Is)" | tee -a "$OUT_ROOT/overnight.log"

  cp -f output/*.csv "$OUT/" 2>/dev/null || true
  cp -f output/french_T_rz_t*.csv "$OUT/trz_snaps/" 2>/dev/null || true
  cp -a output/*.vtp "$OUT/out/" 2>/dev/null || true
  echo "$rc" >"$OUT/exit_code.txt"
}

run_nc
run_stir
echo "[doc247] overnight done $(date -Is)" | tee -a "$OUT_ROOT/overnight.log"
