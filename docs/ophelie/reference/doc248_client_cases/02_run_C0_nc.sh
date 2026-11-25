#!/usr/bin/env bash
# Doc248 C0 — client NC (flat slip + ASSUMED 15 kW EM). Requires 01_generate_particles first.
set -euo pipefail
SRC=${SRC:-/home/ubuntu/SPHinXsys}
BUILD=${BUILD:-$SRC/build}
PACK=$SRC/docs/ophelie/reference/doc248_client_cases
OUT=${OUT:-$PACK/runs/C0_nc}
FROZEN=$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_natural_convection_frozen_q/bin
RELOAD=${RELOAD:-$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_client_nc_glass_relax/bin/reload}
END_TIME=${END_TIME:-60}
WALL_S=${WALL_S:-36000}   # 10 h wall cap via timeout

export ONEAPI_DEVICE_SELECTOR=${ONEAPI_DEVICE_SELECTOR:-cuda:gpu}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-8}
export TBB_NUM_THREADS=${TBB_NUM_THREADS:-8}

[[ -f "$RELOAD/Reload.xml" ]] || { echo "FATAL: no Reload at $RELOAD — run 01_generate_particles.sh"; exit 1; }
[[ -x "$FROZEN/test_3d_ophelie_french_natural_convection_frozen_q" ]] || { echo "FATAL: compile frozen_q"; exit 1; }

mkdir -p "$OUT/out"
cd "$FROZEN"
rm -rf output output_backup && mkdir -p output
cp -f "$PACK/configs/geometry_and_probe_mapping.json" "$OUT/"
cp -f "$PACK/configs/C0_effective.json" "$OUT/"

echo "[doc248][C0] start $(date -Is) end=$END_TIME wall=${WALL_S}s"
set +e
timeout --signal=TERM --kill-after=120 "${WALL_S}s" \
  ./test_3d_ophelie_french_natural_convection_frozen_q \
    --client-nc-geometry --client-nc-material \
    --dp=0.010 --melt-height=0.180 \
    --reload-dir="$RELOAD" \
    --doc223-paper-flat-slip-top \
    --em-control=fixed-power --target-absorbed-power=15000 \
    --thermo-em-coupling=periodic --em-update-interval=0.5 --aind=off \
    --enable-thermal-diffusion \
    --doc224-u-blowup-gate=50 \
    --end-time="$END_TIME" --hydro-settle-time=1.0 \
    --state-record-interval=5.0 \
    --meridional-average-start=1.0 \
    --boundary-face-t-mode=local_dt_dz \
    --probes-timeseries="$OUT/probes_timeseries.csv" \
    --summary-json="$OUT/summary.json" \
    2>&1 | tee "$OUT/run.log"
rc=${PIPESTATUS[0]}
set -e
echo "[doc248][C0] exit=$rc $(date -Is)" | tee -a "$OUT/run.log"
echo "$rc" >"$OUT/exit_code.txt"

cp -f output/french_energy_budget.csv output/french_thermo_em_coupling.csv \
      output/french_spatial_stats.csv output/french_T_rz.csv \
      output/french_meridional_timeseries.csv \
      "$OUT/" 2>/dev/null || true
cp -f output/french_T_rz_t*.csv "$OUT/" 2>/dev/null || true
cp -a output/*.vtp "$OUT/out/" 2>/dev/null || true
exit "$rc"
