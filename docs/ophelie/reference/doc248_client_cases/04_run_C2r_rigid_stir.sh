#!/usr/bin/env bash
# Doc248 C2r — client stirring with prescribed rigid paddle + volumetric heat.
# Prefer this over C2 elastic when chasing ~1400 K (no solid substeps).
# Requires 01_generate_particles.sh (stir) first.
set -euo pipefail
SRC=${SRC:-/home/ubuntu/SPHinXsys}
BUILD=${BUILD:-$SRC/build}
PACK=$SRC/docs/ophelie/reference/doc248_client_cases
OUT=${OUT:-$PACK/runs/C2r_rigid_stir}
BIN=$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_client_stirring_rigid_heated/bin
RELOAD=${RELOAD:-$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_client_stirring_glass_relax/bin/reload}
END_TIME=${END_TIME:-7200}
MAX_WALL_HOURS=${MAX_WALL_HOURS:-48}
OMEGA_RAMP=${OMEGA_RAMP:-0.5}
HEAT_POWER=${HEAT_POWER:-15000}
TARGET_T_K=${TARGET_T_K:-1400}

export ONEAPI_DEVICE_SELECTOR=${ONEAPI_DEVICE_SELECTOR:-cuda:gpu}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-8}
export TBB_NUM_THREADS=${TBB_NUM_THREADS:-8}

[[ -f "$RELOAD/Reload.xml" ]] || { echo "FATAL: no Reload at $RELOAD — run 01_generate_particles.sh"; exit 1; }
[[ -x "$BIN/test_3d_ophelie_client_stirring_rigid_heated" ]] || {
  echo "FATAL: missing rigid binary — cmake --build \$BUILD --target test_3d_ophelie_client_stirring_rigid_heated"
  exit 1
}

mkdir -p "$OUT"
cd "$BIN"
rm -rf output && mkdir -p output
cp -f "$PACK/configs/geometry_and_probe_mapping.json" "$OUT/"
cp -f "$PACK/configs/C2r_rigid_effective.json" "$OUT/"

ROTOR=./data/client_stir_paddle_solver_m.stl
[[ -f "$ROTOR" ]] || ROTOR=./input/client_stir_paddle_solver_m.stl

echo "[doc248][C2r] start $(date -Is) end=$END_TIME max-wall-h=$MAX_WALL_HOURS target_T=$TARGET_T_K P=$HEAT_POWER"
set +e
./test_3d_ophelie_client_stirring_rigid_heated \
  --reload-dir="$RELOAD" \
  --rotor-stl="$ROTOR" \
  --end-time="$END_TIME" --omega-ramp="$OMEGA_RAMP" --c0=8 \
  --heat-power="$HEAT_POWER" --target-T-K="$TARGET_T_K" \
  --state-record-interval=5 --screen-every=20 --csv-every=20 \
  --max-wall-hours="$MAX_WALL_HOURS" \
  --out-dir="$OUT" --summary-json=summary.json \
  2>&1 | tee "$OUT/run.log"
rc=${PIPESTATUS[0]}
set -e
echo "[doc248][C2r] exit=$rc $(date -Is)" | tee -a "$OUT/run.log"
echo "$rc" >"$OUT/exit_code.txt"
cp -f output/*.csv "$OUT/" 2>/dev/null || true
cp -a output/*.vtp "$OUT/" 2>/dev/null || true
exit "$rc"
