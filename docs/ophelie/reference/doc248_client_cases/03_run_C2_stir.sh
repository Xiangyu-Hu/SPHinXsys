#!/usr/bin/env bash
# Doc248 C2 — client stirring elastic F1. Requires 01_generate_particles (stir) first.
# Elastic solid substeps → expect long wall time; raise MAX_WALL_HOURS / END_TIME on GPU box.
set -euo pipefail
SRC=${SRC:-/home/ubuntu/SPHinXsys}
BUILD=${BUILD:-$SRC/build}
PACK=$SRC/docs/ophelie/reference/doc248_client_cases
OUT=${OUT:-$PACK/runs/C2_stir}
BIN=$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_client_stirring_f1_heated_elastic/bin
RELOAD=${RELOAD:-$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_client_stirring_glass_relax/bin/reload}
END_TIME=${END_TIME:-60}          # 10 rev @ 10 rpm; increase for overnight
MAX_WALL_HOURS=${MAX_WALL_HOURS:-24}
OMEGA_RAMP=${OMEGA_RAMP:-0.5}
HEAT_POWER=${HEAT_POWER:-15000}

export ONEAPI_DEVICE_SELECTOR=${ONEAPI_DEVICE_SELECTOR:-cuda:gpu}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-8}
export TBB_NUM_THREADS=${TBB_NUM_THREADS:-8}

[[ -f "$RELOAD/Reload.xml" ]] || { echo "FATAL: no Reload at $RELOAD — run 01_generate_particles.sh"; exit 1; }
[[ -x "$BIN/test_3d_ophelie_client_stirring_f1_heated_elastic" ]] || { echo "FATAL: compile f1"; exit 1; }

mkdir -p "$OUT"
cd "$BIN"
rm -rf output && mkdir -p output
cp -f "$PACK/configs/geometry_and_probe_mapping.json" "$OUT/"
cp -f "$PACK/configs/C2_effective.json" "$OUT/"

ROTOR=./data/client_stir_paddle_solver_m.stl
[[ -f "$ROTOR" ]] || ROTOR=./input/client_stir_paddle_solver_m.stl

echo "[doc248][C2] start $(date -Is) end=$END_TIME max-wall-h=$MAX_WALL_HOURS"
set +e
./test_3d_ophelie_client_stirring_f1_heated_elastic \
  --reload-dir="$RELOAD" \
  --rotor-stl="$ROTOR" \
  --end-time="$END_TIME" --omega-ramp="$OMEGA_RAMP" --c0=8 \
  --heat-power="$HEAT_POWER" --clamp-top-band=0.05 \
  --state-record-interval=2 --screen-every=5 --csv-every=5 \
  --max-wall-hours="$MAX_WALL_HOURS" \
  --out-dir="$OUT" --summary-json=summary.json \
  2>&1 | tee "$OUT/run.log"
rc=${PIPESTATUS[0]}
set -e
echo "[doc248][C2] exit=$rc $(date -Is)" | tee -a "$OUT/run.log"
echo "$rc" >"$OUT/exit_code.txt"
cp -f output/*.csv "$OUT/" 2>/dev/null || true
cp -a output/*.vtp "$OUT/" 2>/dev/null || true
exit "$rc"
