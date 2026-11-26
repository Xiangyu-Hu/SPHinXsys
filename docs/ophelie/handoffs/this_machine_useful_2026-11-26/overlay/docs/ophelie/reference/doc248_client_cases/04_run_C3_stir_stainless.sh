#!/usr/bin/env bash
# Doc248 C3 — client stirring elastic F1 with stainless 304 paddle (E=200 GPa).
# Material-claim path; NOT the engineering soft-hang case (that is C2).
# Requires 01_generate_particles (stir) first. Shares Reload with C2; separate OUT.
#
# Default: short probe window (historical F1 PASS scale). Full-size stainless elastic
# long runs have not historically passed numerical quality — long hang is non-acceptance.
#   END_TIME=60|600 MAX_WALL_HOURS=24 bash 04_run_C3_stir_stainless.sh
set -euo pipefail
SRC=${SRC:-/home/ubuntu/software/SPHinXsys}
BUILD=${BUILD:-$SRC/build-doc248}
PACK=$SRC/docs/ophelie/reference/doc248_client_cases
OUT=${OUT:-$PACK/runs/C3_stir_stainless}
BIN=$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_client_stirring_f1_heated_elastic/bin
RELOAD=${RELOAD:-$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_client_stirring_glass_relax/bin/reload}

END_TIME=${END_TIME:-0.2}
MAX_WALL_HOURS=${MAX_WALL_HOURS:-2}
OMEGA_RAMP=${OMEGA_RAMP:-0.05}
HEAT_POWER=${HEAT_POWER:-15000}
C0=${C0:-8}
# Doc228 S0-scale damping only; do not apply C2 engineering over-damping here.
DAMP_SHAPE=${DAMP_SHAPE:-0.4}
# Tight elastic-offset fuse (m): stainless should stay near rigid; abort early on whip.
E_OFF_ABORT=${E_OFF_ABORT:-0.025}
STATE_RECORD_INTERVAL=${STATE_RECORD_INTERVAL:-0.05}
SCREEN_EVERY=${SCREEN_EVERY:-5}
CSV_EVERY=${CSV_EVERY:-5}

export ONEAPI_DEVICE_SELECTOR=${ONEAPI_DEVICE_SELECTOR:-cuda:gpu}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-8}
export TBB_NUM_THREADS=${TBB_NUM_THREADS:-8}

[[ -f "$RELOAD/Reload.xml" ]] || { echo "FATAL: no Reload at $RELOAD — run 01_generate_particles.sh"; exit 1; }
[[ -x "$BIN/test_3d_ophelie_client_stirring_f1_heated_elastic" ]] || { echo "FATAL: compile f1"; exit 1; }

mkdir -p "$OUT"
cd "$BIN"
rm -rf output && mkdir -p output
cp -f "$PACK/configs/geometry_and_probe_mapping.json" "$OUT/"
cp -f "$PACK/configs/C3_effective.json" "$OUT/"

ROTOR=./data/client_stir_paddle_solver_m.stl
[[ -f "$ROTOR" ]] || ROTOR=./input/client_stir_paddle_solver_m.stl

echo "[doc248][C3] start $(date -Is) end=$END_TIME damp_shape=$DAMP_SHAPE omega_ramp=$OMEGA_RAMP e_off_abort=$E_OFF_ABORT c0=$C0 max-wall-h=$MAX_WALL_HOURS"
echo "[doc248][C3] note: stainless_304_development_preset E=200 GPa (no --youngs-modulus); soft hang is C2"
set +e
# Intentionally omit --youngs-modulus so F1 uses makeClientStainless304Development (2e11 Pa).
./test_3d_ophelie_client_stirring_f1_heated_elastic \
  --reload-dir="$RELOAD" \
  --rotor-stl="$ROTOR" \
  --end-time="$END_TIME" --omega-ramp="$OMEGA_RAMP" --c0="$C0" \
  --solid-damping-shape="$DAMP_SHAPE" \
  --e-off-abort="$E_OFF_ABORT" \
  --heat-power="$HEAT_POWER" --clamp-top-band=0.05 \
  --state-record-interval="$STATE_RECORD_INTERVAL" \
  --screen-every="$SCREEN_EVERY" --csv-every="$CSV_EVERY" \
  --max-wall-hours="$MAX_WALL_HOURS" \
  --out-dir="$OUT" --summary-json=summary.json \
  2>&1 | tee "$OUT/run.log"
rc=${PIPESTATUS[0]}
set -e
echo "[doc248][C3] exit=$rc $(date -Is)" | tee -a "$OUT/run.log"
echo "$rc" >"$OUT/exit_code.txt"
cp -f output/*.csv "$OUT/" 2>/dev/null || true
cp -a output/*.vtp "$OUT/" 2>/dev/null || true
exit "$rc"
