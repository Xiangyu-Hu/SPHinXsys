#!/usr/bin/env bash
# Doc248 C4 — client stirring with pure rigid prescribed paddle + volumetric heat.
# For long engineering look at stirred flow / temperature trend.
# NOT elastic_driven; NOT stainless deformation claim. Soft hang is C2; stainless elastic is C3.
# Requires 01_generate_particles (stir) first. Shares Reload with C2/C3; separate OUT.
#
# Default short probe 5 s; long hang:
#   END_TIME=60|600 MAX_WALL_HOURS=24 bash 05_run_C4_stir_rigid.sh
set -euo pipefail
SRC=${SRC:-/home/ubuntu/software/SPHinXsys}
BUILD=${BUILD:-$SRC/build-doc248}
PACK=$SRC/docs/ophelie/reference/doc248_client_cases
OUT=${OUT:-$PACK/runs/C4_stir_rigid}
BIN=$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_client_stirring_f1_heated_elastic/bin
RELOAD=${RELOAD:-$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_client_stirring_glass_relax/bin/reload}

END_TIME=${END_TIME:-5}
MAX_WALL_HOURS=${MAX_WALL_HOURS:-4}
OMEGA_RAMP=${OMEGA_RAMP:-2}
HEAT_POWER=${HEAT_POWER:-15000}
C0=${C0:-8}
STATE_RECORD_INTERVAL=${STATE_RECORD_INTERVAL:-1}
SCREEN_EVERY=${SCREEN_EVERY:-10}
CSV_EVERY=${CSV_EVERY:-10}

export ONEAPI_DEVICE_SELECTOR=${ONEAPI_DEVICE_SELECTOR:-cuda:gpu}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-8}
export TBB_NUM_THREADS=${TBB_NUM_THREADS:-8}

[[ -f "$RELOAD/Reload.xml" ]] || { echo "FATAL: no Reload at $RELOAD — run 01_generate_particles.sh"; exit 1; }
[[ -x "$BIN/test_3d_ophelie_client_stirring_f1_heated_elastic" ]] || { echo "FATAL: compile f1"; exit 1; }

mkdir -p "$OUT"
cd "$BIN"
rm -rf output && mkdir -p output
cp -f "$PACK/configs/geometry_and_probe_mapping.json" "$OUT/"
cp -f "$PACK/configs/C4_effective.json" "$OUT/"

ROTOR=./data/client_stir_paddle_solver_m.stl
[[ -f "$ROTOR" ]] || ROTOR=./input/client_stir_paddle_solver_m.stl

echo "[doc248][C4] start $(date -Is) end=$END_TIME omega_ramp=$OMEGA_RAMP c0=$C0 heat=$HEAT_POWER max-wall-h=$MAX_WALL_HOURS"
echo "[doc248][C4] note: rigid_prescribed paddle (--rigid-rotor); not C2 soft elastic / not C3 stainless elastic"
set +e
./test_3d_ophelie_client_stirring_f1_heated_elastic \
  --reload-dir="$RELOAD" \
  --rotor-stl="$ROTOR" \
  --rigid-rotor \
  --end-time="$END_TIME" --omega-ramp="$OMEGA_RAMP" --c0="$C0" \
  --heat-power="$HEAT_POWER" --clamp-top-band=0.05 \
  --state-record-interval="$STATE_RECORD_INTERVAL" \
  --screen-every="$SCREEN_EVERY" --csv-every="$CSV_EVERY" \
  --max-wall-hours="$MAX_WALL_HOURS" \
  --out-dir="$OUT" --summary-json=summary.json \
  2>&1 | tee "$OUT/run.log"
rc=${PIPESTATUS[0]}
set -e
echo "[doc248][C4] exit=$rc $(date -Is)" | tee -a "$OUT/run.log"
echo "$rc" >"$OUT/exit_code.txt"
cp -f output/*.csv "$OUT/" 2>/dev/null || true
cp -a output/*.vtp "$OUT/" 2>/dev/null || true
exit "$rc"
