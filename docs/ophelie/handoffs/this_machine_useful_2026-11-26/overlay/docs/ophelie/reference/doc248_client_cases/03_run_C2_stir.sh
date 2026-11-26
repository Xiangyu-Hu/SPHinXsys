#!/usr/bin/env bash
# Doc248 C2 — client stirring elastic F1.
# engineering_soft_hang; stainless material path is C3 (04_run_C3_stir_stainless.sh).
# Default: short engineering stability window with softened E + stronger solid damping.
# Requires 01_generate_particles (stir) first. Raise END_TIME for overnight hang
# only after short-window diverged=0 (not a stainless stiffness claim).
set -euo pipefail
SRC=${SRC:-/home/ubuntu/software/SPHinXsys}
BUILD=${BUILD:-$SRC/build-doc248}
PACK=$SRC/docs/ophelie/reference/doc248_client_cases
OUT=${OUT:-$PACK/runs/C2_stir}
BIN=$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_client_stirring_f1_heated_elastic/bin
RELOAD=${RELOAD:-$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_client_stirring_glass_relax/bin/reload}

# Short-window defaults (stabilize through ω ramp). Override for long hang:
#   END_TIME=60|600 MAX_WALL_HOURS=24 bash 03_run_C2_stir.sh
END_TIME=${END_TIME:-5}
MAX_WALL_HOURS=${MAX_WALL_HOURS:-4}
OMEGA_RAMP=${OMEGA_RAMP:-5}
HEAT_POWER=${HEAT_POWER:-15000}
C0=${C0:-8}
# Softened development modulus (Pa). Not stainless 200 GPa; must label in summary.
YOUNGS=${YOUNGS:-5e8}
# (shape/4)·√(ρE)·dp; raised vs Doc228 S0 0.4 to suppress free-paddle whip after ramp.
DAMP_SHAPE=${DAMP_SHAPE:-1.2}
# Engineering e_off fuse (m); default in binary is max(10*dp,0.05) if unset.
E_OFF_ABORT=${E_OFF_ABORT:-0.05}
STATE_RECORD_INTERVAL=${STATE_RECORD_INTERVAL:-0.5}
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
cp -f "$PACK/configs/C2_effective.json" "$OUT/"

ROTOR=./data/client_stir_paddle_solver_m.stl
[[ -f "$ROTOR" ]] || ROTOR=./input/client_stir_paddle_solver_m.stl

echo "[doc248][C2] start $(date -Is) end=$END_TIME youngs=$YOUNGS damp_shape=$DAMP_SHAPE omega_ramp=$OMEGA_RAMP e_off_abort=$E_OFF_ABORT c0=$C0 max-wall-h=$MAX_WALL_HOURS"
echo "[doc248][C2] note: YOUNGS=$YOUNGS is development_softened_not_stainless unless you pass stainless 2e11"
set +e
./test_3d_ophelie_client_stirring_f1_heated_elastic \
  --reload-dir="$RELOAD" \
  --rotor-stl="$ROTOR" \
  --end-time="$END_TIME" --omega-ramp="$OMEGA_RAMP" --c0="$C0" \
  --youngs-modulus="$YOUNGS" \
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
echo "[doc248][C2] exit=$rc $(date -Is)" | tee -a "$OUT/run.log"
echo "$rc" >"$OUT/exit_code.txt"
cp -f output/*.csv "$OUT/" 2>/dev/null || true
cp -a output/*.vtp "$OUT/" 2>/dev/null || true
exit "$rc"
