#!/usr/bin/env bash
# Doc248 C2 — client stirring elastic F1. Requires 01_generate_particles (stir) first.
#
# Prior blow-up (~t=0.32 s): stainless E=2e11 → ~500 solid substeps/fluid step,
# tip runaway (e_off→0.11, tip≫ωr). Defaults below use ASSUMED soft E + 2 s ramp
# so multi-second FSI can finish. Soft E is NOT a stainless certificate.
# Override: YOUNGS_E=2e11 OMEGA_RAMP=0.5 for short stainless smoke (expect early abort).
set -euo pipefail
SRC=${SRC:-/home/ubuntu/SPHinXsys}
BUILD=${BUILD:-$SRC/build}
PACK=$SRC/docs/ophelie/reference/doc248_client_cases
OUT=${OUT:-$PACK/runs/C2_stir}
BIN=$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_client_stirring_f1_heated_elastic/bin
RELOAD=${RELOAD:-$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_client_stirring_glass_relax/bin/reload}
END_TIME=${END_TIME:-5}           # aim ≥2–5 s elastic FSI with soft E
MAX_WALL_HOURS=${MAX_WALL_HOURS:-24}
OMEGA_RAMP=${OMEGA_RAMP:-2.0}
YOUNGS_E=${YOUNGS_E:-2e8}         # ASSUMED soft NeoHookean for FSI longevity
SOLID_DAMP=${SOLID_DAMP:-0.5}
RPM=${RPM:-10}
HEAT_POWER=${HEAT_POWER:-15000}
STATE_INTERVAL=${STATE_INTERVAL:-0.25}
SEED_DIR=${SEED_DIR:-$SRC/docs/ophelie/reference/standard_data/doc248_seeds}
INIT_T_CSV=${INIT_T_CSV:-$SEED_DIR/SEED_CLIENT_STIR_Fig8a_T_rz.csv}
INIT_Q_CSV=${INIT_Q_CSV:-$SEED_DIR/SEED_CLIENT_STIR_Fig10_Q_rz.csv}

export ONEAPI_DEVICE_SELECTOR=${ONEAPI_DEVICE_SELECTOR:-cuda:gpu}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-8}
export TBB_NUM_THREADS=${TBB_NUM_THREADS:-8}

[[ -f "$RELOAD/Reload.xml" ]] || { echo "FATAL: no Reload at $RELOAD — run 01_generate_particles.sh"; exit 1; }
[[ -x "$BIN/test_3d_ophelie_client_stirring_f1_heated_elastic" ]] || { echo "FATAL: compile f1"; exit 1; }
[[ -f "$INIT_T_CSV" ]] || { echo "FATAL: missing INIT_T_CSV=$INIT_T_CSV"; exit 1; }
[[ -f "$INIT_Q_CSV" ]] || { echo "FATAL: missing INIT_Q_CSV=$INIT_Q_CSV"; exit 1; }

mkdir -p "$OUT"
cd "$BIN"
rm -rf output && mkdir -p output
cp -f "$PACK/configs/geometry_and_probe_mapping.json" "$OUT/"
cp -f "$PACK/configs/C2_effective.json" "$OUT/"

ROTOR=./data/client_stir_paddle_solver_m.stl
[[ -f "$ROTOR" ]] || ROTOR=./input/client_stir_paddle_solver_m.stl

echo "[doc248][C2] start $(date -Is) end=$END_TIME max-wall-h=$MAX_WALL_HOURS P=$HEAT_POWER"
echo "[doc248][C2] ASSUMED soft E=$YOUNGS_E ramp=$OMEGA_RAMP damp=$SOLID_DAMP rpm=$RPM"
echo "[doc248][C2] ASSUMED IC T=$INIT_T_CSV Q=$INIT_Q_CSV"
set +e
./test_3d_ophelie_client_stirring_f1_heated_elastic \
  --reload-dir="$RELOAD" \
  --rotor-stl="$ROTOR" \
  --end-time="$END_TIME" --omega-ramp="$OMEGA_RAMP" --c0=8 \
  --youngs-modulus="$YOUNGS_E" --solid-damping="$SOLID_DAMP" --rpm="$RPM" \
  --heat-power="$HEAT_POWER" --clamp-top-band=0.05 \
  --init-T-rz-csv="$INIT_T_CSV" \
  --init-Q-rz-csv="$INIT_Q_CSV" \
  --state-record-interval="$STATE_INTERVAL" --screen-every=5 --csv-every=5 \
  --max-wall-hours="$MAX_WALL_HOURS" \
  --out-dir="$OUT" --summary-json=summary.json \
  2>&1 | tee "$OUT/run.log"
rc=${PIPESTATUS[0]}
set -e
echo "[doc248][C2] exit=$rc $(date -Is)" | tee -a "$OUT/run.log"
echo "$rc" >"$OUT/exit_code.txt"
cp -f output/*.csv "$OUT/" 2>/dev/null || true
cp -a output/*.vtp "$OUT/" 2>/dev/null || true
# Keep a copy of last rotor VTP name hint for ParaView stress coloring
ls -1 output/Rotor*.vtp 2>/dev/null | tee "$OUT/rotor_vtp_list.txt" >/dev/null || true
exit "$rc"
