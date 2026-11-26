#!/usr/bin/env bash
# Doc248 C5 — client stir params on french_stirring_em (EM + skull/crust + Boussinesq + host-spin).
# EM: fixed-power target 30000 W = generator ledger used as engineering scaling target
#     (power_identity=generator_input_used_as_EM_scaling_target; NOT η=100% claim; NOT C0 15 kW).
# Cooling: thicken skull vs French 5 mm so net heating P_em−P_loss ≈ 18–20 kW at T0
#     (baseline C5 with delta=5 mm had P_loss≈26 kW → net≈4 kW).
# Requires 01_generate_particles (stir) first. Separate OUT from C2/C3/C4.
#
# Default 60 s net-ΔT probe:
#   END_TIME=60 bash 06_run_C5_stir_em30.sh
# Retune cooling without rebuild:
#   SKULL_DELTA=0.028 FREE_CRUST_DELTA=0.010 END_TIME=60 bash 06_run_C5_stir_em30.sh
set -euo pipefail
SRC=${SRC:-/home/ubuntu/software/SPHinXsys}
BUILD=${BUILD:-$SRC/build-doc248}
PACK=$SRC/docs/ophelie/reference/doc248_client_cases
OUT=${OUT:-$PACK/runs/C5_stir_em30}
BIN=$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_stirring_em/bin
RELOAD=${RELOAD:-$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_client_stirring_glass_relax/bin/reload}
RELAX_DATA=$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_client_stirring_glass_relax/bin/data
F1_DATA=$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_client_stirring_f1_heated_elastic/bin/data

END_TIME=${END_TIME:-60}
MAX_WALL_HOURS=${MAX_WALL_HOURS:-10}
C0=${C0:-8}
STATE_RECORD_INTERVAL=${STATE_RECORD_INTERVAL:-5}
SCREEN_EVERY=${SCREEN_EVERY:-10}
CSV_EVERY=${CSV_EVERY:-10}
EM_UPDATE_INTERVAL=${EM_UPDATE_INTERVAL:-0.5}
TARGET_ABSORBED=${TARGET_ABSORBED:-30000}
# Side/bottom skull thickness [m]. 0.030 → h_eff_side≈30 (was 120 @ 5 mm); 0D net≈18–19 kW @ T0.
SKULL_DELTA=${SKULL_DELTA:-0.030}
SKULL_K=${SKULL_K:-1}
FREE_CRUST_DELTA=${FREE_CRUST_DELTA:-0.005}
FREE_CRUST_K=${FREE_CRUST_K:-1}

export ONEAPI_DEVICE_SELECTOR=${ONEAPI_DEVICE_SELECTOR:-cuda:gpu}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-8}
export TBB_NUM_THREADS=${TBB_NUM_THREADS:-8}

[[ -f "$RELOAD/Reload.xml" ]] || { echo "FATAL: no Reload at $RELOAD — run 01_generate_particles.sh"; exit 1; }
[[ -x "$BIN/test_3d_ophelie_french_stirring_em" ]] || { echo "FATAL: compile french_stirring_em"; exit 1; }

ROTOR=${ROTOR:-}
if [[ -z "$ROTOR" ]]; then
  for cand in \
    "$RELAX_DATA/client_stir_paddle_solver_m.stl" \
    "$F1_DATA/client_stir_paddle_solver_m.stl" \
    "$BIN/data/client_stir_paddle_solver_m.stl"
  do
    [[ -f "$cand" ]] && ROTOR=$cand && break
  done
fi
[[ -n "$ROTOR" && -f "$ROTOR" ]] || { echo "FATAL: client paddle STL not found (set ROTOR=)"; exit 1; }

# Resolve relative OUT/RELOAD/ROTOR before cd into $BIN (relative paths would break after cd).
mkdir -p "$OUT"
OUT=$(cd "$OUT" && pwd)
RELOAD=$(cd "$RELOAD" && pwd)
ROTOR=$(cd "$(dirname "$ROTOR")" && pwd)/$(basename "$ROTOR")
mkdir -p "$OUT/out"
cd "$BIN"
rm -rf output && mkdir -p output
cp -f "$PACK/configs/geometry_and_probe_mapping.json" "$OUT/"
cp -f "$PACK/configs/C5_effective.json" "$OUT/"

echo "[doc248][C5] start $(date -Is) end=$END_TIME target_W=$TARGET_ABSORBED max-wall-h=$MAX_WALL_HOURS"
echo "[doc248][C5] cooling: skull_delta=${SKULL_DELTA}m k=${SKULL_K} free_crust_delta=${FREE_CRUST_DELTA}m (aim net heating ~18–20 kW = P_em−P_loss)"
echo "[doc248][C5] note: french_stirring_em + client-stir; EM scaling target=$TARGET_ABSORBED W (generator ledger); not C0 15 kW; not C4 fake heat"
set +e
./test_3d_ophelie_french_stirring_em \
  --client-stir-geometry --client-stir-material \
  --reload-dir="$RELOAD" \
  --rotor-stl="$ROTOR" \
  --rotation-rpm=10 \
  --em-control=fixed-power --target-absorbed-power="$TARGET_ABSORBED" \
  --thermo-em-coupling=periodic --em-update-interval="$EM_UPDATE_INTERVAL" --aind=off \
  --skull-resistance --skull-delta="$SKULL_DELTA" --skull-k="$SKULL_K" \
  --free-crust-resistance --free-crust-delta="$FREE_CRUST_DELTA" --free-crust-k="$FREE_CRUST_K" \
  --end-time="$END_TIME" --max-wall-hours="$MAX_WALL_HOURS" --c0="$C0" \
  --state-record-interval="$STATE_RECORD_INTERVAL" \
  --screen-every="$SCREEN_EVERY" --csv-every="$CSV_EVERY" \
  --summary-json="$OUT/summary.json" \
  2>&1 | tee "$OUT/run.log"
rc=${PIPESTATUS[0]}
set -e
echo "[doc248][C5] exit=$rc $(date -Is)" | tee -a "$OUT/run.log"
echo "$rc" >"$OUT/exit_code.txt"

cp -f output/french_energy_budget.csv output/french_thermo_em_coupling.csv \
      output/french_spatial_stats.csv output/french_T_rz.csv \
      "$OUT/" 2>/dev/null || true
cp -f output/french_T_rz_t*.csv "$OUT/" 2>/dev/null || true
cp -f output/*.csv "$OUT/" 2>/dev/null || true
cp -a output/*.vtp "$OUT/out/" 2>/dev/null || true
exit "$rc"
