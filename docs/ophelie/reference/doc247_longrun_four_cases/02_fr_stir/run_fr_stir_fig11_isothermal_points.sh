#!/usr/bin/env bash
# Fig.11 acceptance points: isothermal / constant-μ stirring (no long EM thermal window).
# Runs 2–3 short rpm or implied-Re points; records late-window Np vs 47/Re.
# Override: END_PT=20 WALL_H=0.5 RPMS="5 10 20"
set -euo pipefail

SRC=/home/ubuntu/SPHinXsys
BUILD=${BUILD:-$SRC/build}
OUT_ROOT=${OUT_ROOT:-$SRC/docs/ophelie/reference/doc247_longrun_four_cases/overnight_6h/acceptance_pack_short/fig11_isothermal}
FR_STIR=$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_stirring_em/bin
RELOAD_FR_STIR=$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_stirring_glass_relax/bin/reload
SEED_STIR=$SRC/docs/ophelie/reference/standard_data/doc247_seeds/SEED_STIR_Fig8a_axisym_ASSUMED_T_rz.csv

# Low-cost Fig.11 pack = 3 cases (5/10/20 rpm). One script run = all three, sequential.
END_PT=${END_PT:-20}
WALL_H=${WALL_H:-1}   # per-case wall budget; 10 s EM stir was ~5 min, 20 s isothermal should fit
RPMS=${RPMS:-"5 10 20"}

export ONEAPI_DEVICE_SELECTOR=${ONEAPI_DEVICE_SELECTOR:-cuda:gpu}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-8}
export TBB_NUM_THREADS=${TBB_NUM_THREADS:-8}

BIN=$FR_STIR/test_3d_ophelie_french_stirring_em
[[ -x "$BIN" ]] || { echo "FATAL: missing $BIN — compile first"; exit 1; }
[[ -f "$RELOAD_FR_STIR/Reload.xml" ]] || { echo "FATAL: missing STIR reload"; exit 1; }

mkdir -p "$OUT_ROOT"
echo "[doc247][STIR-Fig11] start $(date -Is) END_PT=$END_PT RPMS=$RPMS OUT=$OUT_ROOT" \
  | tee -a "$OUT_ROOT/run_fig11.log"

for rpm in $RPMS; do
  OUT=$OUT_ROOT/rpm_${rpm}
  mkdir -p "$OUT/out"
  cd "$FR_STIR"
  rm -rf output && mkdir -p output
  echo "[doc247][STIR-Fig11] rpm=$rpm start $(date -Is)" | tee -a "$OUT_ROOT/run_fig11.log"
  set +e
  ./test_3d_ophelie_french_stirring_em \
    --reload-dir="$RELOAD_FR_STIR" \
    --init-T-rz-csv="$SEED_STIR" \
    --rotation-rpm="$rpm" \
    --constant-thermal-rheology \
    --thermo-em-coupling=off \
    --aind=off \
    --skull-resistance --skull-delta=0.005 \
    --free-crust-resistance --free-crust-delta=0.005 --free-crust-k=1 \
    --end-time="$END_PT" --max-wall-hours="$WALL_H" \
    --state-record-interval=5 \
    --screen-every=10 --csv-every=10 \
    2>&1 | tee "$OUT/run.log"
  rc=${PIPESTATUS[0]}
  set -e
  echo "[doc247][STIR-Fig11] rpm=$rpm exit=$rc $(date -Is)" | tee -a "$OUT_ROOT/run_fig11.log"
  cp -f output/*.csv "$OUT/" 2>/dev/null || true
  echo "$rc" >"$OUT/exit_code.txt"
  # continue remaining rpm points even if one fails (unattended tmux)
done

echo "[doc247][STIR-Fig11] all done $(date -Is)" | tee -a "$OUT_ROOT/run_fig11.log"
# non-zero if any point failed
any_fail=0
for rpm in $RPMS; do
  rc=$(cat "$OUT_ROOT/rpm_${rpm}/exit_code.txt" 2>/dev/null || echo 1)
  [[ "$rc" == "0" ]] || any_fail=1
done
exit "$any_fail"
