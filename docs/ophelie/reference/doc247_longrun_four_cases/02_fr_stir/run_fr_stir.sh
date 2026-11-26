#!/usr/bin/env bash
# French stirring — Fig8a axisym IC + 60kW/10rpm + skull/crust (Doc247 locked).
# Defaults: short acceptance window end-time=10 / max-wall-hours=1.
# Long runs only if explicitly requested: END_STIR=... WALL_STIR_H=...
set -euo pipefail

SRC=/home/ubuntu/SPHinXsys
BUILD=${BUILD:-$SRC/build}
OUT_ROOT=${OUT_ROOT:-$SRC/docs/ophelie/reference/doc247_longrun_four_cases/overnight_6h}
FR_STIR=$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_stirring_em/bin
RELOAD_FR_STIR=$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_stirring_glass_relax/bin/reload
SEED_STIR=$SRC/docs/ophelie/reference/standard_data/doc247_seeds/SEED_STIR_Fig8a_axisym_ASSUMED_T_rz.csv

WALL_STIR_H=${WALL_STIR_H:-1}
END_STIR=${END_STIR:-10}

export ONEAPI_DEVICE_SELECTOR=${ONEAPI_DEVICE_SELECTOR:-cuda:gpu}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-8}
export TBB_NUM_THREADS=${TBB_NUM_THREADS:-8}

BIN=$FR_STIR/test_3d_ophelie_french_stirring_em
[[ -x "$BIN" ]] || { echo "FATAL: missing $BIN — compile first"; exit 1; }
[[ -f "$SEED_STIR" ]] || { echo "FATAL: missing $SEED_STIR"; exit 1; }
[[ -f "$RELOAD_FR_STIR/Reload.xml" ]] || { echo "FATAL: missing STIR reload"; exit 1; }

OUT=$OUT_ROOT/02_FR_STIR_Fig8a_EM
mkdir -p "$OUT/out" "$OUT/trz_snaps" "$OUT_ROOT"
cd "$FR_STIR"
rm -rf output && mkdir -p output

echo "[doc247][STIR] start $(date -Is) end-time=$END_STIR max-wall-hours=$WALL_STIR_H OUT=$OUT" \
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
rc=${PIPESTATUS[0]}
set -e

echo "[doc247][STIR] exit=$rc $(date -Is)" | tee -a "$OUT_ROOT/overnight.log"
cp -f output/*.csv "$OUT/" 2>/dev/null || true
cp -f output/french_T_rz_t*.csv "$OUT/trz_snaps/" 2>/dev/null || true
cp -a output/*.vtp "$OUT/out/" 2>/dev/null || true
echo "$rc" >"$OUT/exit_code.txt"
exit "$rc"
