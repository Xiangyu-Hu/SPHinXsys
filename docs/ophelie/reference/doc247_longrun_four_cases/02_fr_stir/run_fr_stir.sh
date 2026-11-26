#!/usr/bin/env bash
# French stirring — Fig8a T IC + true Fig.10 Q-correct + 60 kW / 10 rpm + skull/crust.
# Paddle: rigid Simbody prescribed rotation (NOT NeoHookean elastic; that is client C2).
# Defaults: short acceptance window end-time=10 / max-wall-hours=1.
# Long runs: END_STIR=... WALL_STIR_H=...
# Rheology: default constant (μ=4); for μ(T) use run_fr_stir_muT_varT10s.sh.
set -euo pipefail

SRC=${SRC:-/home/ubuntu/SPHinXsys}
BUILD=${BUILD:-$SRC/build}
OUT_ROOT=${OUT_ROOT:-$SRC/docs/ophelie/reference/doc247_longrun_four_cases/overnight_6h}
FR_STIR=$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_stirring_em/bin
RELOAD_FR_STIR=${RELOAD:-$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_stirring_glass_relax/bin/reload}
SEED_T=${SEED_T:-$SRC/docs/ophelie/reference/standard_data/doc247_seeds/SEED_STIR_Fig8a_axisym_ASSUMED_T_rz.csv}
SEED_Q=${SEED_Q:-$SRC/docs/ophelie/reference/standard_data/doc247_seeds/SEED_STIR_Fig10_Qshape_ASSUMED_rz.csv}

WALL_STIR_H=${WALL_STIR_H:-1}
END_STIR=${END_STIR:-10}
RPM=${RPM:-10}
PABS=${PABS:-60000}
BLEND=${BLEND:-0.8}
EVERY=${EVERY:-1}

export ONEAPI_DEVICE_SELECTOR=${ONEAPI_DEVICE_SELECTOR:-cuda:gpu}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-8}
export TBB_NUM_THREADS=${TBB_NUM_THREADS:-8}

BIN=$FR_STIR/test_3d_ophelie_french_stirring_em
[[ -x "$BIN" ]] || { echo "FATAL: missing $BIN — compile first"; exit 1; }
[[ -f "$SEED_T" && -f "$SEED_Q" ]] || { echo "FATAL: missing seeds"; exit 1; }
[[ -f "$RELOAD_FR_STIR/Reload.xml" ]] || { echo "FATAL: missing STIR reload"; exit 1; }

OUT=${OUT:-$OUT_ROOT/02_FR_STIR_Fig8a_Fig10Q}
mkdir -p "$OUT/out" "$OUT/trz_snaps" "$OUT_ROOT"
cd "$FR_STIR"
rm -rf output && mkdir -p output

echo "[doc247][STIR] start $(date -Is) end=$END_STIR rpm=$RPM P=$PABS blend=$BLEND OUT=$OUT" \
  | tee -a "$OUT_ROOT/overnight.log"

set +e
./test_3d_ophelie_french_stirring_em \
  --reload-dir="$RELOAD_FR_STIR" \
  --init-T-rz-csv="$SEED_T" \
  --init-Q-rz-csv="$SEED_Q" \
  --fig10-q-correct --fig10-q-correct-blend="$BLEND" --fig10-q-correct-every="$EVERY" \
  --rotation-rpm="$RPM" \
  --em-control=fixed-power --target-absorbed-power="$PABS" \
  --thermo-em-coupling=periodic --em-update-interval=0.5 \
  --aind=off \
  --skull-resistance --skull-delta=0.005 \
  --free-crust-resistance --free-crust-delta=0.005 --free-crust-k=1 \
  --end-time="$END_STIR" --max-wall-hours="$WALL_STIR_H" \
  --state-record-interval=5 \
  --screen-every=10 --csv-every=10 \
  --constant-thermal-rheology \
  2>&1 | tee "$OUT/run.log"
rc=${PIPESTATUS[0]}
set -e

echo "[doc247][STIR] exit=$rc $(date -Is)" | tee -a "$OUT_ROOT/overnight.log"
cp -f output/*.csv "$OUT/" 2>/dev/null || true
cp -f output/french_T_rz_t*.csv "$OUT/trz_snaps/" 2>/dev/null || true
cp -a output/*.vtp "$OUT/out/" 2>/dev/null || true
echo "$rc" >"$OUT/exit_code.txt"
exit "$rc"
