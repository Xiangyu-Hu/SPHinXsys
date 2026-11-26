#!/usr/bin/env bash
# French stir — constant rheology μ≡4 Pa·s @ T0 (10 rpm × 10 s).
# Rigid Simbody paddle (prescribed rotation). ASSUMED engineering demo — not literature PASS.
set -euo pipefail

_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC=${SRC:-$(cd "$_SCRIPT_DIR/../../../../.." && pwd)}
BUILD=${BUILD:-$SRC/build}
OUT_ROOT=${OUT_ROOT:-$SRC/docs/ophelie/reference/doc247_longrun_four_cases/overnight_6h/phase5_fr_stir_muT10s}

FR_STIR=$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_stirring_em/bin
RELOAD=${RELOAD:-$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_stirring_glass_relax/bin/reload}
SEED_T=${SEED_T:-$SRC/docs/ophelie/reference/standard_data/doc247_seeds/SEED_STIR_Fig8a_axisym_ASSUMED_T_rz.csv}
SEED_Q=${SEED_Q:-$SRC/docs/ophelie/reference/standard_data/doc247_seeds/SEED_STIR_Fig10_Qshape_ASSUMED_rz.csv}

END_PT=${END_PT:-10}
WALL_H=${WALL_H:-3}
RPM=${RPM:-10}
PABS=${PABS:-60000}
BLEND=${BLEND:-0.8}
EVERY=${EVERY:-1}
id=${CASE_ID:-MU_const_rpm10}

export ONEAPI_DEVICE_SELECTOR=${ONEAPI_DEVICE_SELECTOR:-cuda:gpu}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-8}
export TBB_NUM_THREADS=${TBB_NUM_THREADS:-8}

BIN=$FR_STIR/test_3d_ophelie_french_stirring_em
[[ -x "$BIN" ]] || { echo "FATAL: missing $BIN — compile first"; exit 1; }
[[ -f "$RELOAD/Reload.xml" ]] || { echo "FATAL: missing Reload at $RELOAD"; exit 1; }
[[ -f "$SEED_T" && -f "$SEED_Q" ]] || { echo "FATAL: missing seeds"; exit 1; }

C=$OUT_ROOT/$id
mkdir -p "$C/out" "$C/trz_snaps" "$OUT_ROOT"
cd "$FR_STIR"
rm -rf output && mkdir -p output

echo "[muT-const] $id rpm=$RPM start $(date -Is)" | tee -a "$OUT_ROOT/run_matrix.log"
set +e
./test_3d_ophelie_french_stirring_em \
  --reload-dir="$RELOAD" \
  --init-T-rz-csv="$SEED_T" \
  --init-Q-rz-csv="$SEED_Q" \
  --fig10-q-correct --fig10-q-correct-blend="$BLEND" --fig10-q-correct-every="$EVERY" \
  --rotation-rpm="$RPM" \
  --em-control=fixed-power --target-absorbed-power="$PABS" \
  --thermo-em-coupling=periodic --em-update-interval=0.5 \
  --aind=off \
  --skull-resistance --skull-delta=0.005 \
  --free-crust-resistance --free-crust-delta=0.005 --free-crust-k=1 \
  --end-time="$END_PT" --max-wall-hours="$WALL_H" \
  --state-record-interval=5 --screen-every=20 --csv-every=10 \
  --constant-thermal-rheology \
  2>&1 | tee "$C/run.log"
rc=${PIPESTATUS[0]}
set -e

echo "[muT-const] $id exit=$rc $(date -Is)" | tee -a "$OUT_ROOT/run_matrix.log"
cp -f output/*.csv "$C/" 2>/dev/null || true
cp -f output/french_T_rz_t*.csv "$C/trz_snaps/" 2>/dev/null || true
cp -a output/*.vtp "$C/out/" 2>/dev/null || true
echo "$rc" >"$C/exit_code.txt"
printf '%s\n' "case_id=$id rpm=$RPM end=$END_PT extras=--constant-thermal-rheology" >"$C/case_meta.txt"
exit "$rc"
