#!/usr/bin/env bash
# Fig.11 soft multi-scheme short sweep (凑点用；非论文 PASS).
# Default: END=12 s, sequential variants with --zero-q (true no-heat).
# Override: END_PT=12 SCHEMES="zq10 zq20 zq5 zq15" WALL_H=0.5
set -euo pipefail

SRC=/home/ubuntu/SPHinXsys
BUILD=${BUILD:-$SRC/build}
OUT_ROOT=${OUT_ROOT:-$SRC/docs/ophelie/reference/doc247_longrun_four_cases/overnight_6h/acceptance_pack_short/fig11_soft_sweep}
FR_STIR=$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_stirring_em/bin
RELOAD_FR_STIR=$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_stirring_glass_relax/bin/reload
SEED_STIR=$SRC/docs/ophelie/reference/standard_data/doc247_seeds/SEED_STIR_Fig8a_axisym_ASSUMED_T_rz.csv

END_PT=${END_PT:-12}
WALL_H=${WALL_H:-0.5}
SCHEMES=${SCHEMES:-"zq10 zq20 zq5 zq15"}

export ONEAPI_DEVICE_SELECTOR=${ONEAPI_DEVICE_SELECTOR:-cuda:gpu}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-8}
export TBB_NUM_THREADS=${TBB_NUM_THREADS:-8}

BIN=$FR_STIR/test_3d_ophelie_french_stirring_em
[[ -x "$BIN" ]] || { echo "FATAL: missing $BIN"; exit 1; }
[[ -f "$RELOAD_FR_STIR/Reload.xml" ]] || { echo "FATAL: missing STIR reload"; exit 1; }

mkdir -p "$OUT_ROOT"
echo "[doc247][Fig11-sweep] start $(date -Is) END_PT=$END_PT SCHEMES=$SCHEMES OUT=$OUT_ROOT" \
  | tee -a "$OUT_ROOT/run_sweep.log"

run_one() {
  local name=$1 rpm=$2
  shift 2
  local extra=("$@")
  # rpm=5 blows at default acoustic_cfl=0.55 (~t=0.02); soft fix is 0.35 (no rebuild needed if binary has --acoustic-cfl=).
  local rpm_i=${rpm%.*}
  if [[ "$rpm_i" -le 5 ]]; then
    extra+=(--acoustic-cfl=0.35)
  fi
  local OUT=$OUT_ROOT/$name
  mkdir -p "$OUT/out"
  cd "$FR_STIR"
  rm -rf output && mkdir -p output
  echo "[doc247][Fig11-sweep] $name rpm=$rpm start $(date -Is) extra=${extra[*]:-}" \
    | tee -a "$OUT_ROOT/run_sweep.log"
  set +e
  ./test_3d_ophelie_french_stirring_em \
    --reload-dir="$RELOAD_FR_STIR" \
    --init-T-rz-csv="$SEED_STIR" \
    --rotation-rpm="$rpm" \
    --constant-thermal-rheology \
    --thermo-em-coupling=off \
    --zero-q \
    --aind=off \
    --skull-resistance --skull-delta=0.005 \
    --free-crust-resistance --free-crust-delta=0.005 --free-crust-k=1 \
    --end-time="$END_PT" --max-wall-hours="$WALL_H" \
    --state-record-interval=5 \
    --screen-every=20 --csv-every=10 \
    --no-state-recording \
    "${extra[@]}" \
    2>&1 | tee "$OUT/run.log"
  local rc=${PIPESTATUS[0]}
  set -e
  echo "[doc247][Fig11-sweep] $name exit=$rc $(date -Is)" | tee -a "$OUT_ROOT/run_sweep.log"
  cp -f output/*.csv "$OUT/" 2>/dev/null || true
  echo "$rc" >"$OUT/exit_code.txt"
  echo "$rpm" >"$OUT/rpm.txt"
}

for s in $SCHEMES; do
  case "$s" in
    zq10) run_one zq10 10 ;;
    zq20) run_one zq20 20 ;;
    zq5)  run_one zq5 5 ;;
    zq15) run_one zq15 15 ;;
    zq8)  run_one zq8 8 ;;
    zq12) run_one zq12 12 ;;
    *) echo "skip unknown scheme $s" | tee -a "$OUT_ROOT/run_sweep.log" ;;
  esac
done

echo "[doc247][Fig11-sweep] all done $(date -Is)" | tee -a "$OUT_ROOT/run_sweep.log"
python3 "$SRC/docs/ophelie/reference/doc247_longrun_four_cases/02_fr_stir/post_fig11_soft_compare.py" \
  --root "$OUT_ROOT" || true
