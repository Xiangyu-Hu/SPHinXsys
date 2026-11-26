#!/usr/bin/env bash
# Phase5 muT 10s matrix on GPU0: const vs varT vs varT@20rpm
set -u
set -o pipefail

cd /home/yongchuan/sphinxsys
export SRC=/home/yongchuan/sphinxsys
export BUILD=$SRC/build
export ONEAPI_DEVICE_SELECTOR=cuda:0
export OMP_NUM_THREADS=8
export TBB_NUM_THREADS=8

BIN=$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_stirring_em/bin
# 基线桨高 Δz=0（与 SW_rpm10 / HZ_z00 同口径）；默认 bin/reload 已不存在
RELOAD=$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_stirring_glass_relax/bin/reload_HZ_z00
HO=$SRC/docs/ophelie/reference/doc247_longrun_four_cases/handoff_4090_fr_stir_phase5
SEED_T=$HO/seeds/SEED_STIR_Fig8a_axisym_ASSUMED_T_rz.csv
SEED_Q=$HO/seeds/SEED_STIR_Fig10_Qshape_ASSUMED_rz.csv
OUT=$SRC/docs/ophelie/reference/doc247_longrun_four_cases/overnight_6h/phase5_fr_stir_muT10s
mkdir -p "$OUT"

if [[ ! -x "$BIN/test_3d_ophelie_french_stirring_em" ]]; then
  echo "ERROR: binary missing: $BIN/test_3d_ophelie_french_stirring_em" >&2
  exit 1
fi
if [[ ! -d "$RELOAD" ]]; then
  echo "ERROR: reload dir missing: $RELOAD" >&2
  exit 1
fi

run_one () {
  local id="$1" rpm="$2"; shift 2
  local extra=("$@")
  local C=$OUT/$id
  mkdir -p "$C/out" "$C/trz_snaps"
  echo "[muT] $id rpm=$rpm start $(date -Is) extra=${extra[*]:-}" | tee -a "$OUT/run_matrix.log"
  cd "$BIN"
  rm -rf output && mkdir -p output
  set +e
  ./test_3d_ophelie_french_stirring_em \
    --reload-dir="$RELOAD" \
    --init-T-rz-csv="$SEED_T" \
    --init-Q-rz-csv="$SEED_Q" \
    --fig10-q-correct --fig10-q-correct-blend=0.8 --fig10-q-correct-every=1 \
    --rotation-rpm="$rpm" \
    --em-control=fixed-power --target-absorbed-power=60000 \
    --thermo-em-coupling=periodic --em-update-interval=0.5 \
    --aind=off \
    --skull-resistance --skull-delta=0.005 \
    --free-crust-resistance --free-crust-delta=0.005 --free-crust-k=1 \
    --end-time=10 --max-wall-hours=1 \
    --state-record-interval=5 --screen-every=20 --csv-every=10 \
    "${extra[@]}" \
    2>&1 | tee "$C/run.log"
  local rc=${PIPESTATUS[0]}
  set -e
  echo "[muT] $id exit=$rc $(date -Is)" | tee -a "$OUT/run_matrix.log"
  cp -f output/*.csv "$C/" 2>/dev/null || true
  cp -f output/french_T_rz_t*.csv "$C/trz_snaps/" 2>/dev/null || true
  cp -a output/*.vtp "$C/out/" 2>/dev/null || true
  echo "$rc" >"$C/exit_code.txt"
  printf '%s\n' "case_id=$id rpm=$rpm end=10 extras=${extra[*]:-}" >"$C/case_meta.txt"
}

echo "[muT] ONEAPI_DEVICE_SELECTOR=$ONEAPI_DEVICE_SELECTOR OUT=$OUT"
nvidia-smi --query-gpu=index,utilization.gpu,memory.used --format=csv || true

# ① 常粘对照
run_one MU_const_rpm10 10 --constant-thermal-rheology

# ② 变粘/变物性基线
run_one MU_varT_rpm10 10 --temperature-dependent-properties

# ③ 变粘 + 更高转速
run_one MU_varT_rpm20 20 --temperature-dependent-properties

echo "[muT] all done $(date -Is)"
