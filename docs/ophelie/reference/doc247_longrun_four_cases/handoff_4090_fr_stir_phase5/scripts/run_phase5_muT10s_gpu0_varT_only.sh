#!/usr/bin/env bash
# Single case: MU_varT_rpm10 on GPU0 (~10 s physical)
set -u
set -o pipefail

cd /home/yongchuan/sphinxsys
export SRC=/home/yongchuan/sphinxsys
export BUILD=$SRC/build
export ONEAPI_DEVICE_SELECTOR=cuda:0
export OMP_NUM_THREADS=8
export TBB_NUM_THREADS=8

BIN=$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_stirring_em/bin
RELOAD=$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_stirring_glass_relax/bin/reload_HZ_z00
# Prefer handoff seeds (true Fig.10 digitization); do not depend on shared doc247_seeds dirty tree.
HO=$SRC/docs/ophelie/reference/doc247_longrun_four_cases/handoff_4090_fr_stir_phase5
SEED_T=$HO/seeds/SEED_STIR_Fig8a_axisym_ASSUMED_T_rz.csv
SEED_Q=$HO/seeds/SEED_STIR_Fig10_Qshape_ASSUMED_rz.csv
OUT=$SRC/docs/ophelie/reference/doc247_longrun_four_cases/overnight_6h/phase5_fr_stir_muT10s
id=MU_varT_rpm10
rpm=10
C=$OUT/$id
mkdir -p "$C/out" "$C/trz_snaps" "$OUT"

if [[ ! -x "$BIN/test_3d_ophelie_french_stirring_em" ]]; then
  echo "ERROR: binary missing: $BIN/test_3d_ophelie_french_stirring_em" >&2
  exit 1
fi
if [[ ! -f "$RELOAD/Reload.xml" ]]; then
  echo "ERROR: Reload.xml missing under $RELOAD" >&2
  exit 1
fi

echo "[muT] ONEAPI_DEVICE_SELECTOR=$ONEAPI_DEVICE_SELECTOR case=$id rpm=$rpm start $(date -Is)" | tee -a "$OUT/run_matrix.log"
nvidia-smi -i 0 --query-gpu=index,utilization.gpu,memory.used --format=csv || true

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
  --end-time=10 --max-wall-hours=3 \
  --state-record-interval=5 --screen-every=20 --csv-every=10 \
  --temperature-dependent-properties \
  2>&1 | tee "$C/run.log"
rc=${PIPESTATUS[0]}
set -e

echo "[muT] $id exit=$rc $(date -Is)" | tee -a "$OUT/run_matrix.log"
cp -f output/*.csv "$C/" 2>/dev/null || true
cp -f output/french_T_rz_t*.csv "$C/trz_snaps/" 2>/dev/null || true
cp -a output/*.vtp "$C/out/" 2>/dev/null || true
echo "$rc" >"$C/exit_code.txt"
printf '%s\n' "case_id=$id rpm=$rpm end=10 extras=--temperature-dependent-properties" >"$C/case_meta.txt"
echo "[muT] done $(date -Is) exit=$rc"
exit "$rc"
