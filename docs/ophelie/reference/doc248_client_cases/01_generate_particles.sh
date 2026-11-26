#!/usr/bin/env bash
# Doc248 — regenerate ALL client particles from scratch (no pre-shipped Reload).
# Deletes existing reload dirs then runs NC + stirring lattice relax.
set -euo pipefail
SRC=${SRC:-/home/ubuntu/SPHinXsys}
BUILD=${BUILD:-$SRC/build}
OUT=${OUT:-$SRC/docs/ophelie/reference/doc248_client_cases/runs}
NC_RELAX=$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_client_nc_glass_relax/bin
STIR_RELAX=$BUILD/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_client_stirring_glass_relax/bin
export ONEAPI_DEVICE_SELECTOR=${ONEAPI_DEVICE_SELECTOR:-cuda:gpu}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-8}
export TBB_NUM_THREADS=${TBB_NUM_THREADS:-8}

need() { [[ -x "$1" ]] || { echo "FATAL: compile first: missing $1"; exit 1; }; }
need "$NC_RELAX/test_3d_ophelie_client_nc_glass_relax"
need "$STIR_RELAX/test_3d_ophelie_client_stirring_glass_relax"

mkdir -p "$OUT/particles_nc" "$OUT/particles_stir"

echo "[doc248][particles] NC relax $(date -Is)"
cd "$NC_RELAX"
rm -rf reload output output_backup client_nc_relax_meta
mkdir -p reload output
./test_3d_ophelie_client_nc_glass_relax \
  --dp=0.010 --melt-height=0.230 \
  --relax-steps=1000 --relax-vtp-every=200 \
  --output-meta-dir="$OUT/particles_nc/meta" \
  2>&1 | tee "$OUT/particles_nc/relax.log"
[[ -f reload/Reload.xml ]] || { echo "FATAL: NC Reload.xml missing"; exit 1; }
cp -a reload "$OUT/particles_nc/"
echo "[doc248][particles] NC Reload OK ($(du -h reload/Reload.xml | awk '{print $1}'))"

echo "[doc248][particles] STIR relax $(date -Is)"
cd "$STIR_RELAX"
rm -rf reload output output_backup
mkdir -p reload output
# Prefer solver-frame paddle STL shipped next to binary
ROTOR=./data/client_stir_paddle_solver_m.stl
[[ -f "$ROTOR" ]] || ROTOR=./input/client_stir_paddle_solver_m.stl
./test_3d_ophelie_client_stirring_glass_relax \
  --bodies=all \
  --dp=0.005 --melt-height=0.230 \
  --relax-steps=1000 --relax-vtp-every=200 \
  --rotor-stl="$ROTOR" \
  2>&1 | tee "$OUT/particles_stir/relax.log"
[[ -f reload/Reload.xml ]] || { echo "FATAL: STIR Reload.xml missing"; exit 1; }
cp -a reload "$OUT/particles_stir/"
echo "[doc248][particles] STIR Reload OK ($(du -h reload/Reload.xml | awk '{print $1}'))"
echo "[doc248][particles] done $(date -Is)"
echo "Use bin/reload (live) or $OUT/particles_*/reload for runs."
