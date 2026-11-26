#!/usr/bin/env bash
# Doc248 — compile client NC + stirring stack (run on build machine).
set -euo pipefail
SRC=${SRC:-/home/ubuntu/software/SPHinXsys}
BUILD=${BUILD:-$SRC/build-doc248}
cmake --build "$BUILD" --target \
  test_3d_ophelie_client_nc_glass_relax \
  test_3d_ophelie_client_stirring_glass_relax \
  test_3d_ophelie_french_natural_convection_frozen_q \
  test_3d_ophelie_client_stirring_f1_heated_elastic \
  test_3d_ophelie_french_stirring_em \
  -j"$(nproc)"
echo "[doc248] compile OK $(date -Is)"
