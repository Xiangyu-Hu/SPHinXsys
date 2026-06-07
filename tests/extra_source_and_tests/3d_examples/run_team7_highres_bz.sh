#!/usr/bin/env bash
# High-resolution TEAM7 native SI workflow (small air, dp=3 mm).
# Step 1 relax is slow (~8× particles vs dp=6 mm); run interactively or in tmux.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
BUILD="${TEAM7_BUILD_DIR:-${ROOT}/build}"
DP="${TEAM7_DP:-0.003}"
AIR="${TEAM7_AIR:-small}"

echo "=== [1/4] Build particle_generation_em + P3 Bz probe ==="
cd "${BUILD}"
cmake .. -DCMAKE_BUILD_TYPE=Release >/dev/null
ninja particle_generation_em test_3d_aphi_ck_team7_native_geometry_bz_reference_probe

PG_BIN="${BUILD}/tests/extra_source_and_tests/3d_examples/particle_generation_em/bin"
P3_BIN="${BUILD}/tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_team7_native_geometry_bz_reference_probe/bin"

echo "=== [2/4] Relax + reload (dp=${DP} m, air=${AIR}) — may take hours ==="
cd "${PG_BIN}"
./particle_generation_em --team7-dp="${DP}" --team7-air="${AIR}"

echo "=== [3/4] Sync reload into P3 test bin (re-configure copies reload/) ==="
cd "${BUILD}"
cmake .. >/dev/null
ninja test_3d_aphi_ck_team7_native_geometry_bz_reference_probe

echo "=== [4/4] Run Bz reference probe ==="
cd "${P3_BIN}"
./test_3d_aphi_ck_team7_native_geometry_bz_reference_probe
echo "Report: ${P3_BIN}/team7_native_bz_reference_report/"
