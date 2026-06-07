#!/usr/bin/env bash
# Stage 10.17 native TEAM7 SI smoke: P0 audit, P1 source, P3 Bz reference (P2 optional).
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
BUILD="${ROOT}/build"
RUN_P2="${RUN_TEAM7_P2:-0}"

cd "${BUILD}"

ninja \
    test_3d_aphi_ck_team7_native_geometry_audit \
    test_3d_aphi_ck_team7_native_source_rhs_audit \
    test_3d_aphi_ck_team7_native_geometry_source_driven_smoke \
    test_3d_aphi_ck_team7_native_geometry_bz_reference_probe

run_one() {
    local test_name="$1"
    local bin="${BUILD}/tests/extra_source_and_tests/3d_examples/${test_name}/bin"
    echo "=== ${test_name} ==="
    (cd "${bin}" && "./${test_name}")
}

run_one test_3d_aphi_ck_team7_native_geometry_audit
run_one test_3d_aphi_ck_team7_native_source_rhs_audit
run_one test_3d_aphi_ck_team7_native_geometry_source_driven_smoke
if [[ "${RUN_P2}" == "1" ]]; then
    ninja test_3d_aphi_ck_team7_native_vacuum_source_b_sanity
    run_one test_3d_aphi_ck_team7_native_vacuum_source_b_sanity
fi
run_one test_3d_aphi_ck_team7_native_geometry_bz_reference_probe

echo "TEAM7 native SI smoke: all requested tests passed."
