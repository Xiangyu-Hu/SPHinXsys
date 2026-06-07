#!/usr/bin/env bash
# Copy SI reload + input + reference into all native TEAM7 test bin directories.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
BUILD="${ROOT}/build"
SRC="${BUILD}/tests/extra_source_and_tests/3d_examples/particle_generation_em/bin"
REF="${ROOT}/tests/extra_source_and_tests/3d_examples/reference_data"

if [[ ! -d "${SRC}/reload" ]]; then
    echo "Missing ${SRC}/reload — run particle_generation_em first." >&2
    exit 1
fi

TESTS=(
    test_3d_aphi_ck_team7_native_geometry_audit
    test_3d_aphi_ck_team7_native_source_rhs_audit
    test_3d_aphi_ck_team7_native_vacuum_source_b_sanity
    test_3d_aphi_ck_team7_native_geometry_bz_reference_probe
    test_3d_aphi_ck_team7_native_geometry_source_driven_smoke
    test_3d_aphi_ck_team7_native_joule_post_smoke
)

for T in "${TESTS[@]}"; do
    DST="${BUILD}/tests/extra_source_and_tests/3d_examples/${T}/bin"
    mkdir -p "${DST}"
    cp -r "${SRC}/reload" "${SRC}/input" "${DST}/"
    if [[ -d "${REF}" ]]; then
        cp -r "${REF}" "${DST}/"
    fi
    echo "synced: ${T}"
done

echo "Done. Re-run cmake configure to pick up reload via team7_native_reload_sync.cmake."
