#!/usr/bin/env bash
# Uniform TEAM7: dp=6 mm, legacy (slightly larger) air box — particle relax optional, EM solve + VTP.
set -euo pipefail

ROOT="/home/yongchuan/sphinxsys"
BUILD="${ROOT}/build"
PGEN_BIN="${BUILD}/tests/extra_source_and_tests/3d_examples/particle_generation_em/bin"
PROBE_BIN="${BUILD}/tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_team7_native_geometry_bz_reference_probe/bin"
CASE_ID="uniform_legacy_dp6"

cmake --build "${BUILD}" --target particle_generation_em test_3d_aphi_ck_team7_native_geometry_bz_reference_probe -j"$(nproc)"

if [[ "${SKIP_RELAX:-0}" != "1" ]]; then
  cd "${PGEN_BIN}"
  ./particle_generation_em --team7-case=uniform-legacy-dp6
fi

cd "${PROBE_BIN}"
echo "=== TEAM7 Bz probe: reload_case=${CASE_ID} preset=uniform-legacy-dp6 ==="
./test_3d_aphi_ck_team7_native_geometry_bz_reference_probe \
  --team7-case=uniform-legacy-dp6 \
  --team7-reload-case="${CASE_ID}" \
  --team7-write-vtp

VTP_DIR="${PROBE_BIN}/output"
ARCHIVE="${PROBE_BIN}/output_uniform_legacy_dp6_vtp"
rm -rf "${ARCHIVE}"
mkdir -p "${ARCHIVE}"
cp -a "${VTP_DIR}"/*_ite_0000000000.vtp "${ARCHIVE}/" 2>/dev/null || true
echo "VTP archived under ${ARCHIVE}"
