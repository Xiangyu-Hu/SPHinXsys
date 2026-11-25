#!/usr/bin/env bash
# Minimal SPDB slot-hold smoke (extra tests only; does not touch src/).
# Default: one 60 s physical step at dp = 1 mm. Not 7200 s.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
BIN="${ROOT}/build-cpu-icpx/tests/extra_source_and_tests/3d_examples/test_3d_spdb_geometry/bin/test_3d_spdb_geometry"
DP="${DP:-0.001}"
PHYSICAL_DT="${PHYSICAL_DT:-60}"
HOLD_STEPS="${HOLD_STEPS:-1}"
STIFFNESS="${STIFFNESS:-1}"
END_TIME="${END_TIME:-0.0015}"
OUT="${OUT:-${ROOT}/docs/spdb/results/slot_hold_dt${PHYSICAL_DT}_n${HOLD_STEPS}}"

if [[ ! -x "${BIN}" ]]; then
  echo "Missing binary: ${BIN}" >&2
  echo "Build with SPHINXSYS_BUILD_EXTRA_SOURCE_AND_TESTS=ON, then:" >&2
  echo "  cmake --build ${ROOT}/build-cpu-icpx --target test_3d_spdb_geometry -j\$(nproc)" >&2
  exit 1
fi

mkdir -p "${OUT}"
cd "${OUT}"
echo "Running slot_hold in ${OUT} (physical=${PHYSICAL_DT}s x ${HOLD_STEPS}; not claimed as 7200 unless n*dt=7200)"
stdbuf -oL -eL "${BIN}" \
  --slot_hold \
  --dp "${DP}" \
  --physical_dt "${PHYSICAL_DT}" \
  --hold_steps "${HOLD_STEPS}" \
  --contact_stiffness_scale "${STIFFNESS}" \
  --end_time "${END_TIME}" \
  | tee run.log

echo "Done. Metrics: ${OUT}/spdb_slot_hold_*.txt"
