#!/usr/bin/env bash
# Optional SYCL USM regression for matrix-free A–phi style binaries (extra_src).
# Compares SYCL-all vs SYCL-all + graph topology USM + scalar field USM (read + write paths).
# Usage (from build tree):
#   bash ../tests/extra_source_and_tests/3d_examples/matrix_free_aphi_sycl_usm_smoke.sh \
#       ./tests/extra_source_and_tests/3d_examples/test_3d_em_aphi_matrix_free_*/bin/<target> [extra binary args...]
set -euo pipefail
BIN="${1:?Usage: $0 /path/to/binary [binary args...]}"
shift || true
EXTRA=("$@")

# Short smoke defaults when driving TEAM7-like binaries (overridable by caller).
export EM_APHI_STAGGERED_SMOKE_DP="${EM_APHI_STAGGERED_SMOKE_DP:-0.1}"
export EM_APHI_STAGGERED_SMOKE_OUTER_ITERS="${EM_APHI_STAGGERED_SMOKE_OUTER_ITERS:-2}"

extract() { grep -oE "$1" | head -1 || true; }

run_case() {
  local name="$1"
  shift
  echo "=== ${name} ==="
  local out
  out="$(env "$@" "${BIN}" "${EXTRA[@]}" 2>&1)" || true
  echo "${out}" | extract 'residual_ay_l2=[^ ]+' || true
  echo "${out}" | extract 'validation_pass=[^ ]+' || true
  echo "${out}" | extract 'matrix_free_sycl_policy_selftest_pass=[^ ]+' || true
}

run_case "SYCL-all" \
  EM_APHI_MATRIX_FREE_USE_SYCL_GRADIENT=1 \
  EM_APHI_MATRIX_FREE_USE_SYCL_HARMONIC_GRADIENT=1 \
  EM_APHI_MATRIX_FREE_USE_SYCL_JACOBI=1 \
  EM_APHI_MATRIX_FREE_USE_SYCL_LAPLACE_RESIDUAL=1

run_case "SYCL-all + graph topology USM + scalar field USM" \
  EM_APHI_MATRIX_FREE_USE_SYCL_GRADIENT=1 \
  EM_APHI_MATRIX_FREE_USE_SYCL_HARMONIC_GRADIENT=1 \
  EM_APHI_MATRIX_FREE_USE_SYCL_JACOBI=1 \
  EM_APHI_MATRIX_FREE_USE_SYCL_LAPLACE_RESIDUAL=1 \
  EM_APHI_MATRIX_FREE_SYCL_GRAPH_TOPOLOGY_USM=1 \
  EM_APHI_MATRIX_FREE_SYCL_FIELD_VALUES_USM=1

echo "Done."
