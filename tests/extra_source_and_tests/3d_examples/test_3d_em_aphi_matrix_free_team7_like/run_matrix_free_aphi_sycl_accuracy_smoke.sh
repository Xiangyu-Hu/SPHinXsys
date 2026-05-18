#!/usr/bin/env bash
# Matrix-free A–phi SYCL regression: compare CPU vs SYCL-all residual_ay_l2 and validation_pass.
# Usage: from build directory,
#   bash ../tests/extra_source_and_tests/3d_examples/test_3d_em_aphi_matrix_free_team7_like/run_matrix_free_aphi_sycl_accuracy_smoke.sh ./tests/.../bin/test_3d_em_aphi_matrix_free_team7_like
set -euo pipefail
BIN="${1:?Usage: $0 /path/to/test_3d_em_aphi_matrix_free_team7_like}"
export EM_APHI_STAGGERED_SMOKE_DP="${EM_APHI_STAGGERED_SMOKE_DP:-0.1}"
export EM_APHI_STAGGERED_SMOKE_OUTER_ITERS="${EM_APHI_STAGGERED_SMOKE_OUTER_ITERS:-2}"

extract() { grep -oE "$1" | head -1; }

run_case() {
  local name="$1"
  shift
  echo "=== ${name} ==="
  local out
  out="$(env "$@" "${BIN}" 2>&1)"
  echo "${out}" | extract 'residual_ay_l2=[^ ]+' || true
  echo "${out}" | extract 'validation_pass=[^ ]+' || true
  echo "${out}" | extract 'sycl_value_gradient_downloads=[^ ]+' || true
  echo "${out}" | extract 'matrix_free_sycl_policy_selftest_pass=[^ ]+' || true
  echo "${out}" | extract 'matrix_free_sycl_policy_selftest_fail=[^ ]+' || true
}

run_case "CPU (SYCL operators off)" \
  EM_APHI_MATRIX_FREE_USE_SYCL_GRADIENT=0 \
  EM_APHI_MATRIX_FREE_USE_SYCL_JACOBI=0 \
  EM_APHI_MATRIX_FREE_USE_SYCL_LAPLACE_RESIDUAL=0

run_case "SYCL-all" \
  EM_APHI_MATRIX_FREE_USE_SYCL_GRADIENT=1 \
  EM_APHI_MATRIX_FREE_USE_SYCL_HARMONIC_GRADIENT=1 \
  EM_APHI_MATRIX_FREE_USE_SYCL_JACOBI=1 \
  EM_APHI_MATRIX_FREE_USE_SYCL_LAPLACE_RESIDUAL=1

run_case "SYCL-all + policy snapshot selftest" \
  EM_APHI_MATRIX_FREE_USE_SYCL_GRADIENT=1 \
  EM_APHI_MATRIX_FREE_USE_SYCL_HARMONIC_GRADIENT=1 \
  EM_APHI_MATRIX_FREE_USE_SYCL_JACOBI=1 \
  EM_APHI_MATRIX_FREE_USE_SYCL_LAPLACE_RESIDUAL=1 \
  EM_APHI_MATRIX_FREE_SYCL_POLICY_SELFTEST=1

run_case "SYCL-all + graph topology USM (CSR + edge weights on device)" \
  EM_APHI_MATRIX_FREE_USE_SYCL_GRADIENT=1 \
  EM_APHI_MATRIX_FREE_USE_SYCL_HARMONIC_GRADIENT=1 \
  EM_APHI_MATRIX_FREE_USE_SYCL_JACOBI=1 \
  EM_APHI_MATRIX_FREE_USE_SYCL_LAPLACE_RESIDUAL=1 \
  EM_APHI_MATRIX_FREE_SYCL_GRAPH_TOPOLOGY_USM=1

run_case "SYCL-all + graph topology USM + scalar field USM (Jacobi / graph-Helmholtz)" \
  EM_APHI_MATRIX_FREE_USE_SYCL_GRADIENT=1 \
  EM_APHI_MATRIX_FREE_USE_SYCL_HARMONIC_GRADIENT=1 \
  EM_APHI_MATRIX_FREE_USE_SYCL_JACOBI=1 \
  EM_APHI_MATRIX_FREE_USE_SYCL_LAPLACE_RESIDUAL=1 \
  EM_APHI_MATRIX_FREE_SYCL_GRAPH_TOPOLOGY_USM=1 \
  EM_APHI_MATRIX_FREE_SYCL_FIELD_VALUES_USM=1

run_case "SYCL-all + ReservedDeviceKrylov backend bind" \
  EM_APHI_MATRIX_FREE_USE_SYCL_GRADIENT=1 \
  EM_APHI_MATRIX_FREE_USE_SYCL_HARMONIC_GRADIENT=1 \
  EM_APHI_MATRIX_FREE_USE_SYCL_JACOBI=1 \
  EM_APHI_MATRIX_FREE_USE_SYCL_LAPLACE_RESIDUAL=1 \
  EM_APHI_MATRIX_FREE_LA_BACKEND=1

echo "Done."
