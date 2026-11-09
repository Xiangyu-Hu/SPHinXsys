#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="${BUILD_DIR:-$(cd -- "${SCRIPT_DIR}/../../.." && pwd)/build-sycl}"
source "${SCRIPT_DIR}/benchmark_run_common.sh"

GPU_SELECTOR="${ONEAPI_DEVICE_SELECTOR:-cuda:gpu}"
GPU_DEVICE="${SPH_BENCH_DEVICE:-${GPU_SELECTOR}}"

run_verification_suite "sycl" "${GPU_DEVICE}" "${GPU_SELECTOR}"
