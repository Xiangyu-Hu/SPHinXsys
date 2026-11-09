#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="${BUILD_DIR:-$(cd -- "${SCRIPT_DIR}/../../.." && pwd)/build-host}"
source "${SCRIPT_DIR}/benchmark_run_common.sh"

# This launcher intentionally removes ONEAPI_DEVICE_SELECTOR. The CPU baseline
# is the native SPHINXSYS_USE_SYCL=OFF host build, not a SYCL CPU device.
CPU_DEVICE="${SPH_BENCH_DEVICE:-host_cpu}"

run_verification_suite "host_tbb" "${CPU_DEVICE}" ""
