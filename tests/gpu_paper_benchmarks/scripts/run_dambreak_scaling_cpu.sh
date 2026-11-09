#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="${BUILD_DIR:-$(cd -- "${SCRIPT_DIR}/../../.." && pwd)/build-host}"
source "${SCRIPT_DIR}/benchmark_run_common.sh"

# Use the native host build for the CPU baseline and remove any inherited SYCL
# selector from the child process environment.
CPU_DEVICE="${SPH_BENCH_DEVICE:-host_cpu}"

run_dambreak_scaling_suite "host_tbb" "${CPU_DEVICE}" ""
