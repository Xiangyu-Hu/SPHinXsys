#!/usr/bin/env bash
# Short GPU smoke over all six paper cases to validate the shared summary schema.
# Runs in the foreground and always returns to the shell when finished
# (no nohup / tail -f / interactive wait).
#
# Example:
#   ./tests/gpu_paper_benchmarks/scripts/run_gpu_schema_smoke.sh
#
set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../../.." && pwd)"
BUILD_DIR="${BUILD_DIR:-${REPO_ROOT}/build-sycl}"

CAMPAIGN_ID="${CAMPAIGN_ID:-schema-smoke-$(date -u +%Y%m%dT%H%M%SZ)}"
RESULT_ROOT="${RESULT_ROOT:-${BUILD_DIR}/benchmark_results/${CAMPAIGN_ID}}"
REPETITIONS="${REPETITIONS:-1}"
END_TIME_MODE="${END_TIME_MODE:-short}"
SHORT_END_TIME="${SHORT_END_TIME:-0.01}"
OUTPUT="${OUTPUT:-off}"
CONFIG="${CONFIG:-Release}"

CAMPAIGN_BUILD_DIR="${BUILD_DIR}"
CAMPAIGN_RESULT_ROOT="${RESULT_ROOT}"
export BUILD_DIR RESULT_ROOT REPETITIONS END_TIME_MODE SHORT_END_TIME OUTPUT CONFIG
export ONEAPI_DEVICE_SELECTOR="${ONEAPI_DEVICE_SELECTOR:-cuda:gpu}"
export SPH_BENCH_DEVICE="${SPH_BENCH_DEVICE:-rtx-2080s}"
export SPH_BENCH_HOST="${SPH_BENCH_HOST:-$(hostname)}"
export SPH_BENCH_OS="${SPH_BENCH_OS:-$(uname -sr)}"

mkdir -p "${RESULT_ROOT}"
CAMPAIGN_LOG="${RESULT_ROOT}/campaign.log"

log() {
    printf '[%s] %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*" | tee -a "${CAMPAIGN_LOG}"
}

cleanup_and_exit() {
    local code="$1"
    log "Smoke finished with exit_code=${code}"
    log "Result root: ${RESULT_ROOT}"
    log "Returning to the shell now (no Ctrl+C needed)."
    exit "${code}"
}

trap 'cleanup_and_exit 130' INT
trap 'cleanup_and_exit 143' TERM

log "Schema smoke id: ${CAMPAIGN_ID}"
log "Build dir:       ${BUILD_DIR}"
log "Result root:     ${RESULT_ROOT}"
log "Repetitions:     ${REPETITIONS}, end_time_mode=${END_TIME_MODE}, short_end_time=${SHORT_END_TIME}, output=${OUTPUT}"
log "GPU selector:    ${ONEAPI_DEVICE_SELECTOR}, device label: ${SPH_BENCH_DEVICE}"

if [[ ! -d "${BUILD_DIR}" ]]; then
    log "ERROR: BUILD_DIR does not exist: ${BUILD_DIR}"
    cleanup_and_exit 1
fi

if [[ -f /opt/intel/oneapi/setvars.sh ]]; then
    log "Sourcing /opt/intel/oneapi/setvars.sh"
    set +euo pipefail
    # shellcheck disable=SC1091
    source /opt/intel/oneapi/setvars.sh --force >>"${CAMPAIGN_LOG}" 2>&1
    setvars_status=$?
    set -euo pipefail
    if [[ "${setvars_status}" -ne 0 ]]; then
        log "WARNING: setvars.sh returned ${setvars_status}; continuing"
    fi
else
    log "WARNING: /opt/intel/oneapi/setvars.sh not found"
fi

BUILD_DIR="${CAMPAIGN_BUILD_DIR}"
RESULT_ROOT="${CAMPAIGN_RESULT_ROOT}"
export BUILD_DIR RESULT_ROOT
unset BIN_DIR || true

source "${SCRIPT_DIR}/benchmark_run_common.sh"

failures=0
log "=== schema smoke: all verification cases (short end time) ==="
run_verification_suite "sycl" "${SPH_BENCH_DEVICE}" "${ONEAPI_DEVICE_SELECTOR}" || failures=1

summary_count="$(find "${RESULT_ROOT}" -name summary.csv 2>/dev/null | wc -l)"
log "Summaries written: ${summary_count}"
log "Inspect one header with:"
log "  head -n 1 ${RESULT_ROOT}/*/summary.csv | head -n 20"
log "Collect with:"
log "  python3 ${SCRIPT_DIR}/collect_benchmark_results.py \\"
log "    ${RESULT_ROOT} \\"
log "    ${REPO_ROOT}/tests/gpu_paper_benchmarks/curated/${CAMPAIGN_ID}_all_runs.csv"

cleanup_and_exit "${failures}"
