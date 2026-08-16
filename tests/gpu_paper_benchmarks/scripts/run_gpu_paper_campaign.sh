#!/usr/bin/env bash
# Formal GPU paper campaign (CUDA SYCL): data intended for manuscript tables/figures.
#
# Coverage (REPETITIONS default 3, END_TIME_MODE=full, OUTPUT=off):
#   Phase 1 — verification
#     poiseuille_2d:        coarse, standard, fine
#     diffusion_neumann_2d: coarse, standard, fine
#     dambreak_3d:          standard
#     fish_fsi_2d:          standard
#     oscillating_beam_2d:  standard
#     twisting_column_3d:   standard
#   Phase 2 — dam-break scaling
#     dambreak_3d:          s1 .. s6  (end_time=20)
#
# Total with defaults: 10 verification configs × 3 + 6 scaling × 3 = 48 runs.
#
# Foreground (recommended; returns to the shell when finished):
#   SPH_BENCH_DEVICE=rtx-2080s \
#     ./tests/gpu_paper_benchmarks/scripts/run_gpu_paper_campaign.sh
#
# Same protocol on other GPUs (only change the device label / campaign id):
#   SPH_BENCH_DEVICE=rtx-3090 ./tests/gpu_paper_benchmarks/scripts/run_gpu_paper_campaign.sh
#   SPH_BENCH_DEVICE=rtx-4090 ./tests/gpu_paper_benchmarks/scripts/run_gpu_paper_campaign.sh
#
set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../../.." && pwd)"
BUILD_DIR="${BUILD_DIR:-${REPO_ROOT}/build-sycl}"

SPH_BENCH_DEVICE="${SPH_BENCH_DEVICE:-rtx-2080s}"
device_slug="$(printf '%s' "${SPH_BENCH_DEVICE}" | tr '[:upper:]' '[:lower:]' | sed -E 's/[^a-z0-9._-]+/-/g; s/^-+//; s/-+$//; s/--+/-/g')"
[[ -n "${device_slug}" ]] || device_slug="gpu"

CAMPAIGN_ID="${CAMPAIGN_ID:-${device_slug}-paper-$(date -u +%Y%m%dT%H%M%SZ)}"
RESULT_ROOT="${RESULT_ROOT:-${BUILD_DIR}/benchmark_results/${CAMPAIGN_ID}}"
REPETITIONS="${REPETITIONS:-3}"
END_TIME_MODE="${END_TIME_MODE:-full}"
OUTPUT="${OUTPUT:-off}"
CONFIG="${CONFIG:-Release}"
WARMUP="${WARMUP:-1}"
COLLECT="${COLLECT:-1}"

CAMPAIGN_BUILD_DIR="${BUILD_DIR}"
CAMPAIGN_RESULT_ROOT="${RESULT_ROOT}"
export BUILD_DIR RESULT_ROOT REPETITIONS END_TIME_MODE OUTPUT CONFIG
export ONEAPI_DEVICE_SELECTOR="${ONEAPI_DEVICE_SELECTOR:-cuda:gpu}"
export SPH_BENCH_DEVICE
export SPH_BENCH_HOST="${SPH_BENCH_HOST:-$(hostname)}"
export SPH_BENCH_OS="${SPH_BENCH_OS:-$(uname -sr)}"

mkdir -p "${RESULT_ROOT}"
CAMPAIGN_LOG="${RESULT_ROOT}/campaign.log"
CURATED_ALL="${REPO_ROOT}/tests/gpu_paper_benchmarks/curated/${CAMPAIGN_ID}_all_runs.csv"
CURATED_STATS="${REPO_ROOT}/tests/gpu_paper_benchmarks/curated/${CAMPAIGN_ID}_repeat_stats.csv"

log() {
    printf '[%s] %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*" | tee -a "${CAMPAIGN_LOG}"
}

finish() {
    local code="$1"
    log "Campaign finished with exit_code=${code}"
    log "Result root: ${RESULT_ROOT}"
    log "Returning to the shell now (no Ctrl+C needed)."
    exit "${code}"
}

trap 'finish 130' INT
trap 'finish 143' TERM

log "Campaign id: ${CAMPAIGN_ID}"
log "Build dir:   ${BUILD_DIR}"
log "Result root: ${RESULT_ROOT}"
log "Repetitions: ${REPETITIONS}, end_time_mode=${END_TIME_MODE}, output=${OUTPUT}"
log "GPU selector: ${ONEAPI_DEVICE_SELECTOR}, device label: ${SPH_BENCH_DEVICE}"
log "Expected completed summaries with defaults: $((10 * REPETITIONS + 6 * REPETITIONS))"

if [[ ! -d "${BUILD_DIR}" ]]; then
    log "ERROR: BUILD_DIR does not exist: ${BUILD_DIR}"
    finish 1
fi

if [[ -f /opt/intel/oneapi/setvars.sh ]]; then
    log "Sourcing /opt/intel/oneapi/setvars.sh"
    set +euo pipefail
    # shellcheck disable=SC1091
    source /opt/intel/oneapi/setvars.sh --force >>"${CAMPAIGN_LOG}" 2>&1
    setvars_status=$?
    set -euo pipefail
    if [[ "${setvars_status}" -ne 0 ]]; then
        log "WARNING: setvars.sh returned ${setvars_status}; continuing with current environment"
    else
        log "setvars.sh completed"
    fi
else
    log "WARNING: /opt/intel/oneapi/setvars.sh not found"
fi

BUILD_DIR="${CAMPAIGN_BUILD_DIR}"
RESULT_ROOT="${CAMPAIGN_RESULT_ROOT}"
export BUILD_DIR RESULT_ROOT
unset BIN_DIR || true

source "${SCRIPT_DIR}/benchmark_run_common.sh"
log "Shared launcher loaded; starting warm-up / suites"
log "Resolved dambreak binary: $(benchmark_executable dambreak_3d gpu_paper_dambreak_3d)"

GPU_SELECTOR="${ONEAPI_DEVICE_SELECTOR}"
GPU_DEVICE="${SPH_BENCH_DEVICE}"
failures=0

if [[ "${WARMUP}" == "1" ]]; then
    log "Warm-up: dambreak_3d standard (end_time=${SHORT_END_TIME}, not counted in campaign totals)"
    warmup_dir="${BUILD_DIR}/tests/gpu_paper_benchmarks/verification/dambreak_3d/bin"
    warmup_bin="${warmup_dir}/gpu_paper_dambreak_3d"
    if [[ ! -x "${warmup_bin}" ]]; then
        log "WARNING: warm-up binary missing: ${warmup_bin}"
    else
        (
            cd "${warmup_dir}"
            env "ONEAPI_DEVICE_SELECTOR=${GPU_SELECTOR}" \
                "SPH_BENCH_BACKEND=sycl" \
                "SPH_BENCH_DEVICE=${GPU_DEVICE}" \
                "${warmup_bin}" --benchmark --resolution=standard \
                --end-time="${SHORT_END_TIME}" --output=off \
                --result-dir="${BUILD_DIR}/benchmark_results/_warmup" \
                --run-id="warmup-$(date -u +%Y%m%dT%H%M%SZ)"
        ) >>"${CAMPAIGN_LOG}" 2>&1 || log "WARNING: warm-up exited non-zero (continuing)"
    fi
fi

log "=== Phase 1/2: verification suite ==="
run_verification_suite "sycl" "${GPU_DEVICE}" "${GPU_SELECTOR}" || failures=1

log "=== Phase 2/2: dam-break scaling s1..s6 ==="
run_dambreak_scaling_suite "sycl" "${GPU_DEVICE}" "${GPU_SELECTOR}" || failures=1

summary_count="$(find "${RESULT_ROOT}" -name summary.csv 2>/dev/null | wc -l)"
status_count="$(find "${RESULT_ROOT}" -name run.status 2>/dev/null | wc -l)"
completed_count="$(find "${RESULT_ROOT}" -name run.status -exec grep -l '^status=completed$' {} + 2>/dev/null | wc -l)"
log "Summaries written: ${summary_count}, run.status files: ${status_count}, completed: ${completed_count}"

if [[ "${COLLECT}" == "1" ]]; then
    log "Collecting curated CSVs"
    if python3 "${SCRIPT_DIR}/collect_benchmark_results.py" \
        "${RESULT_ROOT}" \
        "${CURATED_ALL}" \
        --stats-output "${CURATED_STATS}" \
        --allow-failures >>"${CAMPAIGN_LOG}" 2>&1; then
        log "Wrote ${CURATED_ALL}"
        log "Wrote ${CURATED_STATS}"
    else
        log "WARNING: collect_benchmark_results.py failed; raw results remain under ${RESULT_ROOT}"
        failures=1
    fi
else
    log "Collect skipped (COLLECT=${COLLECT}). Example:"
    log "  python3 ${SCRIPT_DIR}/collect_benchmark_results.py \\"
    log "    ${RESULT_ROOT} \\"
    log "    ${CURATED_ALL} \\"
    log "    --stats-output ${CURATED_STATS}"
fi

finish "${failures}"
