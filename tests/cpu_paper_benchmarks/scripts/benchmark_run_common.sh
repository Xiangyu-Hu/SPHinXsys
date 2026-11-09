#!/usr/bin/env bash
set -Eeuo pipefail

# Shared implementation for the CPU-paper benchmark launchers.

BENCHMARK_SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
BENCHMARK_ROOT="$(cd -- "${BENCHMARK_SCRIPT_DIR}/.." && pwd)"
REPO_ROOT="$(cd -- "${BENCHMARK_ROOT}/../.." && pwd)"

if [[ -z "${RESULT_ROOT:-}" ]]; then
    if [[ -z "${BUILD_DIR:-}" ]]; then
        printf 'ERROR: BUILD_DIR must be set before sourcing benchmark_run_common.sh\n' >&2
        exit 2
    fi
    RESULT_ROOT="${BUILD_DIR}/benchmark_results"
fi
REPETITIONS="${REPETITIONS:-3}"
OUTPUT="${OUTPUT:-off}"
END_TIME_MODE="${END_TIME_MODE:-short}"
SHORT_END_TIME="${SHORT_END_TIME:-0.01}"
CONFIG="${CONFIG-Release}"

case "${OUTPUT}" in
    on|off|true|false|1|0) ;;
    *)
        printf 'ERROR: OUTPUT must be on/off, true/false, or 1/0 (got %q)\n' "${OUTPUT}" >&2
        exit 2
        ;;
esac

case "${END_TIME_MODE}" in
    short|full) ;;
    *)
        printf 'ERROR: END_TIME_MODE must be short or full (got %q)\n' "${END_TIME_MODE}" >&2
        exit 2
        ;;
esac

if [[ ! "${REPETITIONS}" =~ ^[1-9][0-9]*$ ]]; then
    printf 'ERROR: REPETITIONS must be a positive integer (got %q)\n' "${REPETITIONS}" >&2
    exit 2
fi

mkdir -p -- "${RESULT_ROOT}"
RESULT_ROOT="$(cd -- "${RESULT_ROOT}" && pwd)"

GIT_SHORT_HASH="$(git -C "${REPO_ROOT}" rev-parse --short=12 HEAD 2>/dev/null || printf 'unknown')"
GIT_SHORT_HASH="${GIT_SHORT_HASH//$'\n'/}"

safe_slug()
{
    local value="${1,,}"
    value="${value//[^a-z0-9._-]/-}"
    while [[ "${value}" == *--* ]]; do
        value="${value//--/-}"
    done
    value="${value#-}"
    value="${value%-}"
    [[ -n "${value}" ]] || value="unknown"
    printf '%s' "${value}"
}

shell_join()
{
    local quoted
    printf -v quoted '%q ' "$@"
    printf '%s' "${quoted% }"
}

file_contains_literal()
{
    local path="$1"
    local needle="$2"
    local line
    while IFS= read -r line || [[ -n "${line}" ]]; do
        [[ "${line}" == *"${needle}"* ]] && return 0
    done <"${path}"
    return 1
}

case_full_end_time()
{
    case "$1" in
        poiseuille_2d) printf '2.0' ;;
        diffusion_neumann_2d) printf '1.0' ;;
        dambreak_3d) printf '20.0' ;;
        fish_fsi_2d) printf '1.7' ;;
        oscillating_beam_2d) printf '1.0' ;;
        twisting_column_3d) printf '0.5' ;;
        *)
            printf 'ERROR: no full end time is defined for %q\n' "$1" >&2
            return 2
            ;;
    esac
}

selected_end_time()
{
    local case_name="$1"
    if [[ -n "${END_TIME:-}" ]]; then
        printf '%s' "${END_TIME}"
    elif [[ "${END_TIME_MODE}" == "short" ]]; then
        printf '%s' "${SHORT_END_TIME}"
    elif [[ -n "${FULL_END_TIME:-}" ]]; then
        printf '%s' "${FULL_END_TIME}"
    else
        case_full_end_time "${case_name}"
    fi
}

benchmark_bin_directory()
{
    local case_dir="$1"
    # Prefer an explicit launcher override. Intel oneAPI setvars.sh exports a
    # relative BIN_DIR=bin64, which must not hijack executable lookup.
    if [[ -n "${SPH_BENCH_BIN_DIR:-}" ]]; then
        printf '%s' "${SPH_BENCH_BIN_DIR%/}"
        return
    fi
    if [[ -n "${BIN_DIR:-}" && "${BIN_DIR}" == /* && -d "${BIN_DIR}" ]]; then
        printf '%s' "${BIN_DIR%/}"
        return
    fi
    printf '%s/tests/cpu_paper_benchmarks/verification/%s/bin' \
        "${BUILD_DIR%/}" "${case_dir}"
}

benchmark_executable()
{
    local case_dir="$1"
    local target="$2"
    local bin_directory direct_candidate configured_candidate

    bin_directory="$(benchmark_bin_directory "${case_dir}")"

    direct_candidate="${bin_directory}/${target}"
    if [[ -x "${direct_candidate}" ]]; then
        printf '%s' "${direct_candidate}"
        return
    fi

    if [[ -n "${CONFIG}" ]]; then
        configured_candidate="${bin_directory}/${CONFIG}/${target}"
        if [[ -x "${configured_candidate}" ]]; then
            printf '%s' "${configured_candidate}"
            return
        fi
    fi

    printf 'ERROR: benchmark executable %q was not found. Checked %s' \
        "${target}" "${direct_candidate}" >&2
    if [[ -n "${CONFIG}" ]]; then
        printf ' and %s' "${configured_candidate}" >&2
    fi
    printf '. Set SPH_BENCH_BIN_DIR (or an absolute BIN_DIR) to the target directory, or CONFIG to the build configuration (current CONFIG=%q).\n' \
        "${CONFIG}" >&2
    printf '%s' "${direct_candidate}"
}

write_status()
{
    local run_dir="$1"
    local status="$2"
    local exit_code="$3"
    local reason="$4"
    local status_path="${run_dir}/run.status"
    if [[ -e "${status_path}" ]]; then
        printf 'ERROR: refusing to replace existing status file: %s\n' "${status_path}" >&2
        return 1
    fi
    {
        printf 'status=%s\n' "${status}"
        printf 'exit_code=%s\n' "${exit_code}"
        printf 'reason=%s\n' "${reason}"
    } >"${status_path}"
}

run_benchmark()
{
    local case_name="$1"
    local case_dir="$2"
    local target="$3"
    local resolution="$4"
    local repeat="$5"
    local backend="$6"
    local device="$7"
    local selector="$8"
    local executable end_time utc_stamp run_id run_dir temporary_log
    local command_text exit_code status reason log_path status_dir

    executable="$(benchmark_executable "${case_dir}" "${target}")"
    end_time="$(selected_end_time "${case_name}")"
    utc_stamp="$(date -u +%Y%m%dT%H%M%S%NZ)"
    run_id="$(safe_slug "${utc_stamp}-p$$-${GIT_SHORT_HASH}-${backend}-${device}-${case_name}-${resolution}-r$(printf '%02d' "${repeat}")")"
    run_dir="${RESULT_ROOT}/${case_name}/${run_id}"

    if [[ -e "${run_dir}" ]]; then
        mkdir -p -- "${RESULT_ROOT}/.launcher_failures/${case_name}"
        status_dir="$(mktemp -d \
            "${RESULT_ROOT}/.launcher_failures/${case_name}/${run_id}.preexisting.XXXXXXXX")"
        log_path="${status_dir}/run.log"
        printf 'ERROR: refusing to reuse existing run directory: %s\n' \
            "${run_dir}" >"${log_path}"
        printf 'ERROR: refusing to reuse existing run directory: %s\n' \
            "${run_dir}" >&2
        write_status "${status_dir}" "failed" "1" "run_directory_preexisting"
        return 1
    fi

    temporary_log="$(mktemp "${RESULT_ROOT}/.run-log.XXXXXXXX")"
    local -a benchmark_args=(
        "${executable}"
        "--benchmark"
        "--resolution" "${resolution}"
        "--end-time" "${end_time}"
        "--output" "${OUTPUT}"
        "--result-dir" "${RESULT_ROOT}"
        "--run-id" "${run_id}"
    )

    if [[ -n "${OUTPUT_INTERVAL:-}" ]]; then
        benchmark_args+=("--output-interval" "${OUTPUT_INTERVAL}")
    fi

    if [[ "${backend}" == "sycl" ]]; then
        command_text="$(shell_join env \
            "ONEAPI_DEVICE_SELECTOR=${selector}" \
            "SPH_BENCH_BACKEND=${backend}" \
            "SPH_BENCH_DEVICE=${device}" \
            "${benchmark_args[@]}")"
    else
        command_text="$(shell_join env -u ONEAPI_DEVICE_SELECTOR \
            "SPH_BENCH_BACKEND=${backend}" \
            "SPH_BENCH_DEVICE=${device}" \
            "${benchmark_args[@]}")"
    fi

    {
        printf 'run_id: %s\n' "${run_id}"
        printf 'utc_started: %s\n' "${utc_stamp}"
        printf 'working_directory: %s\n' "$(dirname -- "${executable}")"
        printf 'command: %s\n\n' "${command_text}"
    } >"${temporary_log}"

    exit_code=0
    if [[ ! -x "${executable}" ]]; then
        printf 'ERROR: executable is missing or not executable: %s\n' "${executable}" >>"${temporary_log}"
        exit_code=127
        reason="executable_missing"
    elif [[ "${backend}" == "sycl" ]]; then
        (
            cd -- "$(dirname -- "${executable}")"
            env "ONEAPI_DEVICE_SELECTOR=${selector}" \
                "SPH_BENCH_BACKEND=${backend}" \
                "SPH_BENCH_DEVICE=${device}" \
                "${benchmark_args[@]}"
        ) >>"${temporary_log}" 2>&1 || exit_code=$?
        reason="process_exit"
    else
        (
            cd -- "$(dirname -- "${executable}")"
            env -u ONEAPI_DEVICE_SELECTOR \
                "SPH_BENCH_BACKEND=${backend}" \
                "SPH_BENCH_DEVICE=${device}" \
                "${benchmark_args[@]}"
        ) >>"${temporary_log}" 2>&1 || exit_code=$?
        reason="process_exit"
    fi

    status_dir="${run_dir}"
    if [[ "${exit_code}" -ne 0 ]] &&
        file_contains_literal "${temporary_log}" \
            "Refusing to overwrite benchmark run"; then
        mkdir -p -- "${RESULT_ROOT}/.launcher_failures/${case_name}"
        status_dir="$(mktemp -d \
            "${RESULT_ROOT}/.launcher_failures/${case_name}/${run_id}.conflict.XXXXXXXX")"
        reason="run_directory_conflict"
        printf '\nERROR: the benchmark refused a concurrently created run directory; artifacts are isolated here: %s\n' \
            "${status_dir}" >>"${temporary_log}"
    elif [[ ! -d "${run_dir}" ]]; then
        mkdir -p -- "${RESULT_ROOT}/.launcher_failures/${case_name}"
        status_dir="$(mktemp -d \
            "${RESULT_ROOT}/.launcher_failures/${case_name}/${run_id}.failed.XXXXXXXX")"
        printf '\nERROR: the benchmark did not create its run directory; failure artifacts are isolated here: %s\n' \
            "${status_dir}" >>"${temporary_log}"
    fi

    if [[ "${exit_code}" -eq 0 && ! -f "${status_dir}/summary.csv" ]]; then
        exit_code=1
        reason="summary_missing"
        printf '\nERROR: process exited successfully but did not create summary.csv\n' >>"${temporary_log}"
    fi

    if [[ "${exit_code}" -eq 0 ]]; then
        status="completed"
        reason="process_exit"
    else
        status="failed"
    fi

    log_path="${status_dir}/run.log"
    if [[ -e "${log_path}" ]]; then
        printf 'ERROR: refusing to replace existing log file: %s\n' "${log_path}" >&2
        rm -f -- "${temporary_log}"
        return 1
    fi
    mv -- "${temporary_log}" "${log_path}"
    write_status "${status_dir}" "${status}" "${exit_code}" "${reason}"

    printf '%s: %s (%s, repeat %d)\n' \
        "${status^^}" "${case_name}" "${resolution}" "${repeat}"
    return "${exit_code}"
}

run_verification_suite()
{
    local backend="$1"
    local device="$2"
    local selector="$3"
    local failures=0
    local repeat specification case_name case_dir target resolutions resolution
    local -a verification_cases=(
        "poiseuille_2d|poiseuille_2d|cpu_paper_poiseuille_2d|coarse standard fine"
        "diffusion_neumann_2d|diffusion_neumann_2d|cpu_paper_diffusion_neumann_2d|coarse standard fine"
        "dambreak_3d|dambreak_3d|cpu_paper_dambreak_3d|standard"
        "fish_fsi_2d|fish_fsi_2d|cpu_paper_fish_fsi_2d|standard"
        "oscillating_beam_2d|oscillating_beam_2d|cpu_paper_oscillating_beam_2d|standard"
        "twisting_column_3d|twisting_column_3d|cpu_paper_twisting_column_3d|standard"
    )

    for ((repeat = 1; repeat <= REPETITIONS; repeat += 1)); do
        for specification in "${verification_cases[@]}"; do
            IFS='|' read -r case_name case_dir target resolutions <<<"${specification}"
            for resolution in ${resolutions}; do
                run_benchmark "${case_name}" "${case_dir}" "${target}" \
                    "${resolution}" "${repeat}" "${backend}" "${device}" \
                    "${selector}" || failures=1
            done
        done
    done
    return "${failures}"
}

run_dambreak_scaling_suite()
{
    local backend="$1"
    local device="$2"
    local selector="$3"
    local failures=0
    local repeat resolution

    for ((repeat = 1; repeat <= REPETITIONS; repeat += 1)); do
        for resolution in s1 s2 s3 s4 s5 s6; do
            run_benchmark "dambreak_3d" "dambreak_3d" \
                "cpu_paper_dambreak_3d" "${resolution}" "${repeat}" \
                "${backend}" "${device}" "${selector}" || failures=1
        done
    done
    return "${failures}"
}
