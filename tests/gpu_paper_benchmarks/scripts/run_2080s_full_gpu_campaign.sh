#!/usr/bin/env bash
# Compatibility wrapper: formal paper campaign with SPH_BENCH_DEVICE=rtx-2080s.
# Prefer run_gpu_paper_campaign.sh directly for 3090/4090 as well.
set -Eeuo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
export SPH_BENCH_DEVICE="${SPH_BENCH_DEVICE:-rtx-2080s}"
exec "${SCRIPT_DIR}/run_gpu_paper_campaign.sh" "$@"
