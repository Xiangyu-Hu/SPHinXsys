#!/usr/bin/env bash
# Launch slot_hold under nohup so it survives agent / SSH disconnect.
# Writes: run.log, job.pid, hold_status.txt (from the binary), hold_heartbeat.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
BIN="${ROOT}/build-cpu-icpx/tests/extra_source_and_tests/3d_examples/test_3d_spdb_geometry/bin/test_3d_spdb_geometry"
DP="${DP:-0.001}"
PHYSICAL_DT="${PHYSICAL_DT:-60}"
HOLD_STEPS="${HOLD_STEPS:-120}"
STIFFNESS="${STIFFNESS:-1}"
END_TIME="${END_TIME:-0.0015}"
OUT="${OUT:-${ROOT}/docs/spdb/results/slot_hold_dt${PHYSICAL_DT}_n${HOLD_STEPS}}"

if [[ ! -x "${BIN}" ]]; then
  echo "Missing binary: ${BIN}" >&2
  exit 1
fi

mkdir -p "${OUT}"
cd "${OUT}"

if [[ -f job.pid ]] && kill -0 "$(cat job.pid)" 2>/dev/null; then
  echo "Already running: pid=$(cat job.pid) in ${OUT}" >&2
  echo "Check with: ${ROOT}/docs/spdb/scripts/check_slot_hold_status.sh ${OUT}" >&2
  exit 2
fi

# Preserve a previous incomplete log once.
if [[ -f run.log ]]; then
  mv -f run.log "run.log.bak.$(date +%Y%m%d_%H%M%S)"
fi

TARGET_S="$(awk -v dt="${PHYSICAL_DT}" -v n="${HOLD_STEPS}" 'BEGIN{printf "%.0f", dt*n}')"

cat > job_meta.txt <<EOF
started_at=$(date -Is)
out=${OUT}
dp=${DP}
physical_dt=${PHYSICAL_DT}
hold_steps=${HOLD_STEPS}
target_hold_s=${TARGET_S}
bin=${BIN}
note=NOT_claimed_complete_until_finished=1_in_hold_status.txt
EOF

nohup stdbuf -oL -eL "${BIN}" \
  --slot_hold \
  --dp "${DP}" \
  --physical_dt "${PHYSICAL_DT}" \
  --hold_steps "${HOLD_STEPS}" \
  --contact_stiffness_scale "${STIFFNESS}" \
  --end_time "${END_TIME}" \
  > run.log 2>&1 &

echo $! > job.pid
echo "started pid=$(cat job.pid)"
echo "out=${OUT}"
echo "target_hold_s=${TARGET_S}"
echo "monitor: ${ROOT}/docs/spdb/scripts/check_slot_hold_status.sh ${OUT}"
echo "or:     cat ${OUT}/hold_status.txt"
echo "or:     tail -f ${OUT}/run.log"
