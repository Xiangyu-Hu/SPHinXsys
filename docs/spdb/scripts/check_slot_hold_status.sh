#!/usr/bin/env bash
# Inspect a background slot_hold job. Safe after agent disconnect.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
OUT="${1:-${ROOT}/docs/spdb/results/slot_hold_dt60_n120}"

echo "=== slot_hold status ==="
echo "out=${OUT}"

if [[ ! -d "${OUT}" ]]; then
  echo "missing directory" >&2
  exit 1
fi

PID=""
ALIVE=0
if [[ -f "${OUT}/job.pid" ]]; then
  PID="$(cat "${OUT}/job.pid")"
  if kill -0 "${PID}" 2>/dev/null; then
    ALIVE=1
  fi
fi

echo "pid=${PID:-none} alive=${ALIVE}"

if [[ -f "${OUT}/job_meta.txt" ]]; then
  echo "--- job_meta ---"
  cat "${OUT}/job_meta.txt"
fi

if [[ -f "${OUT}/hold_status.txt" ]]; then
  echo "--- hold_status.txt ---"
  cat "${OUT}/hold_status.txt"
else
  echo "hold_status.txt: not written yet (still configuring / early elastic)"
fi

if [[ -f "${OUT}/hold_heartbeat" ]]; then
  echo "--- hold_heartbeat ---"
  cat "${OUT}/hold_heartbeat"
  HB_AGE=$(( $(date +%s) - $(stat -c %Y "${OUT}/hold_heartbeat") ))
  echo "heartbeat_age_s=${HB_AGE}"
  if [[ "${ALIVE}" -eq 1 && "${HB_AGE}" -gt 900 ]]; then
    echo "WARN: process alive but heartbeat older than 15 min"
  fi
fi

if [[ -f "${OUT}/run.log" ]]; then
  echo "--- run.log tail ---"
  tail -n 8 "${OUT}/run.log"
  LAST_ACC="$(grep -E 'accepted hold_physical_s=' "${OUT}/run.log" | tail -n 1 || true)"
  if [[ -n "${LAST_ACC}" ]]; then
    echo "last_accepted_line=${LAST_ACC}"
  fi
fi

METRICS="$(ls -1 "${OUT}"/spdb_slot_hold_*_metrics.txt 2>/dev/null | tail -n 1 || true)"
if [[ -n "${METRICS}" ]]; then
  echo "--- metrics present ---"
  echo "${METRICS}"
  grep -E 'hold_physical_s=|aborted=|wall_time_s=|hold_is_7200s=' "${METRICS}" || true
fi

if [[ "${ALIVE}" -eq 0 ]]; then
  if [[ -n "${METRICS}" ]]; then
    echo "STATE=finished_or_stopped_with_metrics"
  elif grep -q 'accepted hold_physical_s=' "${OUT}/run.log" 2>/dev/null; then
    echo "STATE=stopped_incomplete"
  else
    echo "STATE=not_running"
  fi
else
  echo "STATE=running"
fi
