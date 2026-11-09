#!/usr/bin/env bash
# Continuously refresh LIVE_STATUS.txt from active gp_* logs (for agent/UI polling).
set -euo pipefail
ROOT=/home/yongchuan/SPHinXsys-gpu-paper/docs/gpu_paper_discussion
OUT=$ROOT/logs/LIVE_STATUS.txt
LOGDIR=$ROOT/logs
INTERVAL=${1:-5}

while true; do
  {
    echo "updated: $(date -Is)"
    echo "tmux:"
    tmux ls 2>/dev/null | grep -E '^gp_' || echo "  (none)"
    echo
    for log in "$LOGDIR"/gp_A1_*.log; do
      [[ -f "$log" ]] || continue
      # skip stale finished logs unless recently written
      base=$(basename "$log")
      last=$(grep -E '^N=|^EXIT:|^START:|Passed:|Failed:|Total wall' "$log" | tail -n 3 | tr '\n' ' | ' || true)
      echo "$base: $last"
    done
  } >"$OUT.tmp"
  mv -f "$OUT.tmp" "$OUT"
  sleep "$INTERVAL"
done
