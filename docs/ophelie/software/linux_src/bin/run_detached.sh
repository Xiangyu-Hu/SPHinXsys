#!/usr/bin/env bash
# Detached long run helper — survives terminal disconnect better than bare foreground.
# Usage:
#   ./docs/ophelie/software/linux_src/bin/run_detached.sh \
#     --config docs/ophelie/software/linux_src/examples/mechanical-stirring/case.json \
#     --prepared /writable/path/prepared \
#     --output /writable/path/runs/f0
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOLVER="$HERE/coldcrucible-solver"
LOGDIR="${OPHELIE_DETACH_LOGDIR:-$(pwd)/ophelie_detached_logs}"
mkdir -p "$LOGDIR"
stamp="$(date -u +%Y%m%d_%H%M%S)"
out_log="$LOGDIR/detached_${stamp}.log"
pid_file="$LOGDIR/detached_${stamp}.pid"
nohup "$SOLVER" run "$@" >"$out_log" 2>&1 &
echo $! >"$pid_file"
echo "started pid=$(cat "$pid_file")"
echo "log=$out_log"
echo "stop later: touch <run-dir>/cancel.request"
