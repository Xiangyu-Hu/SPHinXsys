#!/usr/bin/env bash
# Launch French NC / stirring in tmux (single GPU → prefer sequential).
#
# Usage:
#   ./RUN_TMUX_FR_TWO_CASES.sh              # default: one session, NC then STIR
#   ./RUN_TMUX_FR_TWO_CASES.sh both-seq     # same as default
#   ./RUN_TMUX_FR_TWO_CASES.sh nc           # only NC  → session doc247_fr_nc
#   ./RUN_TMUX_FR_TWO_CASES.sh stir         # only STIR → session doc247_fr_stir
#
# Attach / list:
#   tmux attach -t doc247_fr_seq
#   tmux attach -t doc247_fr_nc
#   tmux ls
set -euo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)
NC="$HERE/01_fr_nc/run_fr_nc.sh"
STIR="$HERE/02_fr_stir/run_fr_stir.sh"
MODE=${1:-both-seq}

need() { [[ -x "$1" ]] || { echo "FATAL: not executable $1"; exit 1; }; }
need "$NC"
need "$STIR"
command -v tmux >/dev/null || { echo "FATAL: tmux not found"; exit 1; }

start_session() {
  local name=$1 cmd=$2
  if tmux has-session -t "$name" 2>/dev/null; then
    echo "WARN: session '$name' already exists — skip. Attach: tmux attach -t $name"
    return 0
  fi
  tmux new-session -d -s "$name" -n run "bash -lc '$cmd; echo; echo DONE exit=\$?; exec bash'"
  echo "OK: started tmux session '$name'"
  echo "    attach: tmux attach -t $name"
}

case "$MODE" in
  nc)
    start_session doc247_fr_nc "$NC"
    ;;
  stir)
    start_session doc247_fr_stir "$STIR"
    ;;
  both-seq|"")
    start_session doc247_fr_seq "$NC && $STIR"
    ;;
  *)
    echo "Usage: $0 [both-seq|nc|stir]"
    exit 2
    ;;
esac
