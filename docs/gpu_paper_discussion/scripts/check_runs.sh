#!/usr/bin/env bash
# 重连后快速查看 GPU 论文跑点状态。不启动新任务。
set -euo pipefail

ROOT="/home/yongchuan/SPHinXsys-gpu-paper/docs/gpu_paper_discussion"
LOGDIR="$ROOT/logs"
SESSDIR="$ROOT/sessions"

echo "=== tmux sessions (gp_*) ==="
if tmux ls 2>/dev/null | grep -E '^gp_' ; then
  :
else
  echo "(none)"
fi

echo
echo "=== GPU 1 ==="
nvidia-smi -i 1 --query-gpu=name,uuid,utilization.gpu,memory.used,memory.total --format=csv,noheader 2>/dev/null || echo "nvidia-smi failed"

echo
echo "=== recent logs (mtime) ==="
if [[ -d "$LOGDIR" ]]; then
  find "$LOGDIR" -type f -name '*.log' -printf '%T@ %p\n' 2>/dev/null \
    | sort -nr | head -n 10 | while read -r _ path; do
        echo "--- $path (last 5 lines) ---"
        tail -n 5 "$path" || true
        echo
      done
else
  echo "(no $LOGDIR yet)"
fi

echo "=== sessions meta ==="
if [[ -d "$SESSDIR" ]]; then
  ls -1 "$SESSDIR" 2>/dev/null || true
else
  echo "(no $SESSDIR yet)"
fi
