#!/usr/bin/env bash
set -euo pipefail
PACK=$(cd "$(dirname "$0")/.." && pwd)
ROOT=$(cd "$PACK/../../../.." && pwd)
cd "$ROOT"
failed=0
for p in "$PACK"/patches/nc/*.patch "$PACK"/patches/stir/*.patch; do
  [[ -f "$p" ]] || continue
  echo "[pack] apply $(basename "$(dirname "$p")")/$(basename "$p")"
  if git apply --check "$p" 2>/dev/null; then
    git apply "$p"
  elif git apply --3way "$p"; then
    echo "[pack] 3-way ok"
  else
    echo "[pack] FAIL $p" >&2
    failed=1
  fi
done
exit "$failed"
