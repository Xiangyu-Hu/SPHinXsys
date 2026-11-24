from __future__ import annotations

import json
import sys
import time
from typing import Any


def emit(event: str, **fields: Any) -> None:
    row = {"event": event, "wall_time_s": time.time(), **fields}
    sys.stdout.write(json.dumps(row, ensure_ascii=False) + "\n")
    sys.stdout.flush()


def log(msg: str) -> None:
    sys.stderr.write(msg.rstrip() + "\n")
    sys.stderr.flush()
