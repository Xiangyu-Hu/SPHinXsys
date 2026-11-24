from __future__ import annotations

import json
from pathlib import Path
from typing import Any


def load_case(path: Path) -> dict[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError("case.json must be an object")
    data["_case_path"] = str(path.resolve())
    data["_case_dir"] = str(path.resolve().parent)
    return data


def asset_path(case: dict[str, Any], rel: str) -> Path:
    return (Path(case["_case_dir"]) / rel).resolve()


def stirring_enabled(case: dict[str, Any]) -> bool:
    stirring = case.get("stirring") or {}
    if "enabled" in stirring:
        return bool(stirring["enabled"])
    return bool(case.get("stirring_enabled", False))


def mechanics_mode(case: dict[str, Any]) -> str:
    stirring = case.get("stirring") or {}
    if "mechanics_mode" in stirring:
        return str(stirring["mechanics_mode"])
    return str(case.get("mechanics_mode", "none"))
