from __future__ import annotations

from pathlib import Path
from typing import Any

from . import SCHEMA_VERSION
from .case_io import asset_path, mechanics_mode, stirring_enabled
from .paths import resolve_binary


def _issue(path: str, severity: str, code: str, message: str) -> dict[str, str]:
    return {"path": path, "severity": severity, "code": code, "message": message}


def validate_case(case: dict[str, Any]) -> dict[str, Any]:
    issues: list[dict[str, str]] = []

    schema = case.get("schema_version")
    if schema not in (SCHEMA_VERSION, "client-case-2026-11-20", "ophelie-case-2026-11-21"):
        issues.append(_issue("schema_version", "error", "SCHEMA", f"unsupported schema_version={schema!r}"))

    assembly = case.get("assembly") or {}
    roles = assembly.get("roles") or {}
    required = ["crucible", "coil", "cooling_pipes", "thermocouple"]
    for role in required:
        rel = roles.get(role)
        if not rel:
            issues.append(_issue(f"assembly.roles.{role}", "error", "MISSING_ROLE", f"missing STL role: {role}"))
            continue
        p = asset_path(case, rel)
        if not p.is_file():
            issues.append(_issue(f"assembly.roles.{role}", "error", "MISSING_STL", f"STL not found: {p}"))

    stir_on = stirring_enabled(case)
    mode = mechanics_mode(case)
    if stir_on:
        paddle = roles.get("paddle")
        if not paddle:
            issues.append(_issue("assembly.roles.paddle", "error", "MISSING_PADDLE", "stirring.enabled requires paddle STL"))
        else:
            p = asset_path(case, paddle)
            if not p.is_file():
                issues.append(_issue("assembly.roles.paddle", "error", "MISSING_STL", f"paddle STL not found: {p}"))
        if mode != "elastic_driven":
            issues.append(
                _issue(
                    "stirring.mechanics_mode",
                    "error",
                    "MECHANICS",
                    f"stirring project requires mechanics_mode=elastic_driven, got {mode!r}; rigid fallback forbidden",
                )
            )
        solver_paddle = roles.get("paddle_solver")
        if solver_paddle and not asset_path(case, solver_paddle).is_file():
            issues.append(_issue("assembly.roles.paddle_solver", "error", "MISSING_STL", "paddle_solver STL missing"))
    else:
        if mode not in ("none", "", None):
            issues.append(_issue("stirring.mechanics_mode", "warning", "MECHANICS_IGNORED", "stirring off; mechanics_mode ignored for solve"))

    power = ((case.get("heating") or {}).get("generator_input_power_W")
             or ((case.get("power") or {}).get("generator_input_power_W")))
    if power is not None and float(power) == 30000:
        issues.append(
            _issue(
                "heating.generator_input_power_W",
                "warning",
                "POWER_SCOPE",
                "30 kW is generator input, not melt absorbed power / not verified EM",
            )
        )

    current = ((case.get("heating") or {}).get("coil_current_A")
               or ((case.get("power") or {}).get("coil_current_A")))
    heat_mode = (case.get("heating") or {}).get("mode", "none")
    if heat_mode in ("online_em", "given_current") and (current in (None, 0, "unknown")):
        issues.append(_issue("heating.coil_current_A", "error", "CURRENT_UNKNOWN", "current unknown; cannot enable given-current EM"))

    # binaries present?
    need = ["nc_relax", "nc_run"] if not stir_on else ["stir_relax", "stir_f0"]
    for key in need:
        b = resolve_binary(key)
        if not b.is_file():
            issues.append(_issue("backend.binaries", "error", "MISSING_BIN", f"solver binary missing: {b}"))

    errors = [i for i in issues if i["severity"] == "error"]
    return {
        "ok": len(errors) == 0,
        "schema_version": SCHEMA_VERSION,
        "stirring_enabled": stir_on,
        "mechanics_mode": mode if stir_on else "none",
        "issues": issues,
    }
