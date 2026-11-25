from __future__ import annotations

from pathlib import Path
from typing import Any

from . import SCHEMA_VERSION
from .case_io import asset_path, mechanics_mode, stirring_enabled
from .paths import resolve_binary


def _issue(path: str, severity: str, code: str, message: str) -> dict[str, str]:
    return {"path": path, "severity": severity, "code": code, "message": message}


_ALLOWED_HEAT_MODES = {
    "none",
    "fixed_absorbed_power",
    "prescribed_volume_heat",
    "f1",
    "online_em",
    "given_current",
}

_ALLOWED_NC_PROFILES = {
    "isothermal_hydro",
    "dT20_boussinesq",
    "fixed_absorbed_power_em",
}


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
        dp = (case.get("numerics") or {}).get("dp_m")
        if dp is not None and abs(float(dp) - 0.005) > 1e-9:
            issues.append(
                _issue(
                    "numerics.dp_m",
                    "warning",
                    "DP_NOT_DEFAULT_DENSE",
                    f"stirring engineering default is dp=0.005; case has dp={dp}",
                )
            )
    else:
        if mode not in ("none", "", None):
            issues.append(_issue("stirring.mechanics_mode", "warning", "MECHANICS_IGNORED", "stirring off; mechanics_mode ignored for solve"))

    heating = case.get("heating") or {}
    heat_mode = heating.get("mode", "none")
    if heat_mode not in _ALLOWED_HEAT_MODES:
        issues.append(
            _issue(
                "heating.mode",
                "error",
                "ILLEGAL_HEAT_MODE",
                f"unsupported heating.mode={heat_mode!r}; allowed={sorted(_ALLOWED_HEAT_MODES)}",
            )
        )

    power = heating.get("generator_input_power_W") or ((case.get("power") or {}).get("generator_input_power_W"))
    if power is not None and float(power) == 30000:
        issues.append(
            _issue(
                "heating.generator_input_power_W",
                "warning",
                "POWER_SCOPE",
                "30 kW is generator input, not melt absorbed power / not verified EM",
            )
        )

    current = heating.get("coil_current_A")
    if current is None:
        current = (case.get("power") or {}).get("coil_current_A")
    current_unknown = current in (None, 0, "unknown", "")

    if heat_mode in ("online_em", "given_current") and current_unknown:
        issues.append(
            _issue(
                "heating.coil_current_A",
                "error",
                "CURRENT_UNKNOWN",
                "current unknown; cannot enable given-current / online_em",
            )
        )

    if heat_mode == "fixed_absorbed_power":
        target = heating.get("target_absorbed_power_W")
        if target is None or float(target) <= 0:
            issues.append(
                _issue(
                    "heating.target_absorbed_power_W",
                    "error",
                    "ABSORBED_TARGET",
                    "fixed_absorbed_power requires target_absorbed_power_W > 0",
                )
            )
        if heating.get("power_status") not in ("ASSUMED", "measured", "calibrated"):
            issues.append(
                _issue(
                    "heating.power_status",
                    "error",
                    "POWER_STATUS",
                    "fixed_absorbed_power requires power_status in {ASSUMED,measured,calibrated}",
                )
            )

    run_cfg = case.get("run") or {}
    profile = run_cfg.get("profile")
    if not stir_on:
        if profile not in _ALLOWED_NC_PROFILES:
            issues.append(
                _issue(
                    "run.profile",
                    "error",
                    "ILLEGAL_NC_PROFILE",
                    f"unsupported NC run.profile={profile!r}; allowed={sorted(_ALLOWED_NC_PROFILES)}",
                )
            )
        # Conflict: real EM heating must not use zero-q diagnostic profiles
        if heat_mode == "fixed_absorbed_power" and profile in ("isothermal_hydro", "dT20_boussinesq"):
            issues.append(
                _issue(
                    "run.profile",
                    "error",
                    "HEAT_PROFILE_CONFLICT",
                    "fixed_absorbed_power cannot use isothermal_hydro/dT20_boussinesq (those inject --zero-q)",
                )
            )
        if profile == "fixed_absorbed_power_em" and heat_mode != "fixed_absorbed_power":
            issues.append(
                _issue(
                    "heating.mode",
                    "error",
                    "HEAT_PROFILE_CONFLICT",
                    "run.profile=fixed_absorbed_power_em requires heating.mode=fixed_absorbed_power",
                )
            )
        if heat_mode == "none" and profile == "fixed_absorbed_power_em":
            issues.append(
                _issue(
                    "heating.mode",
                    "error",
                    "HEAT_PROFILE_CONFLICT",
                    "EM profile refuses heating.mode=none (no silent isothermal fallback)",
                )
            )

    # binaries present?
    need = ["nc_relax", "nc_run"] if not stir_on else ["stir_relax", "stir_f0"]
    if stir_on and heat_mode in ("prescribed_volume_heat", "f1"):
        need.append("stir_f1")
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
        "heating_mode": heat_mode,
        "issues": issues,
    }
