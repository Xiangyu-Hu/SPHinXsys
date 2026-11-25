from __future__ import annotations

import hashlib
import json
import math
import shutil
import subprocess
from pathlib import Path
from typing import Any

from .case_io import asset_path, mechanics_mode, stirring_enabled
from .events import emit, log
from .paths import resolve_binary
from .validate import validate_case


def _file_sha256(path: Path) -> str | None:
    if not path.is_file():
        return None
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _fingerprint(case: dict[str, Any]) -> str:
    roles = (case.get("assembly") or {}).get("roles") or {}
    role_hashes = {}
    for role, rel in roles.items():
        if not rel:
            continue
        p = asset_path(case, rel)
        role_hashes[role] = {"path": rel, "sha256": _file_sha256(p)}
    payload = {
        "stirring": stirring_enabled(case),
        "mechanics": mechanics_mode(case),
        "dp": (case.get("numerics") or {}).get("dp_m"),
        "melt_height": (case.get("geometry") or {}).get("melt_height_m"),
        "roles": roles,
        "role_content_sha256": role_hashes,
        "heating_mode": (case.get("heating") or {}).get("mode"),
        "run_profile": (case.get("run") or {}).get("profile"),
    }
    raw = json.dumps(payload, sort_keys=True).encode()
    return hashlib.sha256(raw).hexdigest()[:16]


def _run(cmd: list[str], cwd: Path) -> int:
    log(f"$ {' '.join(cmd)}  (cwd={cwd})")
    proc = subprocess.run(cmd, cwd=str(cwd))
    return proc.returncode


def _audit_reload(reload_dir: Path, expected_dp: float | None) -> dict[str, Any]:
    """Best-effort volume / neighbor-scale check from Reload.xml particle volumes."""
    audit: dict[str, Any] = {"reload_dir": str(reload_dir), "bodies": {}}
    xml = reload_dir / "Reload.xml"
    if not xml.is_file():
        audit["error"] = "Reload.xml missing"
        return audit
    text = xml.read_text(errors="replace")
    # crude Vol_ extract if present
    import re

    vols = [float(x) for x in re.findall(r"<Vol[^>]*>([^<]+)</Vol>", text)]
    if not vols:
        # SPHINXSYS sometimes stores as attribute-less text blocks; keep empty audit
        audit["note"] = "no <Vol> tags parsed; skip Vol^(1/3) check"
        return audit
    mean_v = sum(vols) / len(vols)
    dp_est = mean_v ** (1.0 / 3.0) if mean_v > 0 else None
    audit["n_vol"] = len(vols)
    audit["mean_vol"] = mean_v
    audit["dp_est_vol_cbrt"] = dp_est
    if expected_dp is not None and dp_est is not None:
        audit["dp_expected"] = expected_dp
        audit["dp_rel_err"] = abs(dp_est - expected_dp) / expected_dp
        audit["dp_ok"] = audit["dp_rel_err"] < 0.25  # coarse gate; not a substitute for NN audit
    return audit


def prepare(case: dict[str, Any], output_dir: Path) -> int:
    report = validate_case(case)
    if not report["ok"]:
        emit("error", stage="validate", issues=report["issues"])
        return 2

    fp = _fingerprint(case)
    prepared = output_dir / fp
    prepared.mkdir(parents=True, exist_ok=True)
    emit("start", stage="prepare", fingerprint=fp, output=str(prepared))

    stir = stirring_enabled(case)
    numerics = case.get("numerics") or {}
    geom = case.get("geometry") or {}
    dp = float(numerics.get("dp_m", 0.010 if not stir else 0.005))
    melt_h = float(geom.get("melt_height_m", 0.180))

    if not stir:
        bin_path = resolve_binary("nc_relax")
        work = bin_path.parent
        rc = _run(
            [
                str(bin_path),
                f"--dp={dp}",
                f"--melt-height={melt_h}",
                "--relax-steps=1000",
                "--relax-vtp-every=100",
                "--output-meta-dir=./client_nc_relax_meta",
            ],
            work,
        )
        if rc != 0:
            emit("error", stage="prepare", code=rc, message="nc_relax failed")
            return rc
        reload_src = work / "reload"
    else:
        bin_path = resolve_binary("stir_relax")
        work = bin_path.parent
        # Must actually pass dp to generator (not bodies-only).
        rc = _run(
            [
                str(bin_path),
                "--bodies=all",
                f"--dp={dp}",
                f"--melt-height={melt_h}",
            ],
            work,
        )
        if rc != 0:
            emit("error", stage="prepare", code=rc, message="stir_relax failed")
            return rc
        reload_src = work / "reload"

    if not reload_src.is_dir():
        emit("error", stage="prepare", message=f"reload missing after generate: {reload_src}")
        return 2

    reload_dst = prepared / "reload"
    if reload_dst.exists():
        shutil.rmtree(reload_dst)
    shutil.copytree(reload_src, reload_dst)

    audit = _audit_reload(reload_dst, expected_dp=dp)
    if audit.get("dp_ok") is False:
        emit(
            "error",
            stage="prepare",
            message=f"Reload Vol^(1/3) dp mismatch: {audit}",
        )
        return 3

    roles = (case.get("assembly") or {}).get("roles") or {}
    role_hashes = {k: _file_sha256(asset_path(case, v)) for k, v in roles.items() if v}

    meta = {
        "schema_version": "prepared-2026-11-24",
        "fingerprint": fp,
        "stirring_enabled": stir,
        "mechanics_mode": mechanics_mode(case) if stir else "none",
        "dp_m": dp,
        "melt_height_m": melt_h,
        "reload_dir": "reload",
        "source_binary": str(bin_path),
        "role_content_sha256": role_hashes,
        "reload_audit": audit,
        "note": "Reload copied from solver prepare; fingerprint includes STL content hashes + dp",
    }
    (prepared / "prepared.json").write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
    (output_dir / "latest").write_text(fp + "\n", encoding="utf-8")
    emit("end", stage="prepare", status="completed", fingerprint=fp, prepared=str(prepared), reload_audit=audit)
    return 0
