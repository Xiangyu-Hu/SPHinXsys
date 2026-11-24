from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
from pathlib import Path
from typing import Any

from .case_io import mechanics_mode, stirring_enabled
from .events import emit, log
from .paths import EXTRA_3D, resolve_binary
from .validate import validate_case


def _fingerprint(case: dict[str, Any]) -> str:
    payload = {
        "stirring": stirring_enabled(case),
        "mechanics": mechanics_mode(case),
        "dp": (case.get("numerics") or {}).get("dp_m"),
        "melt_height": (case.get("geometry") or {}).get("melt_height_m"),
        "roles": (case.get("assembly") or {}).get("roles"),
    }
    raw = json.dumps(payload, sort_keys=True).encode()
    return hashlib.sha256(raw).hexdigest()[:16]


def _run(cmd: list[str], cwd: Path) -> int:
    log(f"$ {' '.join(cmd)}  (cwd={cwd})")
    proc = subprocess.run(cmd, cwd=str(cwd))
    return proc.returncode


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
        rc = _run([str(bin_path), "--bodies=all"], work)
        if rc != 0:
            emit("error", stage="prepare", code=rc, message="stir_relax failed")
            return rc
        reload_src = work / "reload"

    reload_dst = prepared / "reload"
    if reload_dst.exists():
        shutil.rmtree(reload_dst)
    shutil.copytree(reload_src, reload_dst)

    meta = {
        "schema_version": "prepared-2026-11-21",
        "fingerprint": fp,
        "stirring_enabled": stir,
        "mechanics_mode": mechanics_mode(case) if stir else "none",
        "dp_m": dp,
        "melt_height_m": melt_h,
        "reload_dir": "reload",
        "source_binary": str(bin_path),
        "note": "Reload copied from solver prepare run; display particles must match this fingerprint",
    }
    (prepared / "prepared.json").write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
    # convenience pointer
    (output_dir / "latest").write_text(fp + "\n", encoding="utf-8")
    emit("end", stage="prepare", status="completed", fingerprint=fp, prepared=str(prepared))
    return 0
