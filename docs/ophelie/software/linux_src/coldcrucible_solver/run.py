from __future__ import annotations

import json
import shutil
import subprocess
import time
from pathlib import Path
from typing import Any

from .case_io import asset_path, mechanics_mode, stirring_enabled
from .events import emit, log
from .paths import resolve_binary
from .validate import validate_case


def _resolve_prepared(prepared_dir: Path) -> Path:
    prep = prepared_dir
    if (prep / "prepared.json").is_file():
        return prep
    latest = prep / "latest"
    if latest.is_file():
        fp = latest.read_text(encoding="utf-8").strip()
        cand = prep / fp
        if (cand / "prepared.json").is_file():
            return cand
    kids = [p for p in prep.iterdir() if p.is_dir() and (p / "prepared.json").is_file()] if prep.is_dir() else []
    if len(kids) == 1:
        return kids[0]
    raise FileNotFoundError(f"prepared.json not found under {prepared_dir}")


def _build_nc_cmd(
    case: dict[str, Any],
    binary: Path,
    reload_dir: Path,
    summary_path: Path,
) -> tuple[list[str], str, dict[str, Any]]:
    numerics = case.get("numerics") or {}
    geom = case.get("geometry") or {}
    run_cfg = case.get("run") or {}
    heating = case.get("heating") or {}
    dp = float(numerics.get("dp_m", 0.010))
    melt_h = float(geom.get("melt_height_m", 0.180))
    end_time = float(run_cfg.get("end_time_s", 1.0))
    state_iv = float(run_cfg.get("state_record_interval_s", 0.25))
    profile = run_cfg.get("profile", "isothermal_hydro")
    heat_mode = heating.get("mode", "none")

    cmd = [
        str(binary),
        "--client-nc-geometry",
        "--client-nc-material",
        f"--dp={dp}",
        f"--melt-height={melt_h}",
        f"--reload-dir={reload_dir}",
        f"--summary-json={summary_path}",
        "--doc227-actual-timestep",
        "--boundary-face-t-mode=local_dt_dz",
    ]
    meta: dict[str, Any] = {
        "profile": profile,
        "heating_mode": heat_mode,
        "zero_q": False,
        "top_surface_model": geom.get("top_surface_model"),
        "top_surface_identity": geom.get("top_surface_identity"),
    }

    if profile == "isothermal_hydro":
        if heat_mode not in ("none",):
            raise ValueError("isothermal_hydro requires heating.mode=none")
        cmd += [
            "--zero-q",
            "--zero-buoyancy",
            f"--end-time={end_time}",
            "--hydro-settle-time=0.15",
            f"--state-record-interval={state_iv}",
        ]
        meta["zero_q"] = True
        phys = "nc_isothermal_hydro_diagnostic"
    elif profile == "dT20_boussinesq":
        if heat_mode not in ("none",):
            raise ValueError("dT20_boussinesq requires heating.mode=none")
        cmd += [
            "--zero-q",
            "--doc223-paper-flat-slip-top",
            "--prescribe-t-radial=20",
            "--freeze-prescribed-t",
            f"--end-time={end_time}",
            "--hydro-settle-time=1.0",
            f"--state-record-interval={state_iv}",
        ]
        meta["zero_q"] = True
        meta["top_surface_model"] = "paper_flat_slip_top"
        phys = "nc_dT20_boussinesq_diagnostic"
    elif profile == "fixed_absorbed_power_em":
        if heat_mode != "fixed_absorbed_power":
            raise ValueError("fixed_absorbed_power_em requires heating.mode=fixed_absorbed_power")
        target = float(heating["target_absorbed_power_W"])
        em_iv = float(heating.get("em_update_interval_s", 0.1))
        aind = str(heating.get("aind", "off"))
        cmd += [
            "--doc223-paper-flat-slip-top",
            "--em-control=fixed-power",
            f"--target-absorbed-power={target}",
            "--thermo-em-coupling=periodic",
            f"--em-update-interval={em_iv}",
            f"--aind={aind}",
            "--enable-thermal-diffusion",
            f"--end-time={end_time}",
            "--hydro-settle-time=1.0",
            f"--state-record-interval={max(state_iv, end_time / 4.0)}",
            "--meridional-average-start=0.5",
        ]
        # Explicitly refuse zero-q
        if "--zero-q" in cmd:
            raise RuntimeError("internal error: zero-q leaked into EM path")
        meta.update(
            {
                "zero_q": False,
                "target_absorbed_power_W": target,
                "power_status": heating.get("power_status", "ASSUMED"),
                "generator_input_power_W": heating.get("generator_input_power_W"),
                "current_status": heating.get("current_status", "unknown"),
                "top_surface_model": "paper_flat_slip_top",
                "top_surface_identity": geom.get("top_surface_identity") or "CLIENT_NC_FLAT_SLIP_ENGINEERING",
            }
        )
        phys = "nc_fixed_absorbed_power_em_ASSUMED"
    else:
        raise ValueError(f"unsupported NC profile={profile}")

    return cmd, phys, meta


def run(case: dict[str, Any], prepared_dir: Path, run_dir: Path) -> int:
    report = validate_case(case)
    if not report["ok"]:
        emit("error", stage="validate", issues=report["issues"])
        return 2

    prepared = _resolve_prepared(prepared_dir).resolve()
    reload_dir = (prepared / "reload").resolve()
    if not reload_dir.is_dir():
        emit("error", stage="run", message=f"reload missing: {reload_dir}")
        return 2
    if not (reload_dir / "Reload.xml").is_file():
        emit("error", stage="run", message=f"Reload.xml missing under {reload_dir}")
        return 2

    # prepared dp consistency
    try:
        prep_meta = json.loads((prepared / "prepared.json").read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        prep_meta = {}
    case_dp = float((case.get("numerics") or {}).get("dp_m", 0.0) or 0.0)
    prep_dp = prep_meta.get("dp_m")
    if prep_dp is not None and case_dp > 0 and abs(float(prep_dp) - case_dp) > 1e-9:
        emit(
            "error",
            stage="run",
            message=f"dp mismatch: case dp={case_dp} vs prepared dp={prep_dp}; refuse silent reuse",
        )
        return 2

    run_dir = run_dir.resolve()
    run_dir.mkdir(parents=True, exist_ok=True)
    events_path = run_dir / "events.ndjson"
    log_path = run_dir / "solver.log"
    summary_path = run_dir / "summary.json"
    effective_path = run_dir / "effective-case.json"

    stir = stirring_enabled(case)
    numerics = case.get("numerics") or {}
    geom = case.get("geometry") or {}
    run_cfg = case.get("run") or {}
    heating = case.get("heating") or {}
    stirring = case.get("stirring") or {}
    dp = float(numerics.get("dp_m", 0.010 if not stir else 0.005))
    melt_h = float(geom.get("melt_height_m", 0.180))
    end_time = float(run_cfg.get("end_time_s", 0.10 if not stir else 0.05))
    omega_ramp = float(stirring.get("omega_ramp_s", run_cfg.get("omega_ramp_s", 0.05)))
    clamp_top = float(stirring.get("clamp_top_band_m", run_cfg.get("clamp_top_band_m", 0.05)))

    if not stir:
        binary = resolve_binary("nc_run")
        try:
            cmd, phys_label, nc_meta = _build_nc_cmd(case, binary, reload_dir, summary_path)
        except ValueError as e:
            emit("error", stage="run", message=str(e))
            return 2
        work = binary.parent
        effective_extra = nc_meta
    else:
        heat_mode = heating.get("mode", "none")
        if heat_mode in ("prescribed_volume_heat", "f1"):
            binary = resolve_binary("stir_f1")
            if not binary.is_file():
                emit("error", stage="run", message=f"F1 binary missing: {binary}")
                return 2
            phys_label = "stir_f1_prescribed_heat_NOT_EM"
        else:
            binary = resolve_binary("stir_f0")
            phys_label = "stir_f0_isothermal_elastic_fsi"
        roles = (case.get("assembly") or {}).get("roles") or {}
        rotor_rel = roles.get("paddle_solver") or roles.get("paddle")
        rotor = asset_path(case, rotor_rel)
        # Prefer project-selected paddle; only fall back to binary/data with explicit warning.
        data_rotor = binary.parent / "data" / "client_stir_paddle_solver_m.stl"
        if rotor.is_file():
            rotor_arg = str(rotor)
            if data_rotor.is_file():
                emit(
                    "warning",
                    stage="run",
                    message="using case paddle_solver STL; binary/data copy ignored for geometry identity",
                )
        elif data_rotor.is_file():
            rotor_arg = str(data_rotor)
            emit(
                "warning",
                stage="run",
                message="case paddle STL missing; falling back to binary/data STL (validate should have caught this)",
            )
        else:
            emit("error", stage="run", message=f"paddle STL missing: {rotor}")
            return 2
        cmd = [
            str(binary),
            f"--reload-dir={reload_dir}",
            f"--rotor-stl={rotor_arg}",
            f"--end-time={end_time}",
            f"--omega-ramp={omega_ramp}",
            f"--clamp-top-band={clamp_top}",
            "--c0=8",
            "--max-wall-hours=2",
            f"--out-dir={run_dir / 'solver_out'}",
            f"--summary-json={summary_path}",
            f"--state-record-interval={float(run_cfg.get('state_record_interval_s', max(end_time / 4.0, 0.01)))}",
        ]
        if "f1" in phys_label:
            heat_w = float(heating.get("prescribed_heat_power_W", 1000.0))
            cmd.append(f"--heat-power={heat_w}")
            emit(
                "warning",
                stage="run",
                message="prescribed volume heat path; NOT real EM / NOT verified 30 kW absorbed",
            )
        work = binary.parent
        effective_extra = {
            "omega_ramp_s": omega_ramp,
            "clamp_top_band_m": clamp_top,
            "rotor_stl": rotor_arg,
            "heating_mode": heat_mode,
        }

    effective = {
        "case_id": case.get("case_id") or case.get("project_name"),
        "stirring_enabled": stir,
        "mechanics_mode": mechanics_mode(case) if stir else "none",
        "dp_m": dp,
        "melt_height_m": melt_h,
        "end_time_s": end_time,
        "heating_mode": heating.get("mode", "none"),
        "material_id": (case.get("material") or {}).get("material_id"),
        "claims_forbidden": [
            "Do not label E2 short window as full revolution",
            "Do not label prescribed heat (F1) as solved EM",
            "30 kW is generator input only",
            "ASSUMED absorbed power is not experimental verification",
        ],
        "prepared_fingerprint": prep_meta.get("fingerprint"),
        "physics": phys_label,
        **effective_extra,
    }
    effective_path.write_text(json.dumps(effective, indent=2) + "\n", encoding="utf-8")

    cancel_path = run_dir / "cancel.request"
    if cancel_path.exists():
        cancel_path.unlink()
    (run_dir / "HOW_TO_CANCEL.txt").write_text(
        "To stop this run from another shell:\n"
        f"  touch {cancel_path}\n"
        "Solver is started in a new session so SSH disconnect should not kill it;\n"
        "still prefer cancel.request over kill -9.\n",
        encoding="utf-8",
    )

    emit("start", stage="run", physics=phys_label, device_selector="cuda:gpu", cmd=cmd, cancel_file=str(cancel_path))
    t0 = time.time()
    cancelled = False
    with events_path.open("a", encoding="utf-8") as evf, log_path.open("w", encoding="utf-8") as lf:
        def _tee_event(obj: dict[str, Any]) -> None:
            line = json.dumps(obj, ensure_ascii=False)
            print(line, flush=True)
            evf.write(line + "\n")
            evf.flush()

        _tee_event({"event": "progress", "stage": "flow", "physical_time_s": 0.0, "message": "solver launching"})
        log(f"$ {' '.join(cmd)}")
        proc = subprocess.Popen(
            cmd,
            cwd=str(work),
            stdout=lf,
            stderr=subprocess.STDOUT,
            start_new_session=True,
        )
        (run_dir / "solver.pid").write_text(str(proc.pid) + "\n", encoding="utf-8")
        while proc.poll() is None:
            if cancel_path.is_file():
                cancelled = True
                _tee_event({"event": "warning", "stage": "run", "message": "cancel.request seen; terminating solver"})
                proc.terminate()
                try:
                    proc.wait(timeout=15)
                except subprocess.TimeoutExpired:
                    proc.kill()
                    proc.wait(timeout=5)
                break
            time.sleep(0.5)
        rc = proc.returncode if proc.returncode is not None else (130 if cancelled else 1)
        wall = time.time() - t0

    if cancelled:
        status = "cancelled"
    elif rc == 0:
        status = "completed"
    else:
        status = "failed"

    frames_dir = run_dir / "frames"
    collected = {"vtp": 0, "csv": 0, "pvd": 0, "sources": []}
    candidate_dirs = [work / "output", run_dir / "solver_out"]
    for src_dir in candidate_dirs:
        if not src_dir.is_dir():
            continue
        collected["sources"].append(str(src_dir))
        frames_dir.mkdir(parents=True, exist_ok=True)
        for pattern, key in (("*.vtp", "vtp"), ("*.csv", "csv"), ("*.pvd", "pvd"), ("*.dat", "csv")):
            for f in src_dir.glob(pattern):
                dest = frames_dir / f.name if key == "vtp" else run_dir / f.name
                try:
                    shutil.copy2(f, dest)
                    collected[key if key in collected else "csv"] = collected.get(key, 0) + 1
                except OSError as e:
                    log(f"copy failed {f}: {e}")

    index: dict[str, Any] = {
        "status": status,
        "exit_code": rc,
        "physics": phys_label,
        "wall_time_s": wall,
        "summary": "summary.json" if summary_path.is_file() else None,
        "effective_case": "effective-case.json",
        "solver_log": "solver.log",
        "claims_scope": effective["claims_forbidden"],
        "cancel_file": "cancel.request",
        "frames_dir": "frames" if frames_dir.is_dir() else None,
        "collected": collected,
        "results_root_note": "outputs are under this run_dir; install/prefix stays read-only",
        "capability_tags": {},
    }
    if summary_path.is_file():
        try:
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
            index["summary_passed"] = summary.get("passed")
            index["physical_time_s"] = summary.get("end_time") or summary.get("physical_time_s")
            tags = index["capability_tags"]
            if summary.get("em_status") == "RUN" or summary.get("P_joule_W"):
                tags["EM_EXECUTED"] = 1 if summary.get("em_status") == "RUN" else 0
                tags["P_joule_W"] = summary.get("P_joule_W")
                tags["P_heat_fluid_W"] = summary.get("P_heat_fluid_W")
                tags["target_power_w"] = summary.get("target_power_w")
            tags["hydro_model"] = summary.get("hydro_model")
            tags["numerical_quality_ok"] = summary.get("numerical_quality_ok")
            tags["reached_end_time"] = summary.get("reached_end_time")
            for k in (
                "ELASTIC_FSI_SHORT_PASS",
                "HEATED_ELASTIC_INTEGRATION_PASS",
                "full_revolution_pass",
                "NUMERICAL_QUALITY_PASS",
            ):
                if k in summary:
                    tags[k] = summary[k]
            if effective.get("power_status") == "ASSUMED":
                tags["power_status"] = "ASSUMED"
                tags["target_absorbed_power_W"] = effective.get("target_absorbed_power_W")
        except json.JSONDecodeError:
            pass

    if frames_dir.is_dir():
        vtps = sorted(frames_dir.glob("*.vtp"))
        index["frames"] = [f"frames/{p.name}" for p in vtps]
        index["frame_count"] = len(vtps)
        bodies: dict[str, list[str]] = {}
        for p in vtps:
            body = p.name.split("_ite_")[0] if "_ite_" in p.name else p.stem
            bodies.setdefault(body, []).append(f"frames/{p.name}")
        index["bodies"] = {k: {"count": len(v), "files": v} for k, v in bodies.items()}
    index["csv_files"] = sorted(p.name for p in run_dir.glob("*.csv"))

    (run_dir / "result-index.json").write_text(json.dumps(index, indent=2) + "\n", encoding="utf-8")

    emit("end", stage="run", status=index["status"], exit_code=rc, wall_time_s=wall, run_dir=str(run_dir), collected=collected)
    return rc
