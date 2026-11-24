from __future__ import annotations

import json
from pathlib import Path
from typing import Any


def inspect_run(run_dir: Path) -> dict[str, Any]:
    run_dir = run_dir.resolve()
    out: dict[str, Any] = {
        "run_dir": str(run_dir),
        "exists": run_dir.is_dir(),
        "files": {},
        "status": "unknown",
        "warnings": [],
    }
    if not run_dir.is_dir():
        out["status"] = "missing"
        return out

    for name in (
        "result-index.json",
        "summary.json",
        "effective-case.json",
        "events.ndjson",
        "solver.log",
        "probes.csv",
        "results.pvd",
    ):
        p = run_dir / name
        out["files"][name] = {"present": p.is_file(), "size": p.stat().st_size if p.is_file() else 0}

    idx_path = run_dir / "result-index.json"
    if idx_path.is_file():
        idx = json.loads(idx_path.read_text(encoding="utf-8"))
        out["result_index"] = idx
        out["status"] = idx.get("status", "unknown")
        if idx.get("full_revolution_pass") == 1:
            out["warnings"].append("index claims full_revolution_pass=1; verify this is not a short-window mislabel")
        if "HEATED_ELASTIC_INTEGRATION_PASS" in idx and idx.get("physics", "").find("NOT_EM") < 0:
            # soft note only if physics label missing
            pass

    summary_path = run_dir / "summary.json"
    if summary_path.is_file():
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        out["summary_keys"] = sorted(summary.keys())
        out["passed"] = summary.get("passed")
        out["physical_time_s"] = summary.get("physical_time_s")
        out["mechanics_mode"] = summary.get("mechanics_mode")
        # anti-overclaim helpers
        if summary.get("full_revolution_pass") == 0:
            out["warnings"].append("full_revolution_pass=0 (short window only)")
        if summary.get("HEATED_ELASTIC_INTEGRATION_PASS") == 1:
            out["warnings"].append("F1 heated integration PASS is prescribed heat, NOT real EM")

    effective_path = run_dir / "effective-case.json"
    if effective_path.is_file():
        eff = json.loads(effective_path.read_text(encoding="utf-8"))
        out["effective"] = {
            "stirring_enabled": eff.get("stirring_enabled"),
            "mechanics_mode": eff.get("mechanics_mode"),
            "heating_mode": eff.get("heating_mode"),
            "claims_forbidden": eff.get("claims_forbidden"),
        }

    frames = sorted(run_dir.glob("frames/*.vtp"))
    if not frames:
        frames = sorted(run_dir.glob("**/*.vtp"))
    out["vtp_count"] = len(frames)
    out["vtp_samples"] = [str(p.relative_to(run_dir)) for p in frames[:12]]
    # group by body prefix for GUI replay lists
    bodies: dict[str, list[str]] = {}
    for p in frames:
        name = p.name
        body = name.split("_ite_")[0] if "_ite_" in name else p.stem
        bodies.setdefault(body, []).append(str(p.relative_to(run_dir)))
    out["bodies"] = {k: {"count": len(v), "first": v[0], "last": v[-1]} for k, v in bodies.items()}
    csvs = sorted([str(p.relative_to(run_dir)) for p in run_dir.glob("*.csv")])
    out["csv_files"] = csvs
    return out


def import_history(src_summary: Path, dest_run: Path, label: str) -> dict[str, Any]:
    """Attach a historical handoff summary as a clearly marked replay-only run."""
    dest_run.mkdir(parents=True, exist_ok=True)
    summary = json.loads(src_summary.read_text(encoding="utf-8"))
    (dest_run / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    index = {
        "status": "history_replay",
        "source": str(src_summary),
        "label": label,
        "not_a_fresh_local_run": True,
        "summary_passed": summary.get("passed"),
        "physical_time_s": summary.get("physical_time_s"),
        "mechanics_mode": summary.get("mechanics_mode"),
        "note": "Historical evidence for GUI/CLI replay demos; do not claim as this-machine fresh PASS unless re-run",
    }
    for k in (
        "ELASTIC_FSI_SHORT_PASS",
        "HEATED_ELASTIC_INTEGRATION_PASS",
        "DRIVEN_ELASTIC_ROTATION_SHORT_PASS",
        "full_revolution_pass",
    ):
        if k in summary:
            index[k] = summary[k]
    (dest_run / "result-index.json").write_text(json.dumps(index, indent=2) + "\n", encoding="utf-8")
    (dest_run / "effective-case.json").write_text(
        json.dumps(
            {
                "history_import": True,
                "source_summary": str(src_summary),
                "claims_forbidden": [
                    "Do not label history replay as fresh GPU run on this host without re-run",
                    "Do not label E2 short window as full revolution",
                    "Do not label F1 prescribed heat as solved EM",
                ],
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    return index
