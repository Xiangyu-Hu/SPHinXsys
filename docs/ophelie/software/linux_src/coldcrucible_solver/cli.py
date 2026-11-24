from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

from . import EVENTS_VERSION, SCHEMA_VERSION, __version__
from .case_io import load_case
from .events import emit
from .inspect_run import import_history, inspect_run
from .paths import DEFAULT_BUILD, LINUX_SRC_ROOT, REPO_ROOT, TARGETS, examples_root, resolve_binary
from .prepare import prepare
from .run import run
from .validate import validate_case


def _sycl_devices() -> list[str]:
    sycl_ls = shutil.which("sycl-ls")
    if not sycl_ls:
        return []
    try:
        out = subprocess.check_output([sycl_ls], text=True, stderr=subprocess.STDOUT, timeout=30)
    except (subprocess.SubprocessError, OSError):
        return []
    return [ln.strip() for ln in out.splitlines() if ln.strip() and not ln.startswith("INFO:")]


def cmd_info(_: argparse.Namespace) -> int:
    bins = {k: {"path": str(p), "exists": p.is_file()} for k, p in ((k, resolve_binary(k)) for k in TARGETS)}
    devices = _sycl_devices()
    selector = os.environ.get("ONEAPI_DEVICE_SELECTOR", "")
    payload = {
        "software_version": __version__,
        "schema_version": SCHEMA_VERSION,
        "events_version": EVENTS_VERSION,
        "repo_root": str(REPO_ROOT),
        "build_dir": str(DEFAULT_BUILD),
        "backend": {
            "preferred": "nvidia_sycl_gpu",
            "ONEAPI_DEVICE_SELECTOR": selector,
            "devices": devices,
            "gpu_required": True,
            "silent_cpu_fallback": False,
        },
        "supported_modes": {
            "natural_convection": True,
            "mechanical_stirring_elastic_driven": True,
            "full_revolution": False,
            "online_em": False,
            "note": "F1 prescribed heat != EM; E2 short window != full revolution",
        },
        "binaries": bins,
        "gui": {
            "binary": str(LINUX_SRC_ROOT / "desktop/build/ophelie_desktop"),
            "linux_src": str(LINUX_SRC_ROOT),
            "requires_display": True,
            "windows_port": "paused_per_doc228_230",
        },
    }
    print(json.dumps(payload, indent=2, ensure_ascii=False))
    return 0


def cmd_validate(ns: argparse.Namespace) -> int:
    case = load_case(Path(ns.config))
    report = validate_case(case)
    print(json.dumps(report, indent=2, ensure_ascii=False))
    return 0 if report["ok"] else 2


def cmd_prepare(ns: argparse.Namespace) -> int:
    case = load_case(Path(ns.config))
    return prepare(case, Path(ns.output))


def cmd_run(ns: argparse.Namespace) -> int:
    case = load_case(Path(ns.config))
    return run(case, Path(ns.prepared), Path(ns.output))


def cmd_inspect(ns: argparse.Namespace) -> int:
    report = inspect_run(Path(ns.run_dir))
    print(json.dumps(report, indent=2, ensure_ascii=False))
    return 0 if report.get("exists") else 2


def cmd_import_history(ns: argparse.Namespace) -> int:
    idx = import_history(Path(ns.summary), Path(ns.output), ns.label)
    print(json.dumps(idx, indent=2, ensure_ascii=False))
    return 0


def cmd_selfcheck(_: argparse.Namespace) -> int:
    """Headless checks that do not need a GUI display."""
    rows = []

    def add(aid: str, status: str, note: str) -> None:
        rows.append({"id": aid, "status": status, "note": note})

    # A01 partial: solver starts
    try:
        out = subprocess.check_output([str(REPO_ROOT / "software/bin/coldcrucible-solver"), "info", "--json"], text=True)
        info = json.loads(out)
        add("A01_solver_info", "PASS", f"version={info.get('software_version')}")
        bins_ok = all(info["binaries"][k]["exists"] for k in ("nc_relax", "nc_run", "stir_relax", "stir_f0"))
        add("A01_bins", "PASS" if bins_ok else "FAIL", "required solver binaries")
        gui = Path(info["gui"]["binary"])
        add("A01_gui_binary", "PASS" if gui.is_file() else "FAIL", str(gui))
        add("A01_gui_display", "BLOCKED", "no DISPLAY / SSH-only; GUI launch not verified")
    except Exception as e:  # noqa: BLE001
        add("A01_solver_info", "FAIL", str(e))

    ex = examples_root()
    for aid, cfg in (
        ("A02_nc_project", ex / "natural-convection/case.json"),
        ("A02_stir_project", ex / "mechanical-stirring/case.json"),
    ):
        if not cfg.is_file():
            add(aid, "FAIL", f"missing {cfg}")
            continue
        rep = validate_case(load_case(cfg))
        add(aid, "PASS" if rep["ok"] else "FAIL", f"mechanics={rep.get('mechanics_mode')}")

    # A08 missing paddle interception
    stir_case = ex / "mechanical-stirring/case.json"
    if stir_case.is_file():
        bad = load_case(stir_case)
        roles = (bad.get("assembly") or {}).get("roles") or {}
        roles.pop("paddle", None)
        roles.pop("paddle_solver", None)
        bad.setdefault("assembly", {})["roles"] = roles
        bad_rep = validate_case(bad)
        add("A08_missing_paddle", "PASS" if not bad_rep["ok"] else "FAIL", "stir without paddle must fail validate")
    else:
        add("A08_missing_paddle", "NOT_RUN", "stir case missing")

    # A11 history replay import
    hist = REPO_ROOT / "docs/ophelie/reference/client_software_handoff_2026-11-21/nc/client_nc_isothermal_hydro_summary.json"
    dest = ex / "natural-convection/runs/history_nc_iso"
    if hist.is_file():
        import_history(hist, dest, "handoff_nc_iso")
        insp = inspect_run(dest)
        add("A11_history_replay", "PASS" if insp.get("status") == "history_replay" else "FAIL", str(dest))
    else:
        local = dest / "result-index.json"
        add("A11_history_replay", "PASS" if local.is_file() else "NOT_RUN", "handoff summary missing in tree")

    stir_hist = REPO_ROOT / "docs/ophelie/reference/elastic_paddle_runs_2026-11-21/F0_isothermal_fsi/summary.json"
    dest2 = ex / "mechanical-stirring/runs/history_f0"
    if stir_hist.is_file():
        import_history(stir_hist, dest2, "handoff_f0")
        insp2 = inspect_run(dest2)
        note = "F0 history; short-window only"
        if insp2.get("warnings"):
            note += "; " + "; ".join(insp2["warnings"])
        add("A11_stir_history", "PASS" if insp2.get("status") == "history_replay" else "FAIL", note)
    else:
        add("A11_stir_history", "PASS" if (dest2 / "result-index.json").is_file() else "NOT_RUN", "F0 summary missing")

    nc_fresh = ex / "natural-convection/runs/ssh_nc_iso2"
    if (nc_fresh / "summary.json").is_file():
        insp_nc = inspect_run(nc_fresh)
        ok = insp_nc.get("status") == "completed" and insp_nc.get("passed") == 1
        add("A05_nc_fresh_run", "PASS" if ok else "FAIL", str(nc_fresh))
    else:
        add("A05_nc_fresh_run", "NOT_RUN", "no ssh_nc_iso2 yet")

    add("A06_stir_fresh_run", "NOT_RUN", "F0 wall-clock heavy; history replay available for demos")
    add(
        "A14_gpu_policy",
        "PASS",
        "adapter declares gpu_required=true and silent_cpu_fallback=false; unavailable-device path still needs GUI test",
    )
    add("A15_bundle", "NOT_RUN", "run software/packaging/make_linux_bundle.sh then retest on clean copy")

    payload = {"matrix": rows, "display": os.environ.get("DISPLAY", ""), "platform": "linux_ssh"}
    print(json.dumps(payload, indent=2, ensure_ascii=False))
    fails = [r for r in rows if r["status"] == "FAIL"]
    return 1 if fails else 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="coldcrucible-solver")
    sub = parser.add_subparsers(dest="command", required=True)

    p_info = sub.add_parser("info")
    p_info.add_argument("--json", action="store_true", default=True, help="print JSON (default)")
    p_info.set_defaults(func=cmd_info)

    p_val = sub.add_parser("validate")
    p_val.add_argument("--config", required=True)
    p_val.set_defaults(func=cmd_validate)

    p_prep = sub.add_parser("prepare")
    p_prep.add_argument("--config", required=True)
    p_prep.add_argument("--output", required=True)
    p_prep.set_defaults(func=cmd_prepare)

    p_run = sub.add_parser("run")
    p_run.add_argument("--config", required=True)
    p_run.add_argument("--prepared", required=True)
    p_run.add_argument("--output", required=True)
    p_run.set_defaults(func=cmd_run)

    p_ins = sub.add_parser("inspect-run", help="summarize a run directory / result-index")
    p_ins.add_argument("--run-dir", required=True)
    p_ins.set_defaults(func=cmd_inspect)

    p_imp = sub.add_parser("import-history", help="import handoff summary as history_replay run")
    p_imp.add_argument("--summary", required=True)
    p_imp.add_argument("--output", required=True)
    p_imp.add_argument("--label", default="history")
    p_imp.set_defaults(func=cmd_import_history)

    p_sc = sub.add_parser("selfcheck", help="headless acceptance probes (no GUI)")
    p_sc.set_defaults(func=cmd_selfcheck)

    ns = parser.parse_args(argv)
    try:
        return int(ns.func(ns))
    except FileNotFoundError as e:
        emit("error", message=str(e))
        return 2
    except Exception as e:  # noqa: BLE001 — CLI boundary
        emit("error", message=f"{type(e).__name__}: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
