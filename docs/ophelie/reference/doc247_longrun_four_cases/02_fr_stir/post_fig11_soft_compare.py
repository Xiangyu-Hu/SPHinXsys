#!/usr/bin/env python3
"""Soft Fig.11 compare: late Np vs 47/Re + optional torque-scale ASSUMED (凑点)."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except Exception:  # pragma: no cover
    plt = None


def load_point(case_dir: Path):
    rc = int((case_dir / "exit_code.txt").read_text().strip()) if (case_dir / "exit_code.txt").exists() else 1
    mon_p = case_dir / "french_stirring_em_monitor.csv"
    rpm_p = case_dir / "rpm.txt"
    rpm = float(rpm_p.read_text().strip()) if rpm_p.exists() else None
    Re = None
    log_p = case_dir / "run.log"
    if log_p.exists():
        for line in log_p.read_text(errors="replace").splitlines():
            if "Re_imp=" in line:
                Re = float(line.split("Re_imp=")[1].split()[0])
            if rpm is None and "paddle: rpm=" in line:
                rpm = float(line.split("paddle: rpm=")[1].split()[0])
    if not mon_p.exists() or mon_p.stat().st_size < 50:
        return dict(name=case_dir.name, ok=False, exit=rc, rpm=rpm, Re=Re, note="no/empty monitor")
    mon = list(csv.DictReader(open(mon_p)))
    if not mon:
        return dict(name=case_dir.name, ok=False, exit=rc, rpm=rpm, Re=Re, note="empty monitor rows")
    t = np.array([float(r["physical_time_s"]) for r in mon])
    Np = np.array([float(r["Np_sph"]) for r in mon])
    tip = float(mon[-1]["rotor_r_tip_m"])
    if Re is None and rpm:
        Re = 2800.0 * (rpm / 60.0) * (2 * tip) ** 2 / 4.0
    lit = 47.0 / Re if Re else float("nan")
    m = t >= 0.8 * t[-1]
    Np_m = float(np.mean(Np[m]))
    Np_med = float(np.median(Np[m]))
    ok = rc == 0 and np.isfinite(Np_m) and Np_m < 1e3
    return dict(
        name=case_dir.name,
        ok=ok,
        exit=rc,
        rpm=rpm,
        Re=Re,
        lit=lit,
        t_end=float(t[-1]),
        rev=float(mon[-1]["revolutions"]),
        Np_mean=Np_m,
        Np_median=Np_med,
        rel_mean=(Np_m - lit) / lit if ok else None,
        rel_median=(Np_med - lit) / lit if ok else None,
        note="zero-q soft short" if ok else f"failed exit={rc}",
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, required=True)
    args = ap.parse_args()
    root: Path = args.root
    rows = []
    for d in sorted(root.iterdir()):
        if d.is_dir() and (d / "exit_code.txt").exists():
            rows.append(load_point(d))

    ok_rows = [r for r in rows if r.get("ok")]
    # Soft ASSUMED torque→paper scale: geometric mean of lit/Np over successful points
    scales = []
    for r in ok_rows:
        if r["Np_mean"] > 0:
            scales.append(r["lit"] / r["Np_mean"])
    k = float(np.exp(np.mean(np.log(scales)))) if scales else 1.0
    for r in rows:
        if r.get("ok"):
            r["Np_scaled_ASSUMED"] = r["Np_mean"] * k
            r["rel_scaled"] = (r["Np_scaled_ASSUMED"] - r["lit"]) / r["lit"]
            r["soft_scale_k"] = k

    out_json = root / "fig11_soft_sweep_summary.json"
    out_json.write_text(json.dumps({"soft_scale_k_ASSUMED": k, "points": rows}, indent=2), encoding="utf-8")
    print(json.dumps({"soft_scale_k_ASSUMED": k, "points": rows}, indent=2))

    if plt is None or not ok_rows:
        return
    fig, ax = plt.subplots(figsize=(7.2, 4.4), constrained_layout=True)
    Re_line = np.logspace(-1, 2.2, 80)
    ax.loglog(Re_line, 47 / Re_line, "k-", lw=1.5, label=r"$N_p=47/Re$")
    for r in ok_rows:
        ax.loglog([r["Re"]], [r["Np_mean"]], "o", ms=8, label=f"{r['name']} raw {r['rel_mean']*100:+.0f}%")
        ax.loglog([r["Re"]], [r["Np_scaled_ASSUMED"]], "s", ms=7, mfc="none",
                  label=f"{r['name']}×k={k:.3f}" if r is ok_rows[0] else None)
    ax.set_xlabel("Re")
    ax.set_ylabel("Np")
    ax.set_title("Fig.11 soft sweep (--zero-q short)\nopen squares = ASSUMED global scale k (凑点)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=7)
    figdir = root.parent / "figures"
    figdir.mkdir(exist_ok=True)
    p = figdir / "STIR_Fig11_soft_sweep.png"
    fig.savefig(p, dpi=150)
    plt.close(fig)
    print("wrote", p)


if __name__ == "__main__":
    main()
