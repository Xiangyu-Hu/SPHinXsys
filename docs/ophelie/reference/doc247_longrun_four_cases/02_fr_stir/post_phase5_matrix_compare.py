#!/usr/bin/env python3
"""Phase5 matrix post: MATRIX.md + optional T/Q compare PNGs at selected times.

Soft engineering compare only (ASSUMED / PLACEHOLDER Fig.10). Not literature PASS.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path

import numpy as np

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except Exception:  # pragma: no cover
    plt = None


CASE_ORDER = [
    "SW_rpm05",
    "SW_rpm10",
    "SW_rpm20",
    "EM_P56",
    "EM_P60",
    "EM_P70",
    "ST_d3",
    "ST_d5",
    "ST_d8",
]


def _load_rz_csv(path: Path, value_keys=("T", "T_K", "Q_norm", "Q_Wpm3", "Q", "value")):
    if not path.exists() or path.stat().st_size < 20:
        return None
    with path.open() as f:
        rdr = csv.DictReader(f)
        if not rdr.fieldnames:
            return None
        hdr = {h.lower(): h for h in rdr.fieldnames}
        r_key = next((hdr[k] for k in ("r_m", "r", "r_over_r") if k in hdr), None)
        z_key = next((hdr[k] for k in ("z_m", "z", "z_over_h") if k in hdr), None)
        v_key = next((hdr[k] for k in value_keys if k.lower() in hdr), None)
        if v_key is None:
            # try exact case from keys
            for k in value_keys:
                if k in rdr.fieldnames:
                    v_key = k
                    break
        if not (r_key and z_key and v_key):
            return None
        rows = list(rdr)
    if not rows:
        return None
    r = np.array([float(row[r_key]) for row in rows], dtype=float)
    z = np.array([float(row[z_key]) for row in rows], dtype=float)
    v = np.array([float(row[v_key]) for row in rows], dtype=float)
    return dict(r=r, z=z, v=v, path=str(path))


def _grid_from_scatter(r, z, v):
    ru = np.unique(np.round(r, 12))
    zu = np.unique(np.round(z, 12))
    if ru.size * zu.size != r.size:
        # irregular / folded — return scatter only
        return None, None, None
    R, Z = np.meshgrid(ru, zu, indexing="xy")
    V = np.full(R.shape, np.nan)
    r_index = {float(x): i for i, x in enumerate(ru)}
    z_index = {float(x): i for i, x in enumerate(zu)}
    for ri, zi, vi in zip(r, z, v):
        V[z_index[float(np.round(zi, 12))], r_index[float(np.round(ri, 12))]] = vi
    return ru, zu, V


def _corr_rmse(a, b):
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < 8:
        return None, None
    aa, bb = a[m], b[m]
    aa = aa - aa.mean()
    bb = bb - bb.mean()
    denom = np.linalg.norm(aa) * np.linalg.norm(bb)
    corr = float(np.dot(aa, bb) / denom) if denom > 0 else float("nan")
    rmse = float(np.sqrt(np.mean((a[m] - b[m]) ** 2)))
    return corr, rmse


def _interp_seed_to(seed, r_t, z_t, unit_box=False, R=0.25, H=0.21):
    """Nearest bilinear via matplotlib-free griddata-ish: remap seed axes then sample."""
    sr, sz, sv = seed["r"], seed["z"], seed["v"]
    if unit_box or (sr.max() <= 1.0 + 1e-6 and sz.max() <= 1.0 + 1e-6):
        # target physical → unit
        rn = np.clip(r_t / max(R, 1e-12), 0, 1)
        zn = np.clip(z_t / max(H, 1e-12), 0, 1)
        # build regular seed grid in unit
        ru, zu, V = _grid_from_scatter(sr, sz, sv)
        if ru is None:
            return None
        # map target to indices
        out = np.empty_like(rn)
        for i, (x, y) in enumerate(zip(rn, zn)):
            # bilinear
            ir = np.searchsorted(ru, x)
            iz = np.searchsorted(zu, y)
            ir1 = min(max(ir, 1), len(ru) - 1)
            iz1 = min(max(iz, 1), len(zu) - 1)
            ir0, iz0 = ir1 - 1, iz1 - 1
            r0, r1 = ru[ir0], ru[ir1]
            z0, z1 = zu[iz0], zu[iz1]
            tr = (x - r0) / (r1 - r0) if r1 > r0 else 0.0
            tz = (y - z0) / (z1 - z0) if z1 > z0 else 0.0
            out[i] = (
                (1 - tr) * (1 - tz) * V[iz0, ir0]
                + tr * (1 - tz) * V[iz0, ir1]
                + (1 - tr) * tz * V[iz1, ir0]
                + tr * tz * V[iz1, ir1]
            )
        return out
    # physical seed: same bilinear on physical axes
    ru, zu, V = _grid_from_scatter(sr, sz, sv)
    if ru is None:
        return None
    out = np.empty_like(r_t)
    for i, (x, y) in enumerate(zip(r_t, z_t)):
        x = min(max(x, ru[0]), ru[-1])
        y = min(max(y, zu[0]), zu[-1])
        ir = np.searchsorted(ru, x)
        iz = np.searchsorted(zu, y)
        ir1 = min(max(ir, 1), len(ru) - 1)
        iz1 = min(max(iz, 1), len(zu) - 1)
        ir0, iz0 = ir1 - 1, iz1 - 1
        r0, r1 = ru[ir0], ru[ir1]
        z0, z1 = zu[iz0], zu[iz1]
        tr = (x - r0) / (r1 - r0) if r1 > r0 else 0.0
        tz = (y - z0) / (z1 - z0) if z1 > z0 else 0.0
        out[i] = (
            (1 - tr) * (1 - tz) * V[iz0, ir0]
            + tr * (1 - tz) * V[iz0, ir1]
            + (1 - tr) * tz * V[iz1, ir0]
            + tr * tz * V[iz1, ir1]
        )
    return out


def _parse_time_from_name(p: Path):
    m = re.search(r"french_T_rz_t([0-9]+\.[0-9]+)\.csv$", p.name)
    return float(m.group(1)) if m else None


def _pick_snap(case_dir: Path, t_want: float):
    snaps = list((case_dir / "trz_snaps").glob("french_T_rz_t*.csv"))
    snaps += list(case_dir.glob("french_T_rz_t*.csv"))
    if not snaps:
        # fall back to french_T_rz.csv as "end"
        end = case_dir / "french_T_rz.csv"
        return end if end.exists() else None, None
    best = None
    best_dt = 1e99
    best_t = None
    for p in snaps:
        t = _parse_time_from_name(p)
        if t is None:
            continue
        dt = abs(t - t_want)
        if dt < best_dt:
            best_dt, best, best_t = dt, p, t
    return best, best_t


def _monitor_tail(case_dir: Path):
    mon = case_dir / "french_stirring_em_monitor.csv"
    if not mon.exists():
        return {}
    rows = list(csv.DictReader(open(mon)))
    if not rows:
        return {}
    last = rows[-1]
    out = {}
    for k in (
        "physical_time_s",
        "revolutions",
        "U_max",
        "U_rms",
        "Np_sph",
        "T_mean",
        "T_max",
        "wall_clock_h",
    ):
        if k in last:
            try:
                out[k] = float(last[k])
            except Exception:
                pass
    return out


def _case_row(case_dir: Path, seed_t, seed_q, times, R_melt=0.25, H_melt=0.21):
    name = case_dir.name
    rc = int((case_dir / "exit_code.txt").read_text().strip()) if (case_dir / "exit_code.txt").exists() else 1
    meta = {}
    for key, fname in (("rpm", "rpm.txt"), ("P_abs", "p_abs.txt"), ("skull_delta", "skull_delta.txt")):
        p = case_dir / fname
        if p.exists():
            try:
                meta[key] = float(p.read_text().strip())
            except Exception:
                pass
    mon = _monitor_tail(case_dir)
    q_path = case_dir / "french_q_rz.csv"
    q_field = _load_rz_csv(q_path, value_keys=("Q_Wpm3", "Q", "Q_norm", "value"))
    q_corr = q_rmse = None
    if q_field and seed_q is not None:
        q_seed_on = _interp_seed_to(
            seed_q, q_field["r"], q_field["z"], unit_box=True, R=R_melt, H=H_melt
        )
        if q_seed_on is not None:
            # compare shapes after normalizing each to mean 1
            a = q_field["v"]
            b = q_seed_on
            if np.nanmean(a) > 0 and np.nanmean(b) > 0:
                a = a / np.nanmean(a)
                b = b / np.nanmean(b)
            q_corr, q_rmse = _corr_rmse(a, b)

    t_metrics = {}
    for tw in times:
        snap, t_got = _pick_snap(case_dir, tw)
        t_field = _load_rz_csv(snap, value_keys=("T", "T_K", "Temperature", "value")) if snap else None
        corr = rmse = None
        if t_field and seed_t is not None:
            # Fig8a seed is physical T; sim T is physical
            t_seed_on = _interp_seed_to(
                seed_t, t_field["r"], t_field["z"], unit_box=False, R=R_melt, H=H_melt
            )
            if t_seed_on is not None:
                corr, rmse = _corr_rmse(t_field["v"], t_seed_on)
        t_metrics[tw] = dict(t_got=t_got, corr=corr, rmse=rmse, snap=str(snap) if snap else None)

    ok = rc == 0 and mon.get("physical_time_s", 0) >= 0.9 * max(times) * 0.5  # soft
    return dict(
        name=name,
        ok=ok,
        exit=rc,
        meta=meta,
        mon=mon,
        q_corr=q_corr,
        q_rmse=q_rmse,
        t_metrics=t_metrics,
    )


def _plot_panel(rows, times, out_png: Path, kind: str):
    if plt is None:
        return
    groups = {
        "SW": [r for r in rows if r["name"].startswith("SW_")],
        "EM": [r for r in rows if r["name"].startswith("EM_")],
        "ST": [r for r in rows if r["name"].startswith("ST_")],
    }
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.2), sharey=True)
    for ax, (gname, grows) in zip(axes, groups.items()):
        for r in grows:
            xs, ys = [], []
            for tw in times:
                m = r["t_metrics"].get(tw, {})
                if m.get("corr") is not None:
                    xs.append(tw)
                    ys.append(m["corr"])
            if xs:
                ax.plot(xs, ys, "o-", label=r["name"])
        ax.set_title(gname)
        ax.set_xlabel("t [s]")
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=7)
    axes[0].set_ylabel("corr(T, Fig8a seed)")
    fig.suptitle(f"Phase5 soft T corr vs Fig8a ({kind}) — ASSUMED", fontsize=10)
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=140)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--seed-t", type=Path, required=True)
    ap.add_argument("--seed-q", type=Path, required=True)
    ap.add_argument("--times", type=str, default="20,40,60")
    ap.add_argument("--R", type=float, default=0.25)
    ap.add_argument("--H", type=float, default=0.21)
    args = ap.parse_args()
    times = [float(x) for x in args.times.split(",") if x.strip()]

    seed_t = _load_rz_csv(args.seed_t, value_keys=("T", "T_K", "value_normalized", "value"))
    seed_q = _load_rz_csv(args.seed_q, value_keys=("Q_norm", "Q_Wpm3", "Q", "value"))
    placeholder = False
    meta_q = args.seed_q.with_name(args.seed_q.name.replace("_rz.csv", ".meta.json"))
    if not meta_q.exists():
        meta_q = args.seed_q.with_suffix(".meta.json")
    # conventional name
    meta_q2 = args.seed_q.parent / "SEED_STIR_Fig10_Qshape_ASSUMED.meta.json"
    for mp in (meta_q, meta_q2):
        if mp.exists():
            try:
                placeholder = bool(json.loads(mp.read_text()).get("PLACEHOLDER_NOT_FIG10", False))
            except Exception:
                pass

    rows = []
    for name in CASE_ORDER:
        d = args.root / name
        if d.is_dir():
            rows.append(_case_row(d, seed_t, seed_q, times, R_melt=args.R, H_melt=args.H))
    # also any other case dirs with exit_code
    known = {r["name"] for r in rows}
    for d in sorted(args.root.iterdir()):
        if d.is_dir() and d.name not in known and (d / "exit_code.txt").exists():
            rows.append(_case_row(d, seed_t, seed_q, times, R_melt=args.R, H_melt=args.H))

    fig_dir = args.root / "figures_compare"
    _plot_panel(rows, times, fig_dir / "phase5_T_corr_vs_Fig8a.png", "matrix")

    lines = []
    lines.append("# Phase5 French stir matrix — MATRIX")
    lines.append("")
    lines.append(f"- root: `{args.root}`")
    lines.append(f"- times: {times}")
    lines.append(f"- seed T: `{args.seed_t.name}`")
    lines.append(f"- seed Q: `{args.seed_q.name}`" + (" **PLACEHOLDER_NOT_FIG10**" if placeholder else ""))
    lines.append("- identity: ASSUMED engineering compare; **not** literature PASS")
    lines.append("")
    lines.append("| case | exit | t_end | rev | Np | T_mean | Q_shape_corr | Tcorr@20 | Tcorr@40 | Tcorr@60 |")
    lines.append("|------|------|-------|-----|----|--------|--------------|----------|----------|----------|")
    for r in rows:
        mon = r["mon"]
        tm = r["t_metrics"]

        def fmtc(tw):
            c = tm.get(tw, {}).get("corr")
            return f"{c:.3f}" if c is not None and math.isfinite(c) else "—"

        lines.append(
            "| {name} | {ex} | {te} | {rev} | {np} | {tmn} | {qc} | {c20} | {c40} | {c60} |".format(
                name=r["name"],
                ex=r["exit"],
                te=f"{mon.get('physical_time_s', float('nan')):.2f}" if mon else "—",
                rev=f"{mon.get('revolutions', float('nan')):.2f}" if mon else "—",
                np=f"{mon.get('Np_sph', float('nan')):.3g}" if mon else "—",
                tmn=f"{mon.get('T_mean', float('nan')):.1f}" if mon else "—",
                qc=f"{r['q_corr']:.3f}" if r["q_corr"] is not None else "—",
                c20=fmtc(20.0) if 20.0 in times else fmtc(times[0]),
                c40=fmtc(40.0) if 40.0 in times else (fmtc(times[1]) if len(times) > 1 else "—"),
                c60=fmtc(60.0) if 60.0 in times else (fmtc(times[-1]) if times else "—"),
            )
        )
    lines.append("")
    lines.append("## Notes")
    lines.append("")
    lines.append("- SW/EM/ST are one-at-a-time sweeps (others at baseline 10 rpm / 60 kW / skull=5 mm).")
    lines.append("- Q corr is shape-only (mean-normalized) vs Fig.10 seed; with PLACEHOLDER seed, treat as smoke metric.")
    lines.append("- Figures: `figures_compare/phase5_T_corr_vs_Fig8a.png`")
    md = args.root / "MATRIX.md"
    md.write_text("\n".join(lines) + "\n")
    summary = {"root": str(args.root), "placeholder_fig10": placeholder, "rows": rows}
    (args.root / "matrix_summary.json").write_text(
        json.dumps(summary, indent=2, default=lambda o: None if isinstance(o, float) and not math.isfinite(o) else o)
        + "\n"
    )
    print(f"[phase5-post] wrote {md} cases={len(rows)} placeholder_fig10={int(placeholder)}")


if __name__ == "__main__":
    main()
