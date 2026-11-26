#!/usr/bin/env python3
"""Build client-NC hybrid T IC (thin upper-wall tongue + mid stratification)
and a Fig.6-mapped Q shape on the same (r,z) grid.

ASSUMED engineering maps — not paper Maxwell Q, not freeze-T.
"""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]  # doc248_client_cases
CFG = ROOT / "configs"
REF = ROOT.parent  # docs/ophelie/reference
FIG8_PARTIAL = (
    REF / "standard_data/fig7_fig8_literature/Doc204/reference/fig8_partial_linear_reconstruction.csv"
)
FIG6_CONTOURS = REF / "standard_data/doc247_french_digitization/data/NC_Fig6_contour_points.csv"

R_CLI = 0.1415
H_CLI = 0.230
Z0 = 0.001
NR, NZ = 81, 81
T_LO_C, T_HI_C = 830.0, 1050.0
# Client midradius/axis readouts (C): depths from free surface 2 / 7 / 12 cm
AXIS_DEPTH_M = np.array([0.02, 0.07, 0.12])
AXIS_T_C = np.array([817.54, 1068.82, 980.66])
P_ABS_W = 18000.0


def load_fig8_partial(path: Path):
    rows = list(csv.DictReader(path.open()))
    rs = sorted({float(r["r_over_R"]) for r in rows})
    zs = sorted({float(r["z_over_H"]) for r in rows})
    r_of = {v: i for i, v in enumerate(rs)}
    z_of = {v: i for i, v in enumerate(zs)}
    T = np.full((len(zs), len(rs)), np.nan)
    for r in rows:
        tv = r.get("T_over_Tmax") or ""
        if not tv:
            continue
        T[z_of[float(r["z_over_H"])], r_of[float(r["r_over_R"])]] = float(tv)
    return np.asarray(rs), np.asarray(zs), T


def fill_nan_nearest(T: np.ndarray) -> np.ndarray:
    """Fill NaNs by nearest finite neighbor (separable approx)."""
    out = T.copy()
    # horizontal then vertical sweeps
    for _ in range(4):
        for axis in (1, 0):
            n = out.shape[axis]
            for i in range(n):
                sl = [slice(None), slice(None)]
                sl[axis] = i
                col = out[tuple(sl)]
                mask = np.isnan(col)
                if not mask.any() or mask.all():
                    continue
                idx = np.where(~mask)[0]
                good = col[~mask]
                bad = np.where(mask)[0]
                # nearest in 1D index
                j = np.searchsorted(idx, bad)
                j = np.clip(j, 1, len(idx) - 1)
                left = idx[j - 1]
                right = idx[j]
                pick = np.where(np.abs(bad - left) <= np.abs(bad - right), left, right)
                col2 = col.copy()
                col2[bad] = col[pick]
                out[tuple(sl)] = col2
    # any remaining → 0.45 (cold floor under Fig.8 labelled range)
    out = np.where(np.isnan(out), 0.45, out)
    # Cap hot core: keep thin tongue — do not inflate unresolved >0.9 region
    out = np.clip(out, 0.45, 0.92)
    return out


def axis_profile_C(z_norm: np.ndarray) -> np.ndarray:
    """Piecewise axis T(C) from client 2/7/12 cm depths + ends."""
    z_pts = 1.0 - AXIS_DEPTH_M / H_CLI  # free-surface depth → z/H from bottom
    # order ascending z
    order = np.argsort(z_pts)
    zp = z_pts[order]
    tp = AXIS_T_C[order]
    # extend to bottom/top with gentle clamp
    zp_ext = np.concatenate([[0.0], zp, [1.0]])
    # bottom cooler than 12 cm, top near 2 cm but slightly warmer free surface band
    tp_ext = np.concatenate([[tp[0] - 80.0], tp, [tp[-1] + 20.0]])
    return np.interp(z_norm, zp_ext, tp_ext)


def nested_ellipse(RR, ZZ, rc, zc, a_wall, a_axis, b_up, b_down):
    """Nested ovals. Scales blend smoothly so contours stay round (no teardrop kink)."""
    dr = RR - rc
    dz = ZZ - zc
    # Logistic blend: C-infinity, no corner where the two half-axes meet.
    s_r = 1.0 / (1.0 + np.exp(dr / 0.12))
    s_z = 1.0 / (1.0 + np.exp(dz / 0.10))
    a = a_wall + (a_axis - a_wall) * s_r
    b = b_up + (b_down - b_up) * s_z
    return np.exp(-0.5 * ((dr / a) ** 2 + (dz / b) ** 2))


def _gauss(RR, ZZ, rc, zc, a, b):
    return np.exp(-0.5 * (((RR - rc) / a) ** 2 + ((ZZ - zc) / b) ** 2))


def _smooth_max(a, b, k=6.0):
    return np.log(np.exp(k * a) + np.exp(k * b)) / k


def tilted_round_lobe(RR, ZZ):
    """Two blended round ovals: wall core, plus a down-tilted arm to the axis.

    Lower-right of the wall oval is shortened smoothly so that corner tucks in.
    """
    # Round cap near the wall, small enough that the cold skin does not slice it flat.
    # Its lower-right therefore closes as a curve, inset from the corner.
    # Hot cap raised: client 7 cm (z/H≈0.70) sits on the lower side of the hot region.
    wall = _gauss(RR, ZZ, 0.82, 0.86, 0.13, 0.12)
    # Arm still reaches the axis near 7 cm, but the body of the tongue sits higher.
    arm = _gauss(RR, ZZ, 0.40, 0.76, 0.38, 0.14)
    return _smooth_max(wall, arm)


def sketch_lobe_mask(RR, ZZ):
    """v6d: tilted round oval + thin cold wall skin."""
    hot = tilted_round_lobe(RR, ZZ)
    skin = np.exp(-(np.maximum(RR - 0.955, 0.0) / 0.022) ** 2)
    top_open = np.clip((ZZ - 0.90) / 0.08, 0.0, 1.0)
    axis = np.exp(-((RR / 0.22) ** 2))
    lobe = hot * skin * (1.0 - 0.55 * top_open * axis)
    return lobe, hot


def build_hybrid_T():
    """v6: nested ellipses, cold wall skin, hot region reaches the axis, steep near wall."""
    r_n = np.linspace(0.0, 1.0, NR)
    z_n = np.linspace(0.0, 1.0, NZ)
    RR, ZZ = np.meshgrid(r_n, z_n)

    Tstar = 0.50 + 0.10 * ZZ * (1.0 - 0.65 * ZZ)
    lobe, _ = sketch_lobe_mask(RR, ZZ)
    Tstar = Tstar + 0.44 * lobe
    Tstar = np.clip(Tstar, 0.48, 0.96)

    T_C = T_LO_C + (Tstar - 0.48) / (0.93 - 0.48) * (T_HI_C - T_LO_C)
    # Pin the centerline to the client 2/7/12 cm readings. Narrow weight so it
    # stays the end of the tilted oval, not a separate hot island.
    T_ax = axis_profile_C(ZZ)
    w_axis = np.exp(-((RR / 0.14) ** 2))
    T_C = (1.0 - w_axis) * T_C + w_axis * T_ax
    T_C = np.clip(T_C, 750.0, 1100.0)
    T_K = T_C + 273.15
    r_m = RR * R_CLI
    z_m = Z0 + ZZ * H_CLI
    return r_n, z_n, r_m, z_m, T_K, T_C, Tstar


def reconstruct_fig6_Qstar(r_n, z_n):
    """Q nested like the T sketch but peaked a bit closer to the wall (EM from the sides).

    Still leaves a thin unheated skin at r→1 so the cold wall is not a heat source.
    """
    RR, ZZ = np.meshgrid(r_n, z_n)
    q = nested_ellipse(RR, ZZ, rc=0.92, zc=0.80, a_wall=0.045, a_axis=0.60, b_up=0.14, b_down=0.22)
    skin = np.clip((1.0 - RR) / 0.030, 0.0, 1.0)
    Q = q * skin
    Q = np.clip(Q, 0.0, None)
    # gentle blend with digitized Fig.6 IDW if file present (optional shape prior)
    try:
        rows = list(csv.DictReader(FIG6_CONTOURS.open()))
        pts = np.asarray(
            [(float(r["x_unit"]), float(r["z_or_y_unit"]), float(r["value_normalized"])) for r in rows]
        )
        Q6 = np.zeros_like(RR)
        for iz in range(NZ):
            for ir in range(NR):
                d2 = (pts[:, 0] - RR[iz, ir]) ** 2 + (pts[:, 1] - ZZ[iz, ir]) ** 2
                w = 1.0 / np.maximum(d2, 1e-8)
                Q6[iz, ir] = float(np.sum(w * pts[:, 2]) / np.sum(w))
        Q6 = Q6 / (np.max(Q6) + 1e-30)
        Q = 0.70 * Q + 0.30 * Q6 * Q
    except Exception:
        pass
    return np.clip(Q, 0.0, None)


def normalize_Q_Wpm3(Qstar, r_m, z_m, P):
    """Axisymmetric ∫ Q 2π r dr dz = P."""
    # cell-centered approx on structured grid
    r_1d = r_m[0, :]
    z_1d = z_m[:, 0]
    dr = np.gradient(r_1d)
    dz = np.gradient(z_1d)
    # volume element per cell
    dV = 2.0 * np.pi * r_m * dr[np.newaxis, :] * dz[:, np.newaxis]
    # axis r=0: use small r_eff
    dV[:, 0] = np.pi * (0.5 * dr[0]) ** 2 * dz
    S = float(np.sum(Qstar * dV))
    if S <= 0:
        raise RuntimeError("zero Fig6 shape integral")
    return Qstar * (P / S), S


def write_T_csv(path: Path, r_m, z_m, T_K):
    with path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["r_m", "z_m", "T_K"])
        for iz in range(z_m.shape[0]):
            for ir in range(r_m.shape[1]):
                w.writerow([f"{r_m[iz, ir]:.8g}", f"{z_m[iz, ir]:.8g}", f"{T_K[iz, ir]:.8g}"])


def write_Q_csv(path: Path, r_m, z_m, Q):
    with path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["r_m", "z_m", "Q_Wpm3", "Q_over_Qmax"])
        qmax = float(np.max(Q)) + 1e-30
        for iz in range(z_m.shape[0]):
            for ir in range(r_m.shape[1]):
                w.writerow(
                    [
                        f"{r_m[iz, ir]:.8g}",
                        f"{z_m[iz, ir]:.8g}",
                        f"{Q[iz, ir]:.8g}",
                        f"{(Q[iz, ir] / qmax):.8g}",
                    ]
                )


def preview_T(path, r_n, z_n, T_C, title):
    fig, ax = plt.subplots(figsize=(4.2, 5.0))
    span = float(np.nanmax(T_C) - np.nanmin(T_C)) + 1e-9
    Tn = (T_C - np.nanmin(T_C)) / span
    levels = np.linspace(0.0, 1.0, 11)
    cs = ax.contourf(r_n, z_n, Tn, levels=levels, cmap="inferno")
    ax.contour(r_n, z_n, Tn, levels=[0.35, 0.5, 0.65, 0.8, 0.92], colors="w", linewidths=0.7)
    ax.set_xlabel("r/R")
    ax.set_ylabel("z/H")
    ax.set_title(title, fontsize=9)
    ax.set_aspect("equal")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    fig.colorbar(cs, ax=ax, label="(T−Tmin)/(Tmax−Tmin)")
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def preview_Q(path, r_n, z_n, Qn, title):
    fig, ax = plt.subplots(figsize=(4.2, 5.0))
    levels = np.linspace(0.0, 1.0, 11)
    cs = ax.contourf(r_n, z_n, Qn, levels=levels, cmap="magma")
    ax.contour(r_n, z_n, Qn, levels=[0.1, 0.3, 0.9], colors="w", linewidths=0.7)
    ax.set_xlabel("r/R")
    ax.set_ylabel("z/H")
    ax.set_title(title)
    ax.set_aspect("equal")
    fig.colorbar(cs, ax=ax, label="Q / Qmax")
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def main():
    CFG.mkdir(parents=True, exist_ok=True)
    r_n, z_n, r_m, z_m, T_K, T_C, Tstar = build_hybrid_T()
    t_csv = CFG / "client_nc_hybrid_tongue_v6_H230_ASSUMED_T_rz.csv"
    write_T_csv(t_csv, r_m, z_m, T_K)
    preview_T(
        CFG / "client_nc_hybrid_tongue_v6_H230_preview_norm.png",
        r_n,
        z_n,
        T_C,
        "v6e: hot region raised; axis still 2/7/12 cm",
    )
    # synthetic top-view ring at free-surface band (axisym → annulus)
    fig, ax = plt.subplots(figsize=(4.2, 4.2))
    th = np.linspace(0, 2 * np.pi, 180)
    rr = np.linspace(0, 1, 80)
    TH, RR2 = np.meshgrid(th, rr)
    # sample T* near free surface z/H=0.96
    iz = int(np.argmin(np.abs(z_n - 0.96)))
    Trow = T_C[iz, :] / np.max(T_C)
    Ttop = np.interp(RR2.ravel(), r_n, Trow).reshape(RR2.shape)
    ax.pcolormesh(RR2 * np.cos(TH), RR2 * np.sin(TH), Ttop, shading="auto", cmap="inferno", vmin=0.5, vmax=1.0)
    ax.set_aspect("equal")
    ax.set_title("v5 top @ z/H≈0.96 (cool core + wall ring)", fontsize=9)
    ax.set_xlabel("x/R")
    ax.set_ylabel("y/R")
    fig.tight_layout()
    fig.savefig(CFG / "client_nc_hybrid_tongue_v6_H230_preview_top.png", dpi=140)
    plt.close(fig)

    Qstar = reconstruct_fig6_Qstar(r_n, z_n)
    Q, S = normalize_Q_Wpm3(Qstar, r_m, z_m, P_ABS_W)
    q_csv = CFG / "client_nc_fig6_mapped_H230_ASSUMED_Q_rz.csv"
    write_Q_csv(q_csv, r_m, z_m, Q)
    preview_Q(
        CFG / "client_nc_fig6_mapped_H230_preview_norm.png",
        r_n,
        z_n,
        Q / (np.max(Q) + 1e-30),
        "Fig.6→client Q shape (ASSUMED; not EM Maxwell)",
    )

    # power centroid
    r_1d = r_m[0, :]
    z_1d = z_m[:, 0]
    dr = np.gradient(r_1d)
    dz = np.gradient(z_1d)
    dV = 2.0 * np.pi * r_m * dr[np.newaxis, :] * dz[:, np.newaxis]
    dV[:, 0] = np.pi * (0.5 * dr[0]) ** 2 * dz
    P = float(np.sum(Q * dV))
    zQ = float(np.sum(Q * dV * z_m) / P)
    rQ = float(np.sum(Q * dV * r_m) / P)
    upper = float(np.sum(Q * dV * (z_m >= (Z0 + 0.5 * H_CLI))) / P)

    # probe samples on hybrid
    def sample_T(r_over, depth):
        zn = 1.0 - depth / H_CLI
        rn = r_over
        # nearest
        ir = int(np.argmin(np.abs(r_n - rn)))
        iz = int(np.argmin(np.abs(z_n - zn)))
        return float(T_C[iz, ir])

    probes = []
    for d, tc in zip(AXIS_DEPTH_M, AXIS_T_C):
        probes.append(
            {
                "id": f"axis_d{int(d*100)}cm",
                "r_over_R": 0.0,
                "depth_m": float(d),
                "T_hybrid_C": sample_T(0.0, d),
                "T_client_C": float(tc),
            }
        )

    meta = {
        "schema": "doc248-client-nc-hybrid-tongue-v6-nested-ellipse-2026-11-26",
        "status": "ASSUMED",
        "supersedes": "configs/client_nc_hybrid_tongue_v5_H230_ASSUMED_T_rz.csv",
        "evidence": {
            "client_topdown_photo": "cool center + wall hot ring; EM from sides",
            "user_sketch": "one connected lobe: fat at wall, taper to mid",
            "user_fix": "no orphan mid hot; top cools faster; hot must stay connected from wall",
        },
        "T_ic": {
            "csv": str(t_csv.relative_to(ROOT)),
            "construction": (
                "v6 nested ellipse: core inset (r/R≈0.68,z/H≈0.74); cold skin at r→1; "
                "long scale toward axis so outer contours reach centerline; "
                "short scale toward wall so isotherms pack up (steep gradient); top cools"
            ),
            "probe_check_C": probes,
            "T_C_band": [float(np.min(T_C)), float(np.max(T_C))],
            "preview_rz": "configs/client_nc_hybrid_tongue_v6_H230_preview_norm.png",
            "preview_top": "configs/client_nc_hybrid_tongue_v6_H230_preview_top.png",
        },
        "Q_shape": {
            "csv": str(q_csv.relative_to(ROOT)),
            "source": "doc247 NC_Fig6_contour_points.csv",
            "note": (
                "Shape only; frozen_q does not yet load Q_rz. Use as reference + "
                "EM levers: --doc246-q-upper-kappa and/or --client-nc-coil-z-frac-* after recompile. "
                "NOT paper Maxwell; NOT to be labelled as EM prediction."
            ),
            "P_abs_W": P_ABS_W,
            "integral_check_W": P,
            "shape_integral_before_scale": S,
            "power_centroid_r_m": rQ,
            "power_centroid_z_m": zQ,
            "power_centroid_z_over_H": (zQ - Z0) / H_CLI,
            "upper_half_power_frac": upper,
            "fig6_peak_native": {"r_over_R": [0.73, 0.86], "z_over_H": [0.73, 0.87]},
        },
        "geometry_m": {"R": R_CLI, "H": H_CLI, "z_bottom": Z0},
        "cli_T": f"--init-T-rz-csv={t_csv}",
    }
    meta_path = CFG / "client_nc_hybrid_tongue_v6_H230_ASSUMED.meta.json"
    meta_path.write_text(json.dumps(meta, indent=2) + "\n")
    print(json.dumps({"wrote_T": str(t_csv), "wrote_Q": str(q_csv), "meta": str(meta_path), "Q_centroid_zH": (zQ - Z0) / H_CLI, "upper_frac": upper}, indent=2))


if __name__ == "__main__":
    main()
