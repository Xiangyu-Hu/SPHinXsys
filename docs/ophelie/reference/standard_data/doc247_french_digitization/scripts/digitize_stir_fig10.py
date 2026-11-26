#!/usr/bin/env python3
"""Digitize Jacoutot et al. CES 2008 Fig.10 (axisymmetric Q after coupling).

Manual visual contour traces on the cropped reference frame → unit-box seed
`SEED_STIR_Fig10_Qshape_ASSUMED_rz.csv` for `--fig10-q-correct`.

ASSUMED reconstruction (not author numerical field). Regenerates contour CSV,
QA overlays, and seed + meta (clears PLACEHOLDER_NOT_FIG10).
"""
from __future__ import annotations

import csv
import json
import shutil
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from scipy.interpolate import griddata
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
SEED_DIR = ROOT.parent / "doc247_seeds"
FIGDIR = ROOT / "figures"
DATADIR = ROOT / "data"
FIGDIR.mkdir(parents=True, exist_ok=True)
DATADIR.mkdir(parents=True, exist_ok=True)

# Plot interior of STIR_Fig10_reference.png (symmetry axis left, wall right,
# free surface top, bottom at z=0).
L, T, R, B = 28, 8, 616, 542

# Contour levels from literature Fig.10 labels; pixel coords on reference crop.
TRACES = [
    (
        0.9,
        [
            (598, 250),
            (595, 230),
            (590, 215),
            (580, 205),
            (568, 205),
            (558, 215),
            (552, 230),
            (550, 250),
            (552, 270),
            (558, 288),
            (568, 298),
            (580, 300),
            (590, 292),
            (595, 275),
            (598, 250),
        ],
    ),
    (
        0.7,
        [
            (608, 255),
            (605, 220),
            (598, 195),
            (585, 178),
            (568, 170),
            (548, 172),
            (530, 185),
            (520, 205),
            (515, 230),
            (512, 255),
            (515, 280),
            (522, 305),
            (535, 325),
            (555, 335),
            (575, 335),
            (590, 325),
            (600, 305),
            (606, 280),
            (608, 255),
        ],
    ),
    (
        0.5,
        [
            (612, 260),
            (610, 210),
            (602, 175),
            (585, 150),
            (560, 138),
            (530, 140),
            (500, 155),
            (480, 180),
            (468, 210),
            (462, 245),
            (462, 280),
            (470, 315),
            (485, 345),
            (510, 365),
            (540, 375),
            (570, 372),
            (590, 355),
            (602, 330),
            (608, 300),
            (612, 260),
        ],
    ),
    (
        0.3,
        [
            (614, 265),
            (612, 200),
            (602, 155),
            (580, 125),
            (545, 110),
            (500, 112),
            (455, 130),
            (420, 160),
            (400, 200),
            (390, 245),
            (390, 290),
            (400, 335),
            (420, 375),
            (455, 405),
            (500, 420),
            (545, 422),
            (580, 405),
            (600, 370),
            (610, 325),
            (614, 265),
        ],
    ),
    (
        0.1,
        [
            (616, 270),
            (614, 185),
            (600, 130),
            (570, 95),
            (520, 78),
            (460, 80),
            (400, 100),
            (350, 135),
            (320, 180),
            (300, 230),
            (295, 285),
            (305, 340),
            (330, 390),
            (370, 430),
            (420, 455),
            (480, 468),
            (540, 465),
            (580, 440),
            (600, 395),
            (610, 340),
            (616, 270),
        ],
    ),
]

QMAX = 6.64e6
NR, NZ = 81, 81


def main() -> None:
    ref = FIGDIR / "STIR_Fig10_reference.png"
    if not ref.exists():
        raise SystemExit(f"missing {ref} — crop Fig.10 frame from PDF first")
    arr = np.array(Image.open(ref))

    fig, ax = plt.subplots(figsize=(7, 6))
    ax.imshow(arr)
    for val, pts in TRACES:
        p = np.array(pts)
        ax.plot(p[:, 0], p[:, 1], ".-", lw=1.0, ms=3, label=f"{val:g}")
    ax.legend(loc="upper left", fontsize=8)
    ax.set_title("STIR Fig.10 manual contour traces (ASSUMED digitization)")
    fig.tight_layout()
    qa = FIGDIR / "STIR_Fig10_trace_QA.png"
    fig.savefig(qa, dpi=160)
    plt.close(fig)

    rows = []
    for ci, (val, pts) in enumerate(TRACES):
        for order, (px, py) in enumerate(pts):
            rows.append(
                dict(
                    curve_id=ci,
                    point_order=order,
                    pixel_x=px,
                    pixel_y=py,
                    x_unit=(px - L) / (R - L),
                    z_or_y_unit=(B - py) / (B - T),
                    value_normalized=val,
                    Q_Wpm3=val * QMAX,
                    source_type="manual_visual_contour_trace",
                    coordinate_uncertainty_pixels=5,
                    value_uncertainty_estimate=0.05,
                )
            )
    cont_path = DATADIR / "STIR_Fig10_contour_points.csv"
    with cont_path.open("w", newline="") as f:
        wri = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        wri.writeheader()
        wri.writerows(rows)

    pts = np.array([[r["x_unit"], r["z_or_y_unit"]] for r in rows])
    v = np.array([r["value_normalized"] for r in rows])
    ru = np.linspace(0, 1, NR)
    zu = np.linspace(0, 1, NZ)
    RR, ZZ = np.meshgrid(ru, zu, indexing="xy")
    lin = griddata(pts, v, (RR, ZZ), method="linear")
    extrap = np.isnan(lin)
    near = griddata(pts, v, (RR, ZZ), method="nearest")
    dist = cKDTree(pts).query(np.column_stack([RR.ravel(), ZZ.ravel()]))[0].reshape(RR.shape)
    field = np.where(extrap, near * np.exp(-dist / 0.12), lin)
    field = np.clip(field, 0.0, 1.0)
    field[:, : max(1, NR // 5)] = np.minimum(field[:, : max(1, NR // 5)], 0.05)

    seed_path = SEED_DIR / "SEED_STIR_Fig10_Qshape_ASSUMED_rz.csv"
    with seed_path.open("w", newline="") as f:
        wri = csv.DictWriter(
            f, fieldnames=["r_over_R", "z_over_H", "Q_norm", "nearest_extrapolation"]
        )
        wri.writeheader()
        for iz in range(NZ):
            for ir in range(NR):
                wri.writerow(
                    dict(
                        r_over_R=f"{ru[ir]:.8f}",
                        z_over_H=f"{zu[iz]:.8f}",
                        Q_norm=f"{field[iz, ir]:.8f}",
                        nearest_extrapolation=int(extrap[iz, ir]),
                    )
                )

    meta = {
        "schema": "RECONSTRUCTED_ASSUMED_INPUT",
        "PLACEHOLDER_NOT_FIG10": False,
        "source": "Manual visual digitization of Jacoutot et al. CES 2008 Fig.10",
        "pdf": "docs/ophelie/Numerical modeling of coupled phenomena in a mechanically stirred molten-glass bath heated by induction.pdf",
        "figure": "Fig.10 dimensionless axisymmetric power source after coupling, t=1000s",
        "literature_Qmax_Wpm3": QMAX,
        "literature_Ptot_W": 60000.0,
        "literature_N_rpm": 10.0,
        "unit_box": True,
        "axes": "r_over_R from symmetry axis to wall; z_over_H from bottom to free surface",
        "contour_levels_digitized": [0.1, 0.3, 0.5, 0.7, 0.9],
        "grid": {"nr": NR, "nz": NZ},
        "NOT_EM_prediction": True,
        "ASSUMED": True,
        "note": "Sparse contour reconstruction + linear/nearest fill with axis-side decay prior. Replaces NC-Fig6 PLACEHOLDER.",
        "regenerate": "python3 docs/ophelie/reference/standard_data/doc247_french_digitization/scripts/digitize_stir_fig10.py",
        "contour_csv": str(cont_path.relative_to(SEED_DIR.parents[2])),
        "qa_overlay": str(qa.relative_to(SEED_DIR.parents[2])),
        "use": "--init-Q-rz-csv=...SEED_STIR_Fig10_Qshape_ASSUMED_rz.csv --fig10-q-correct --fig10-q-correct-blend=0.8",
    }
    (SEED_DIR / "SEED_STIR_Fig10_Qshape_ASSUMED.meta.json").write_text(
        json.dumps(meta, indent=2) + "\n"
    )

    fig, ax = plt.subplots(figsize=(5.2, 4.6))
    pcm = ax.pcolormesh(RR, ZZ, field, shading="auto", cmap="inferno", vmin=0, vmax=1)
    for val, pts_px in TRACES:
        p = np.array(pts_px, float)
        xu = (p[:, 0] - L) / (R - L)
        zu = (B - p[:, 1]) / (B - T)
        ax.plot(xu, zu, "c.--", lw=0.8, ms=2, alpha=0.9)
    ax.set_xlabel("r/R")
    ax.set_ylabel("z/H")
    ax.set_title("STIR Fig.10 seed Q/Qmax (ASSUMED digitization)")
    ax.set_aspect("equal")
    fig.colorbar(pcm, ax=ax, label="Q/Qmax")
    fig.tight_layout()
    outp = FIGDIR / "STIR_Fig10_INITIAL_GUESS.png"
    fig.savefig(outp, dpi=160)
    plt.close(fig)

    # refresh digitization_summary row if present
    summary = DATADIR / "digitization_summary.csv"
    if summary.exists():
        existing = list(csv.DictReader(summary.open()))
        existing = [r for r in existing if r.get("figure") != "STIR_Fig10"]
        existing.append(
            dict(
                figure="STIR_Fig10",
                points=str(len(rows)),
                curves=str(len(TRACES)),
                grid_points=str(NR * NZ),
                extrapolated_fraction=f"{float(extrap.mean()):.6f}",
            )
        )
        with summary.open("w", newline="") as f:
            wri = csv.DictWriter(
                f,
                fieldnames=[
                    "figure",
                    "points",
                    "curves",
                    "grid_points",
                    "extrapolated_fraction",
                ],
            )
            wri.writeheader()
            wri.writerows(existing)

    print(f"[fig10] seed={seed_path} points={len(rows)} Qnorm=[{field.min():.3f},{field.max():.3f}]")
    print(f"[fig10] QA={qa}")


if __name__ == "__main__":
    main()
