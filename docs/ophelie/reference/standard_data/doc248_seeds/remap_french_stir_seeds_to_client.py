#!/usr/bin/env python3
"""Remap French STIR Fig.8a T + Fig.10 Q seeds onto client stir geometry (ASSUMED).

French Fig.8a: R=0.25 m, H=0.21 m, physical r_m/z_m/T_K
French Fig.10: unit box r_over_R/z_over_H/Q_norm
Client STIR:   R=0.1415 m, H=0.230 m, z_bottom=0.001 m

Output physical grids (r_m,z_m) on the client melt AABB for --init-T-rz-csv / --init-Q-rz-csv.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

SRC = Path("/home/ubuntu/SPHinXsys")
IN_DIR = SRC / "docs/ophelie/搅拌caase起始"
OUT_DIR = SRC / "docs/ophelie/reference/standard_data/doc248_seeds"
FIG_DIR = OUT_DIR / "figures"

# French source geometry (Fig.8a meta)
R_FR = 0.25
H_FR = 0.21
# Client STIR (geometry_and_probe_mapping.json)
R_CLI = 0.1415
H_CLI = 0.230
Z_BOTTOM = 0.001
NR = NZ = 81


def load_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open() as f:
        rows = list(csv.DictReader(f))
    return list(rows[0].keys()), rows


def grid_from_physical_T(rows: list[dict[str, str]]):
    r = np.array([float(x["r_m"]) for x in rows])
    z = np.array([float(x["z_m"]) for x in rows])
    T = np.array([float(x["T_K"]) for x in rows])
    ru = np.unique(np.round(r, 12))
    zu = np.unique(np.round(z, 12))
    V = np.full((len(zu), len(ru)), np.nan)
    rmap = {float(x): i for i, x in enumerate(ru)}
    zmap = {float(x): i for i, x in enumerate(zu)}
    for ri, zi, vi in zip(r, z, T):
        V[zmap[float(np.round(zi, 12))], rmap[float(np.round(ri, 12))]] = vi
    return ru, zu, V


def grid_from_unit_Q(rows: list[dict[str, str]]):
    r = np.array([float(x["r_over_R"]) for x in rows])
    z = np.array([float(x["z_over_H"]) for x in rows])
    Q = np.array([float(x["Q_norm"]) for x in rows])
    ru = np.unique(np.round(r, 12))
    zu = np.unique(np.round(z, 12))
    V = np.full((len(zu), len(ru)), np.nan)
    rmap = {float(x): i for i, x in enumerate(ru)}
    zmap = {float(x): i for i, x in enumerate(zu)}
    for ri, zi, vi in zip(r, z, Q):
        V[zmap[float(np.round(zi, 12))], rmap[float(np.round(ri, 12))]] = vi
    return ru, zu, V


def bilinear(ru, zu, V, rq, zq):
    rq = float(np.clip(rq, ru[0], ru[-1]))
    zq = float(np.clip(zq, zu[0], zu[-1]))
    ir = int(np.searchsorted(ru, rq))
    iz = int(np.searchsorted(zu, zq))
    ir1 = min(max(ir, 1), len(ru) - 1)
    iz1 = min(max(iz, 1), len(zu) - 1)
    ir0, iz0 = ir1 - 1, iz1 - 1
    r0, r1 = ru[ir0], ru[ir1]
    z0, z1 = zu[iz0], zu[iz1]
    tr = (rq - r0) / (r1 - r0) if r1 > r0 else 0.0
    tz = (zq - z0) / (z1 - z0) if z1 > z0 else 0.0
    return (
        (1 - tr) * (1 - tz) * V[iz0, ir0]
        + tr * (1 - tz) * V[iz0, ir1]
        + (1 - tr) * tz * V[iz1, ir0]
        + tr * tz * V[iz1, ir1]
    )


def write_rz_csv(path: Path, r_axis, z_axis, field_name: str, values_iz_ir):
    with path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["r_m", "z_m", field_name])
        for iz, z in enumerate(z_axis):
            for ir, r in enumerate(r_axis):
                w.writerow([f"{r:.8f}", f"{z:.8f}", f"{float(values_iz_ir[iz, ir]):.8g}"])


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    t_path = IN_DIR / "SEED_STIR_Fig8a_axisym_ASSUMED_T_rz.csv"
    q_path = IN_DIR / "SEED_STIR_Fig10_Qshape_ASSUMED_rz.csv"
    _, t_rows = load_csv(t_path)
    _, q_rows = load_csv(q_path)
    ru_t, zu_t, VT = grid_from_physical_T(t_rows)
    ru_q, zu_q, VQ = grid_from_unit_Q(q_rows)

    # Client physical axes (absolute solver z)
    r_cli = np.linspace(0.0, R_CLI, NR)
    z_cli = np.linspace(Z_BOTTOM, Z_BOTTOM + H_CLI, NZ)
    T_cli = np.empty((NZ, NR))
    Q_cli = np.empty((NZ, NR))
    for iz, z in enumerate(z_cli):
        zn = (z - Z_BOTTOM) / H_CLI
        for ir, r in enumerate(r_cli):
            rn = r / R_CLI
            # Fig.8a source in physical FR meters
            T_cli[iz, ir] = bilinear(ru_t, zu_t, VT, rn * R_FR, zn * H_FR)
            # Fig.10 source already unit box
            Q_cli[iz, ir] = bilinear(ru_q, zu_q, VQ, rn, zn)

    out_t = OUT_DIR / "SEED_CLIENT_STIR_Fig8a_T_rz.csv"
    out_q = OUT_DIR / "SEED_CLIENT_STIR_Fig10_Q_rz.csv"
    write_rz_csv(out_t, r_cli, z_cli, "T_K", T_cli)
    write_rz_csv(out_q, r_cli, z_cli, "Q_norm", Q_cli)

    meta_t = {
        "schema": "RECONSTRUCTED_ASSUMED_INPUT",
        "ASSUMED": True,
        "NOT_literature_PASS": True,
        "source_csv": str(t_path),
        "source_figure": "Jacoutot CES 2008 STIR Fig.8a (axisym INITIAL_GUESS)",
        "source_R_m": R_FR,
        "source_H_m": H_FR,
        "client_R_m": R_CLI,
        "client_H_m": H_CLI,
        "client_z_bottom_m": Z_BOTTOM,
        "mapping": "r=(r_fr/R_fr)*R_cli; z=z_bottom+(z_fr/H_fr)*H_cli; absolute Kelvin kept",
        "nr": NR,
        "nz": NZ,
        "Tmin": float(np.nanmin(T_cli)),
        "Tmax": float(np.nanmax(T_cli)),
        "Tmean": float(np.nanmean(T_cli)),
        "use": "--init-T-rz-csv=.../SEED_CLIENT_STIR_Fig8a_T_rz.csv",
        "radial_origin": "crucible axis (not eccentric paddle axis)",
    }
    meta_q = {
        "schema": "RECONSTRUCTED_ASSUMED_INPUT",
        "ASSUMED": True,
        "NOT_literature_PASS": True,
        "NOT_EM_prediction": True,
        "source_csv": str(q_path),
        "source_figure": "Jacoutot CES 2008 STIR Fig.10 Q shape",
        "source_unit_box": True,
        "client_R_m": R_CLI,
        "client_H_m": H_CLI,
        "client_z_bottom_m": Z_BOTTOM,
        "mapping": "r=r_over_R*R_cli; z=z_bottom+z_over_H*H_cli; runtime scale ∫Q dV → --heat-power",
        "nr": NR,
        "nz": NZ,
        "Q_norm_min": float(np.nanmin(Q_cli)),
        "Q_norm_max": float(np.nanmax(Q_cli)),
        "Q_norm_mean": float(np.nanmean(Q_cli)),
        "use": "--init-Q-rz-csv=.../SEED_CLIENT_STIR_Fig10_Q_rz.csv with --heat-power=",
        "radial_origin": "crucible axis (not eccentric paddle axis)",
    }
    (OUT_DIR / "SEED_CLIENT_STIR_Fig8a_T.meta.json").write_text(json.dumps(meta_t, indent=2) + "\n")
    (OUT_DIR / "SEED_CLIENT_STIR_Fig10_Q.meta.json").write_text(json.dumps(meta_q, indent=2) + "\n")

    # QA figures
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(9.5, 8.0), constrained_layout=True)
    Rmg, Zmg = np.meshgrid(ru_t, zu_t, indexing="xy")
    pcm = axes[0, 0].pcolormesh(Rmg, Zmg, VT, shading="auto", cmap="inferno")
    axes[0, 0].set_title("French Fig.8a T (source)")
    axes[0, 0].set_xlabel("r [m]")
    axes[0, 0].set_ylabel("z [m]")
    axes[0, 0].set_aspect("equal")
    fig.colorbar(pcm, ax=axes[0, 0], fraction=0.046, pad=0.04, label="T [K]")

    Rmg, Zmg = np.meshgrid(r_cli, z_cli, indexing="xy")
    pcm = axes[0, 1].pcolormesh(Rmg, Zmg, T_cli, shading="auto", cmap="inferno")
    axes[0, 1].set_title(f"Client remap T  R={R_CLI} H={H_CLI}")
    axes[0, 1].set_xlabel("r [m]")
    axes[0, 1].set_ylabel("z [m]")
    axes[0, 1].set_aspect("equal")
    fig.colorbar(pcm, ax=axes[0, 1], fraction=0.046, pad=0.04, label="T [K]")

    Rmg, Zmg = np.meshgrid(ru_q, zu_q, indexing="xy")
    pcm = axes[1, 0].pcolormesh(Rmg, Zmg, VQ, shading="auto", cmap="magma")
    axes[1, 0].set_title("French Fig.10 Q_norm (source unit box)")
    axes[1, 0].set_xlabel("r/R")
    axes[1, 0].set_ylabel("z/H")
    axes[1, 0].set_aspect("equal")
    fig.colorbar(pcm, ax=axes[1, 0], fraction=0.046, pad=0.04, label="Q_norm")

    Rmg, Zmg = np.meshgrid(r_cli, z_cli, indexing="xy")
    pcm = axes[1, 1].pcolormesh(Rmg, Zmg, Q_cli, shading="auto", cmap="magma")
    axes[1, 1].set_title("Client remap Q_norm")
    axes[1, 1].set_xlabel("r [m]")
    axes[1, 1].set_ylabel("z [m]")
    axes[1, 1].set_aspect("equal")
    fig.colorbar(pcm, ax=axes[1, 1], fraction=0.046, pad=0.04, label="Q_norm")
    fig.suptitle("ASSUMED client STIR IC remap (NOT literature PASS)", fontsize=12)
    qa = FIG_DIR / "CLIENT_STIR_Fig8a_Fig10_remap_QA.png"
    fig.savefig(qa, dpi=150)

    print("wrote", out_t)
    print("wrote", out_q)
    print("wrote", qa)
    print("T range", meta_t["Tmin"], meta_t["Tmax"], "r_max", r_cli[-1], "z", z_cli[0], z_cli[-1])
    print("Q range", meta_q["Q_norm_min"], meta_q["Q_norm_max"])


if __name__ == "__main__":
    main()
