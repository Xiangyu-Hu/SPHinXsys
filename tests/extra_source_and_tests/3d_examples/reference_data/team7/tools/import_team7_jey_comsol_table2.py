#!/usr/bin/env python3
"""Convert COMSOL multiturn_coil_asymmetric_conductor_table2.txt to TEAM7 Jey reference CSV.

COMSOL table columns (no header):
  x [mm], Jy(x,72,19) at 50Hz, Jy(x,72,19) at 200Hz

Values use sign(real)*abs convention (scalar snapshot per frequency), not full phasor.
For L2 validation we map:
  phase0 = value at 50 Hz (or 200 Hz column when building that frequency block)
  phase90 = 0 (until quadrature measured data is available)

Usage:
  python3 import_team7_jey_comsol_table2.py /path/to/multiturn_coil_asymmetric_conductor_table2.txt \\
      -o ../TEAM7_Jey_A1_surface_reference_Am2.csv
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_comsol_table2(path: Path) -> list[tuple[float, float, float]]:
    rows: list[tuple[float, float, float]] = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = [p.strip() for p in line.replace("\t", ",").split(",") if p.strip()]
        if len(parts) < 3:
            continue
        rows.append((float(parts[0]), float(parts[1]), float(parts[2])))
    return rows


def write_team7_jey_csv(rows: list[tuple[float, float, float]], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "No",
                "x_mm",
                "Jy_50Hz_phase0_Am2",
                "Jy_50Hz_phase90_Am2",
                "Jy_200Hz_phase0_Am2",
                "Jy_200Hz_phase90_Am2",
            ]
        )
        f.write(
            "# Converted from COMSOL table2; phase90 columns set to 0 (scalar sign*abs snapshots).\n"
        )
        for i, (x_mm, jy50, jy200) in enumerate(rows, start=1):
            writer.writerow([i, x_mm, jy50, 0, jy200, 0])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("table2_txt", type=Path)
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "TEAM7_Jey_A1_surface_reference_Am2.csv",
    )
    args = parser.parse_args()
    rows = read_comsol_table2(args.table2_txt)
    if not rows:
        raise SystemExit(f"No data rows in {args.table2_txt}")
    write_team7_jey_csv(rows, args.output)
    print(f"Wrote {len(rows)} probes to {args.output}")


if __name__ == "__main__":
    main()
