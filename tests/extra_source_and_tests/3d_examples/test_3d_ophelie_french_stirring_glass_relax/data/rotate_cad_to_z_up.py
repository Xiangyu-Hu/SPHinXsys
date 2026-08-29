#!/usr/bin/env python3
"""Rotate the y-axis stirring CAD into the z-up frame used by the ophelie French tooling.

The shared French EM/thermal helpers (coil stacking, cylinder level sets, Robin/radiation
face tagging) all assume a cylinder whose axis is +z with the free surface on top. The CAD
delivered for this case is a cylinder about -y with the free surface at y = 0, so the STLs
are rewritten once here instead of generalising the shared headers.

Mapping (proper rotation of +90 deg about x, then a shift so the melt bottom sits at z = 0):

    X =  x_cad
    Y = -z_cad
    Z =  y_cad + H          with H = 210 mm, the melt height

Result: glass occupies z in [0, 210] mm, free surface at z = 210 mm, axis = +z.
"""
import struct
import sys
from pathlib import Path

MELT_HEIGHT_MM = 210.0

DATA_DIR = Path(__file__).resolve().parent
CONVERSIONS = [
    ("glass-y.stl", "glass-z.stl"),
    ("stirring_paddle_2.stl", "stirring_paddle_2_z.stl"),
]


def rotate(v):
    x, y, z = v
    return (x, -z, y + MELT_HEIGHT_MM)


def rotate_normal(v):
    x, y, z = v
    return (x, -z, y)


def convert(src: Path, dst: Path):
    with src.open("rb") as f:
        header = f.read(80)
        count = struct.unpack("<I", f.read(4))[0]
        facets = [struct.unpack("<12fH", f.read(50)) for _ in range(count)]

    lo = [float("inf")] * 3
    hi = [float("-inf")] * 3
    with dst.open("wb") as f:
        f.write(header[:80].ljust(80, b"\0"))
        f.write(struct.pack("<I", count))
        for v in facets:
            n = rotate_normal(v[0:3])
            verts = [rotate(v[3 + 3 * k: 6 + 3 * k]) for k in range(3)]
            for p in verts:
                for d in range(3):
                    lo[d] = min(lo[d], p[d])
                    hi[d] = max(hi[d], p[d])
            f.write(struct.pack("<12fH", *n, *verts[0], *verts[1], *verts[2], v[12]))

    print(f"{src.name} -> {dst.name}: {count} triangles")
    print(f"  bbox min (mm): ({lo[0]:.3f}, {lo[1]:.3f}, {lo[2]:.3f})")
    print(f"  bbox max (mm): ({hi[0]:.3f}, {hi[1]:.3f}, {hi[2]:.3f})")


def main():
    for src_name, dst_name in CONVERSIONS:
        src = DATA_DIR / src_name
        if not src.exists():
            print(f"missing {src}", file=sys.stderr)
            return 1
        convert(src, DATA_DIR / dst_name)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
