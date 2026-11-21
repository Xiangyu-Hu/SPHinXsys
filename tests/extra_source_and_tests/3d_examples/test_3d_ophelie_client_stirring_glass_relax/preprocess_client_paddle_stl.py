#!/usr/bin/env python3
"""Convert client stirring paddle STL (CAD mm, +Y up) to solver metres (+Z up).

Translate crucible outer-bottom axis to origin using glass AABB from Doc230 audit,
then (x,y,z)_cad -> (x,-z,y), then mm->m.
"""
from __future__ import annotations

import argparse
import struct
from pathlib import Path

# From makeClientStirringAssemblyCatalog / client audit
GLASS_MIN = (32.739, 3.732, 1029.0)
GLASS_MAX = (317.739, 253.732, 1314.0)


def origin_cad_mm() -> tuple[float, float, float]:
    return (
        0.5 * (GLASS_MIN[0] + GLASS_MAX[0]),
        GLASS_MIN[1],
        0.5 * (GLASS_MIN[2] + GLASS_MAX[2]),
    )


def cad_mm_to_solver_m(p: tuple[float, float, float], origin: tuple[float, float, float]) -> tuple[float, float, float]:
    x = p[0] - origin[0]
    y = p[1] - origin[1]
    z = p[2] - origin[2]
    return (x * 1.0e-3, -z * 1.0e-3, y * 1.0e-3)


def cad_n_to_solver(n: tuple[float, float, float]) -> tuple[float, float, float]:
    # Same linear map without translation/scale on direction, then renormalize.
    sx, sy, sz = (n[0], -n[2], n[1])
    norm = (sx * sx + sy * sy + sz * sz) ** 0.5
    if norm <= 0.0:
        return (0.0, 0.0, 1.0)
    return (sx / norm, sy / norm, sz / norm)


def convert(src: Path, dst: Path) -> None:
    origin = origin_cad_mm()
    with src.open("rb") as f:
        header = f.read(80)
        (count,) = struct.unpack("<I", f.read(4))
        facets = [struct.unpack("<12fH", f.read(50)) for _ in range(count)]

    lo = [float("inf")] * 3
    hi = [float("-inf")] * 3
    dst.parent.mkdir(parents=True, exist_ok=True)
    with dst.open("wb") as f:
        f.write(header[:80].ljust(80, b"\0"))
        f.write(struct.pack("<I", count))
        for v in facets:
            n = cad_n_to_solver(v[0:3])
            verts = [
                cad_mm_to_solver_m(v[3 + 3 * k : 6 + 3 * k], origin)  # type: ignore[index]
                for k in range(3)
            ]
            for p in verts:
                for d in range(3):
                    lo[d] = min(lo[d], p[d])
                    hi[d] = max(hi[d], p[d])
            f.write(struct.pack("<12fH", *n, *verts[0], *verts[1], *verts[2], v[12]))

    print(f"{src} -> {dst}")
    print(f"  triangles={count}")
    print(f"  origin_cad_mm={origin}")
    print(f"  aabb_solver_m min={tuple(lo)} max={tuple(hi)}")


def main() -> None:
    here = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--src",
        type=Path,
        default=here / "data" / "stirring paddle.STL",
    )
    ap.add_argument(
        "--dst",
        type=Path,
        default=here / "data" / "client_stir_paddle_solver_m.stl",
    )
    args = ap.parse_args()
    convert(args.src, args.dst)


if __name__ == "__main__":
    main()
