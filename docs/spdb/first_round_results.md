# SPDB first-round results (P0 + P1 + elastic smoke)

**Date:** 2026-10-15
**Branch:** `feature/spdb-thermo-viscoplastic` @ `4e686ee35`
**Case:** `tests/extra_source_and_tests/3d_examples/test_3d_spdb_geometry`（当时写在 `tests/3d_examples`，2026-10-18 已迁走）
**Material:** `LinearElasticSolid` **synthetic** — **NOT** 950 °C 20Cr13 prediction.

## Particle / volume table (dp = 1 mm)

| Body / mode | Particles | Theory V [m³] | Discrete V [m³] | Rel. error |
|-------------|-----------|---------------|-----------------|------------|
| Lower, with slot | 24280 | 2.41327e-05 | 2.428e-05 | 0.610% |
| Upper, with slot | 25280 | 2.51327e-05 | 2.528e-05 | 0.586% |
| Lower, no slot | 25280 | 2.51327e-05 | 2.528e-05 | 0.586% |
| Upper, no slot | 25280 | 2.51327e-05 | 2.528e-05 | 0.586% |

- Nominal slot volume 1000 mm³ = 1.0e-6 m³ → particle delta lower/upper with slot = **1000** (matches `dp³`).
- Lower particles inside nominal slot AABB: **0** (slot opens upward).
- Key sizes at lattice centers: `r_max≈19.96 mm`, body height span `19 mm` (half-spacing inset from ±20 mm ends — expected for lattice).
- P1 volume acceptance suggestion ≤2%: **passed** at coarse spacing.

## Runs executed

| Run | Command | Log |
|-----|---------|-----|
| geometry_only + lower slot | `./test_3d_spdb_geometry --geometry_only` | `results/geometry_only_with_slot.log` |
| geometry_only no slot | `./test_3d_spdb_geometry --geometry_only --no_slot` | `results/geometry_only_no_slot.log` |
| elastic no slot (smoke) | `./test_3d_spdb_geometry --elastic_no_slot` | `results/elastic_no_slot.log` (~8.2 s wall) |

VTK: `results/output_slot/`, `results/output_no_slot/`, `results/output_elastic_no_slot/`.
Sections (dimensions annotated): `results/section_xz_with_slot.svg`, `results/section_yz_with_slot.svg`.

## Elastic baseline note

Short dynamic-relaxation-style smoke (`end_time=2e-5 s`, light gravity on upper, bottom holder, contact + damping).
**Not** a validated quasi-static 10 MPa hold, **not** analytical compression certification, **not** SYCL.

## Next-round gaps

1. Proper platen BC + nominal 10 MPa force control; force/penetration checks.
2. Analytical low-load axial compression comparison (no slot).
3. Temperature field / thermal expansion (still prescribed T(t), no E(T) yet).
4. Creep / viscoplastic constitutive (missing in library).
5. Slot model under hold; self-contact near closure.
6. CPU↔SYCL twin; finer spacings 0.5 / 0.25 mm.
7. Electro-thermal (P4) and config/LLM (P5) — out of scope until mechanics verified.
