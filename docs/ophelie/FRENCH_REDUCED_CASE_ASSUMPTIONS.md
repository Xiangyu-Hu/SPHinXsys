# French Reduced Case — Literature vs Assumptions

This case is a **French-paper-inspired reduced geometry**, not an exact reconstruction of the CEA Marcoule experimental CAD.

Primary reference: Jacoutot et al., *Numerical modelling of natural convection in molten glass heated by induction*, Chem. Eng. Process. 47 (2008) 449–455.

## Literature-based (from paper)

| Item | Value | Source |
|------|-------|--------|
| Cold crucible diameter | 650 mm → `glass_radius = 0.325 m` | Introduction |
| Frequency | 300 kHz | Introduction |
| Generator rating | 400 kW | Introduction |
| Current | sinusoidal | Introduction |
| EM solver concept | OPHELIE integral method, Biot–Savart (Eq. 11) | §2.2 |
| Glass mesh | volume mesh in molten bath | §2.2 |
| σ at 1473 K | 16 S/m | Table 1 |
| ρ at 1473 K | 2750 kg/m³ | Table 1 |
| cp at 1473 K | 1150 J/(kg·K) | Table 1 |
| λ (k) at 1473 K | 4 W/(m·K) | Table 1 |
| Example total Joule power | 50 kW (natural convection case) | §3 |
| Thermo/flow domain | glass bath only | §2.1 |

Stirring paper (separate): similar 650 mm / 300 kHz / 400 kW; example **P_tot = 60 kW**, stirrer 10 rev/min.

## NOT in paper (reduced-case assumptions)

These are **CLI parameters with reasonable defaults** — do not claim they come from the French CAD:

| Parameter | Default | Notes |
|-----------|---------|-------|
| `glass_height` | 0.50 m | Fig. 1 labels h but gives no numeric value |
| `coil_radius` | 0.40 m | outside glass, parametric |
| `coil_z_min/max` | 0.05 / 0.45 m | auto from glass height + inset |
| `coil_num_loops` | 8 | parametric |
| `coil_segments_per_loop` | 256 | quadrature resolution |
| `dp` | 0.02 m | SPH resolution |
| `crucible_wall_thickness` | 0.02 m | visual only |
| `ampere_turns` | 8 | absolute scale before power normalization |
| Thermal material (one-way) | reduced: ρ=2500, cp=1200, k=1, T₀=300 | regression default; see `--thermal-material=literature` for Table 1 @ 1473 K |

## What we model

```text
analytic GlassBody cylinder (GeometricShapeCylinder)
    + code-generated multiloop circular line source (no coil SPH particles)
        → Biot–Savart A/B on glass
            → optional PhiImag
                → E / J / JouleHeat
                    → optional power normalization to 50 kW (or 60 kW)
```

Use **`--literature-mode`** for Jacoutot OPHELIE form without post-hoc field scaling (see [FRENCH_LITERATURE_MODE.md](FRENCH_LITERATURE_MODE.md)).

## What we do NOT model

- Segmented cold crucible EM (surface skin-depth on walls)
- STL / particle relaxation geometry
- Water cooling, σ(T), thermal feedback
- Stirrer, skull layer thickness, full 3D crucible periodic mesh
- Full PDE A–φ solver, self-induction

## Optional visualization bodies

- `--coil-visual`: thin annular cylinder at coil radius (no EM)
- `--crucible-visual`: hollow cylindrical wall shell (no EM, future thermal BC)

## When real French geometry arrives

Replace only:

1. `GlassBody` shape (STL or refined analytic)
2. `CoilSource` loop distribution (spiral / real ampere-turns layout)
3. Visual wall / coil geometry

EM pipeline (`Biot → phi → E/J/Q → scaling → VTP`) stays unchanged.
