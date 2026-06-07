# test_3d_aphi_ck_cold_crucible_demo

## Purpose

This case demonstrates the current A-phi electromagnetic post-processing pipeline for a cold-crucible-like geometry. **It is not a quantitative validation case.**

This case is a pipeline demonstration / smoke test, not quantitative physical validation.

## What is demonstrated

- analytic / impressed vector potential **A** around a crucible-like domain (toroidal coil profile in x-y)
- **B = curl A** diagnostic using the corrected `LinearCorrectionMatrix` + `LinearGradient` route
- frequency-domain **E = -i omega A - grad phi** with **phi = 0**
- induced **J = sigma E** in conductive melt / wall regions
- **Joule** heat source `q = 0.5 sigma |E|^2`
- **H = nu B**
- ParaView **VTP** output under `output/` (scalar magnitudes included)

## What is not demonstrated yet

- fully validated coil-current electromagnetic solve (GMRES with impressed RHS)
- real cold crucible slit geometry
- multibody Contact split (single-body region tagging for MVP-A)
- wall dummy / support-complete boundary treatment
- quantitative TEAM7 / cold-crucible benchmark accuracy
- fluid-thermal two-way coupling

## Build and run

```bash
cd build
ninja test_3d_aphi_ck_cold_crucible_demo
./tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_cold_crucible_demo/bin/test_3d_aphi_ck_cold_crucible_demo
```

Expect: `demonstration_only=1 passed=1` and VTP files under `output/`.

## Suggested ParaView visualization

1. Load VTP from `output/`.
2. Color by `MaterialRegionId` (0=air, 1=wall, 2=melt, 3=coil).
3. Color melt particles by `JouleHeatSource` or `JMagnitude`.
4. Color air particles by `BMagnitude` with opacity 0.2–0.4.
5. Glyph `CurrentDensityReal` or `JMagnitude` on melt if desired.

## Material regions (demo defaults)

| Region | sigma | nu |
|--------|-------|-----|
| air / coil | 1e-8 | 1 |
| melt | 1.0 | 1 |
| crucible wall | 0.1 | 1 |

`omega = 100` (nondimensional demo value).
