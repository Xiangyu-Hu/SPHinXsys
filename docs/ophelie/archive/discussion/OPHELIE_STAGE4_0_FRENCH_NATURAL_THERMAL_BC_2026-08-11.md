# Stage 4.0 — French Natural production thermal BC (2026-08-11)

Binding plan: `docs/ophelie/OPHELIE_STAGE3_4_CLOSEOUT_AND_FRENCH_NATURAL_NEXT_PLAN_2026-08-11.md` §§6–7.

## What landed

| Item | Detail |
|------|--------|
| Dirichlet kept as regression | `--thermal-bc=dirichlet-regression` (default with `--thermal-diffusion`) |
| Production BC | `--thermal-bc=french-natural` |
| Side wall | Robin `q = h_side (T − T_c)`, default `h_side=300` (within journal 100–400 range) |
| Bottom | Robin `q = h_bottom (T − T_c)`, default `h_bottom=35` |
| Free top | `q = h_free (T − T_a) + ε σ_SB (T⁴ − T_ar⁴)`, defaults `h_free=20`, `ε=0.8` |
| Diagnostics | `wall_loss_side_W`, `wall_loss_bottom_W`, `free_conv_loss_W`, `free_rad_loss_W`, `total_heat_loss_W` |

Surface fluxes are applied on a shell of thickness `φ_boundary_distance_factor * dp` via volumetric sink `q/δ`.

## Smoke run (cwd=`build/`, normal terminal)

```bash
RELAX_RELOAD=~/SPHinXsysSYCL/build/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_natural_glass_relax/bin/reload
BIN=./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_complex_joule_to_heat_one_way/bin/test_3d_ophelie_french_complex_joule_to_heat_one_way

cmake --build . --target test_3d_ophelie_french_complex_joule_to_heat_one_way -j$(nproc)

$BIN \
  --reload-dir="$RELAX_RELOAD" \
  --dp=0.015 --glass-radius=0.25 --glass-height=0.185 \
  --frequency=282000 --sigma=16 \
  --coil-radius=0.285 --coil-num-loops=7 --coil-segments-per-loop=64 \
  --coil-z-min=-0.0225 --coil-z-max=0.2075 \
  --ophelie-edge-flux-complex=1 \
  --self-induction --self-induction-max-iter=25 --self-induction-relax=0.3 \
  --self-induction-tol=1e-3 --self-induction-phi-tol=2e-4 \
  --target-power=50000 \
  --thermal-bc=french-natural --use-literature-thermal=1 \
  --thermal-dt=0.1 --thermal-steps=5
```

Expect: `passed=1`, `thermal_bc=french_natural_production`, all face loss channels `>0`, `P_joule_W≈50000`.

## Next

Stage 4.1: freeze `Q_50kW(x)` and open Boussinesq natural convection (no EM update yet).
