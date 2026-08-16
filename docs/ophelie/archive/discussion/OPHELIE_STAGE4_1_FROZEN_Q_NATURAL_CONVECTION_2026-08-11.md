# test_3d_ophelie_french_natural_convection_frozen_q

Stage 4.1 smoke: freeze `Q_50kW` from EM + A_glass Picard, then WCSPH natural convection with Boussinesq buoyancy and French Natural thermal BC. **No EM update during flow.**

## Run (cwd = `build/`, normal terminal)

```bash
cmake --build . --target test_3d_ophelie_french_natural_convection_frozen_q -j$(nproc)

RELAX_RELOAD=~/SPHinXsysSYCL/build/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_natural_glass_relax/bin/reload
BIN=./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_natural_convection_frozen_q/bin/test_3d_ophelie_french_natural_convection_frozen_q

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
  --end-time=0.02
```

## Gates

| Check | Criterion |
|-------|-----------|
| EM | `|P−50kW|/50kW < 1%` |
| Flow | `U_max`, `T_mean` finite |
| Thermal BC | side/bottom/free losses `> 0` |
| Buoyancy | `U_max > 0` |

## Notes

- Side+bottom wall particles = no-slip; open top ≈ free-slip (no top viscous wall).
- First slice: frozen-Q Lagrangian heating + Natural Robin/radiation; conduction Laplace deferred.
- Not a steady Ra/Pr validation — smoke that thermal-fluid chain runs.
