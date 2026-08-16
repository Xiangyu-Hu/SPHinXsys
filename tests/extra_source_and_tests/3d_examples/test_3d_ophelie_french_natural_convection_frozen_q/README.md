# test_3d_ophelie_french_natural_convection_frozen_q

Stage 4.1: frozen `Q_50kW` → WCSPH + Boussinesq + French Natural thermal BC (no EM update in the flow loop).

See `docs/ophelie/archive/discussion/OPHELIE_STAGE4_1_FROZEN_Q_NATURAL_CONVECTION_2026-08-11.md`.

## Run (cwd=`build/`)

```bash
cmake --build . --target test_3d_ophelie_french_natural_convection_frozen_q -j$(nproc)

RELAX_RELOAD=~/SPHinXsysSYCL/build/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_natural_glass_relax/bin/reload

./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_natural_convection_frozen_q/bin/test_3d_ophelie_french_natural_convection_frozen_q \
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
