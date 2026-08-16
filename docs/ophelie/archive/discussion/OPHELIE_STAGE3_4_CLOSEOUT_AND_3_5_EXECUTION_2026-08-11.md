# Stage 3.4 closeout + Stage 3.5 execution notes (2026-08-11)

Binding plan: `docs/ophelie/OPHELIE_STAGE3_4_CLOSEOUT_AND_FRENCH_NATURAL_NEXT_PLAN_2026-08-11.md`.

## Accepted GPT decisions (no blockers)

1. Dual-layer φ gate: preferred `2.0e-4`, engineering hard `2.5e-4` (warn in between).
2. Production Natural σ: journal Table-1 anchor 16 @1473 K + thesis III.12 shape, clip `[1e-6, 30]` via `--sigma-law=french-natural`.
3. Raw thesis III.12 remains sensitivity only.
4. Short Stage 3.5 = 50 kW target-power + A_glass sanity, then leave EM for thermal/fluid.

## Code landed

| Item | Where |
|------|--------|
| `σ_nat(T)=16·10^(3179.8·(1/1473−1/T))` clipped | `electromagnetic_ophelie_french_material_laws.h` |
| `--sigma-law=french-natural` (aliases: `cep2008`, `natural`) | joule-to-heat test |
| Dual-layer φ acceptance + `phi_gate_warn` | joule-to-heat test |
| `--phi-pcg-tol=` / `--phi-pcg-max-iter=` | `electromagnetic_ophelie_cli.h` |
| `--target-power[=W]` calibrate I then re-Picard | joule-to-heat test (`|P−P_t|/P_t < 1e-3`) |

## Re-run recipes (normal terminal, cwd=`build/`)

```bash
RELAX_RELOAD=~/SPHinXsysSYCL/build/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_natural_glass_relax/bin/reload
BIN=./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_complex_joule_to_heat_one_way/bin/test_3d_ophelie_french_complex_joule_to_heat_one_way
GEO="--reload-dir=$RELAX_RELOAD --dp=0.015 --glass-radius=0.25 --glass-height=0.185 --frequency=282000 --sigma=16 --coil-radius=0.285 --coil-num-loops=7 --coil-segments-per-loop=64 --coil-z-min=-0.0225 --coil-z-max=0.2075 --ophelie-edge-flux-complex=1"

# 3.4b engineering closeout (expect phi_gate_warn=1 if residual ~2.1e-4)
$BIN $GEO --self-induction --self-induction-max-iter=25 --self-induction-relax=0.3 \
  --self-induction-tol=1e-3 --self-induction-phi-tol=2e-4 \
  --sigma-t --sigma-law=thesis-iii12 --sigma-t-max-iter=8 --sigma-t-relax=0.3 --sigma-t-tol=0.05 \
  --sigma-t-seed-delta-k=100 --thermal-dt=0.1 --thermal-steps=1

# Optional tighter-PCG diagnostic (does not block)
$BIN $GEO --self-induction --self-induction-max-iter=25 --self-induction-relax=0.3 \
  --self-induction-tol=1e-3 --self-induction-phi-tol=2e-4 \
  --phi-pcg-tol=5e-5 --phi-pcg-max-iter=20000 \
  --sigma-t --sigma-law=thesis-iii12 --sigma-t-max-iter=8 --sigma-t-relax=0.3 --sigma-t-tol=0.05 \
  --sigma-t-seed-delta-k=100 --thermal-dt=0.1 --thermal-steps=1

# Stage 3.5 — 50 kW + A_glass (fixed σ=16 or french-natural seed)
$BIN $GEO --self-induction --self-induction-max-iter=25 --self-induction-relax=0.3 \
  --self-induction-tol=1e-3 --self-induction-phi-tol=2e-4 \
  --target-power=50000 --thermal-dt=0.1 --thermal-steps=1
```

Expect Stage 3.5: `power_ok=1`, `P_joule_W≈50000`, `phi_ok=1` (possibly `phi_gate_warn=1`).

## After 3.5

Do **not** stay on EM residuals. Next: Stage 4.0 production thermal BC (Robin/radiation), then frozen-Q natural convection.
