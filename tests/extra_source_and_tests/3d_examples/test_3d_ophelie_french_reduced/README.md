# test_3d_ophelie_french_reduced

French-paper-inspired reduced cold-crucible induction case using **analytic SPHinXsys geometry** (cylinder level-set) with **particle relaxation and body fitting** (same SYCL-CK pipeline as TEAM7).

See [docs/ophelie/FRENCH_REDUCED_CASE_ASSUMPTIONS.md](../../../../docs/ophelie/FRENCH_REDUCED_CASE_ASSUMPTIONS.md) for literature vs assumed parameters.

## Model

- **GlassBody**: `GeometricShapeCylinder`, conductive glass (σ registered on particles), **relaxed to cylinder surface**
- **CoilSource**: multiloop circular **line quadrature** (current moment I·dl) on **code-generated centerline segments** — no coil conductor particles for EM
- **EM**: Biot–Savart → optional PhiImag → E/J/JouleHeat → optional 50 kW power scaling
- **Optional**: `--coil-visual`, `--crucible-visual` (ParaView only, no EM; also relaxed when present)

## Recommended workflow

```bash
cd ~/SPHinXsysSYCL/build
BIN=./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_reduced/bin/test_3d_ophelie_french_reduced

# 1) SYCL body-fitted relax → build/reload/Reload.xml
$BIN --glass-radius=0.325 --glass-height=0.50 --dp=0.02 \
  --relax=1 --state_recording=0 --relax-log-every=50

# 2) Literature OPHELIE (φ + physical I, no field scaling) — see docs/ophelie/FRENCH_LITERATURE_MODE.md
$BIN --glass-radius=0.325 --glass-height=0.50 --dp=0.02 \
  --reload=1 --literature-mode --state_recording=1

# 2b) Same with device GMRES (validated: test_3d_ophelie_phi_device_vector_ops --gmres-parity)
$BIN --glass-radius=0.325 --glass-height=0.50 --dp=0.02 \
  --reload=1 --literature-mode --phi-gmres-device-ops=1 --state_recording=1
# Optional full device Krylov: add --phi-gmres-device-krylov=1

# 3) Demo / thermal handoff (scaled Q to 50 kW)
$BIN --glass-radius=0.325 --glass-height=0.50 --dp=0.02 \
  --reload=1 --no-phi --state_recording=1
```

## A_ind one-way diagnostic (route 2)

After literature φ chain is accepted, measure **||A_ind||/||A_coil||** without feeding A_ind back:

```bash
cd ~/SPHinXsysSYCL/build
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_aind_diagnostic/bin/test_3d_ophelie_french_aind_diagnostic --reload=1

# Picard self-induction (experimental, not literature_passed)
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_self_induction_picard/bin/test_3d_ophelie_french_self_induction_picard --reload=1 --self-induction-max-iter=5

# Or on the main French reduced driver (multiloop + VTP output):
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_reduced/bin/test_3d_ophelie_french_reduced --reload=1 --literature-mode --self-induction --state_recording=0
```

Uses same `DivSigmaGrad` + `DivSigmaA` as literature-mode. Picard `--self-induction` remains **experimental** and is not part of `literature_passed`.

## Literature-aligned run (Jacoutot 2008 natural convection)

```bash
$BIN \
  --glass-radius=0.325 --glass-height=0.50 --dp=0.02 \
  --coil-radius=0.40 --coil-num-loops=8 --coil-segments-per-loop=256 \
  --frequency=300000 --sigma=16 --target-power=50000 \
  --reload=1 --no-phi --state_recording=1
```

## Relax / reload CLI

| Flag | Meaning |
|------|---------|
| *(default)* | Lattice → **relax** → EM |
| `--relax=1` | Relax only; write `Reload.xml` and exit |
| `--reload=1` | Load `GlassBody` from `reload/Reload.xml`; skip relax |
| `--reload-dir=PATH` | Folder containing `Reload.xml` |
| `--skip-relax` | Lattice only (fast debug, poor surface fit) |
| `--literature-mode` | φ ON, no field scaling, coil-current calibrate → `literature_passed` |
| `--no-literature-calibrate` | With literature-mode: keep `--coil-current` |
| `--divj-l2-red-min=1.25` | Literature acceptance gate on divJ_L2 |
| `--phi-gauge-penalty=1.0` | φ gauge penalty λ |
| `--phi-lhs-operator=div-sigma-grad` | Solver LHS: `div(σ∇φ)` (literature default) |
| `--phi-rhs-operator=div-sigma-a` | Solver RHS from A_src: `-ω·div(σA)` (literature-mode **always** forces this; flag only affects non-literature runs) |

RHS A/B comparison: `test_3d_ophelie_phi_biot_rhs_solvability` (`--reload=1` reads particles from `build/reload`); CSV → `./output/ophelie_phi_rhs_solvability.csv`.
| `--phi-gmres-max-outer-iter=120` | GMRES outer iterations |
| `--phi-gmres-restart=80` | GMRES restart dimension |
| `--phi-gmres-eq-res-tol=0.42` | Stop GMRES when `phi_eq_res_vol` below this (0=off) |
| `--phi-eq-res-gate=0.65` | Literature `phi_res` acceptance on `phi_eq_res_vol` |
| `--relax-steps=400` | Relaxation iterations (default 400) |
| `--relax-log-every=50` | Progress log interval |
| `--relax-vtp-every=100` | Relax snapshot VTP interval (`0` = off) |

`Reload.xml` 只应出现在 **build 树** 的 `.../bin/reload/`（由 `--relax=1` 写出，`--reload=1` 读取）。**不要在源 case 目录下放置 `reload/`、`bin/`、`output/`。** `--relax=1` 与 `--reload=1` 使用的几何参数（`--glass-radius`、`--dp` 等）必须一致。

## Geometry / EM CLI

```text
--glass-radius=0.325
--glass-height=0.50
--glass-center=0,0,0.25
--dp=0.02

--coil-radius=0.40
--coil-z-min=0.05
--coil-z-max=0.45
--coil-num-loops=8          (alias: --coil-turns=)
--coil-segments-per-loop=256
--ampere-turns=8
--coil-current=1.0

--frequency=300000
--sigma=16
--target-power=50000
--no-power-scaling
--no-phi
--phi-solver=PCG

--coil-visual
--crucible-visual
--crucible-wall-thickness=0.02
```

## VTP fields (GlassBody)

`Sigma`, `ASrcReal/Imag`, `BSrcReal/Imag`, `EImag`, `JImag`, `JouleHeat`, `DivJImag`, optional `PhiImag`, `GradPhiImag`.

Relax progress VTP: `build/output/GlassBody_ite_*.vtp` (when `--relax-vtp-every>0`).

## divJ diagnostics

Console prints both relative and L2 metrics:

```text
divJ_L0 / divJ_phi / divJ_red
divJ_L2_L0 / divJ_L2_phi / divJ_L2_red
```

Prefer **divJ_L2** when comparing phi correction on induction cases.
