# French Reduced — Literature Mode (Jacoutot 2008 OPHELIE)

Aligns the **electromagnetic model** with Jacoutot et al. (2008) as far as geometry remains reduced.

## Literature vs demo

| | **Literature mode** | **Demo mode** (default) |
|--|---------------------|-------------------------|
| CLI | `--literature-mode` | (no flag) |
| φ correction | **ON** (overrides `--no-phi`) | default ON, can `--no-phi` |
| Power | **No post-hoc field scaling** | `P_scaled` via `field_scale` |
| Coil current | Optional **calibrate** to `P_raw ≈ target` | `I=1 A/loop` then scale fields |
| φ solver | Defaults to **GMRES** if Jacobi | PCG default in params |
| Pass criteria | `literature_passed` (divJ + φ + P_raw) | finite fields + scaled power |

Geometry (cylinder, filament coil) is still **reduced-case** — see [FRENCH_REDUCED_CASE_ASSUMPTIONS.md](FRENCH_REDUCED_CASE_ASSUMPTIONS.md).

## Recommended command

```bash
cd ~/SPHinXsysSYCL/build
BIN=./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_reduced/bin/test_3d_ophelie_french_reduced

# 1) relax once
$BIN --glass-radius=0.325 --glass-height=0.50 --dp=0.02 \
  --relax=1 --relax-steps=400 --state_recording=0

# 2) literature EM (φ + physical current, no field scaling)
$BIN --glass-radius=0.325 --glass-height=0.50 --dp=0.02 \
  --reload=1 --literature-mode --state_recording=1
```

## Physics (matches paper §2.2 form)

```math
\hat{\mathbf{j}} = -\sigma\left(\nabla\hat{V} + i\omega\hat{\mathbf{A}}_{src}\right),\quad
Q = \frac{|\hat{\mathbf{j}}|^2}{2\sigma}
```

Parameters: **300 kHz**, **σ = 16 S/m**, **P ≈ 50 kW** (via coil-current calibration, not `field_scale`).

## Acceptance (`literature_passed=1`)

1. φ solver residual within tolerance  
2. `P_raw` within **5%** of `--target-power` after coil-current calibration  
3. `divJ_L2_L0 / divJ_L2_phi ≥ --divj-l2-red-min` (default **1.25**; literature target **≥ 1.5**)  
4. A/B/E/J/Q finite and positive  

Console prints:

```text
[ophelie] literature_acceptance: ... literature_passed=0|1
test_3d_ophelie_french_reduced ... mode=literature ... passed=0|1
```

## Extra CLI

```text
--literature-mode              OPHELIE-aligned profile
--no-literature-calibrate      Keep I_per_loop; do not scale current to target P
--divj-l2-red-min=1.25         Gate for divJ_L2 improvement
--phi-gauge-penalty=1.0        φ null-space penalty (sweep 0.1, 1, 10)
--phi-solver=GMRES|PCG
--target-power=50000
--no-power-scaling             (set automatically by --literature-mode)
```

## Phi operator defaults (closed)

```text
LHS = DivSigmaGrad
RHS = DivSigmaA
LegacyFlux = deprecated diagnostic only (≈ −DivSigmaA)
eq_res_vol gate = 0.65 on reload Biot A (discrete floor ~0.57)
```

## Coil source model

EM uses **code-generated multiloop filament centerline** (`I·dl` quadrature), **not** coil volume/surface particles.  
`CoilVisualBody` is visualization only.

## Next stage (route 2)

See [`OPHELIE_NEXT_STAGE_DEVICE_GMRES_AIND_PLAN.md`](../../OPHELIE_NEXT_STAGE_DEVICE_GMRES_AIND_PLAN.md):

1. **One-way A_ind diagnostic** — `test_3d_ophelie_french_aind_diagnostic`  
2. Picard `--self-induction` (experimental, separate from `literature_passed`)  
3. device-resident GMRES (parallel engineering track)

`literature_passed=1` means **reduced internal acceptance** (source + φ + divJ), **not** full OPHELIE with self-induction.

Paste: full console line + `mode=literature` + `literature_acceptance` + VTP screenshot of `JouleHeat` / `DivJImag`.
