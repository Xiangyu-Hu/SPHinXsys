# OPHELIE Complex Edge-Flux — Run Data and Fix Record

> Companion: `OPHELIE_PHI_DEVELOPMENT_LOG.md` §4.17–§4.18 
> Updated: 2026-06-02

## 1. Executable paths (under build directory)

```bash
BIN=./tests/extra_source_and_tests/3d_examples

$BIN/test_3d_ophelie_french_aind_diagnostic/bin/test_3d_ophelie_french_aind_diagnostic
$BIN/test_3d_ophelie_edge_flux_scaling/bin/test_3d_ophelie_edge_flux_scaling
$BIN/test_3d_ophelie_french_self_induction_picard/bin/test_3d_ophelie_french_self_induction_picard
$BIN/test_3d_ophelie_french_reduced/bin/test_3d_ophelie_french_reduced
```

**Do not use** abbreviated `tests/.../` paths.

---

## 2. First run data (reload, n=20500, GPU)

### 2.1 A_ind one-way + feedback — **passed**

Command:
```bash
test_3d_ophelie_french_aind_diagnostic --reload=1 \
 --ophelie-current-form=edge-flux --ophelie-edge-flux-complex=1
```

| Metric | coil-only | feedback |
|------|-----------|----------|
| P_complex (W) | 7.27882 | **8.64382** (+18.8%) |
| max_J_imag | 64.07 | 64.07 |
| max_J_real | 0 | **21.86** |
| A_ind_imag/A_src_real | — | 0.433 |
| phi_eq_res_vol | 3.22e-4 | ~same order |
| passed | 1 | 1 |

**Conclusion**: Stage 3.5 one-way feedback works; `K[J_imag]→A_ind_imag` activates real chain.

**Known**: Before feedback real level0≈0 → `edge_res_red_real=-1` (n/a; numerator guard added).

### 2.2 I² scaling — **ratios passed, overall fail (test fixed)**

| Suite | P(0.5I)/P(1I) | P(2I)/P(1I) | target_P gate |
|------|---------------|-------------|---------------|
| imag_only | 0.25 ✓ | 4.0 ✓ | fail (P_ref≈1.7e-9 W) |
| complex | 0.25 ✓ | 4.0 ✓ | fail |

**Conclusion**: Stage 3.4 I² law holds on MMS box; absolute `target_P` calibration applies only to French reload, not MMS.

**Fix (2026-06-02)**: `test_3d_ophelie_edge_flux_scaling` gates only I² ratios; `target_P` is diagnostic output.

### 2.3 Picard — **false pass (CLI fixed)**

First run used `--ophelie-edge-flux-complex=1`, but the test **did not** call `filterOphelieTestCommandLine`, and **unconditionally** `applyFrenchLiteratureMode`:

- Actual path: imag-only `solvePhiImag` + grad_phi/EJQ
- Log: `self_induction` (not `self_induction_complex`)
- `phi_eq_res_vol≈0.29` (GMRES did not converge)
- `J_rel≈3e-6` false convergence (J barely moved)
- `passed=1` **invalid**

**Fix (2026-06-02)**:
- Picard test hooks `filterFrenchReducedCommandLine` + `filterOphelieTestCommandLine`
- Literature preset only when `--literature-mode`
- Output `current_form` / `edge_flux_complex` / `complex_picard_path` / `max_J_{real,imag}`
- Complex path adds soft gate `phi_eq_res_vol < 0.01`

---

## 3. Recommended rerun commands (after fix)

```bash
cd ~/SPHinXsysSYCL/build

# A_ind (accepted; regression)
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_aind_diagnostic/bin/test_3d_ophelie_french_aind_diagnostic \
 --reload=1 --ophelie-current-form=edge-flux --ophelie-edge-flux-complex=1

# I² scaling (expect edge_flux_scaling_passed=1)
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_edge_flux_scaling/bin/test_3d_ophelie_edge_flux_scaling

# Complex Picard (expect self_induction_complex + phi_gmres_real)
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_self_induction_picard/bin/test_3d_ophelie_french_self_induction_picard \
 --reload=1 --ophelie-current-form=edge-flux --ophelie-edge-flux-complex=1 --self-induction
```

---

## 4. Rerun results (after fix, 2026-06-02, build/reload, GPU)

**Working directory**: `cd ~/SPHinXsysSYCL/build` (reload reads `./reload/Reload.xml`)

### 4.1 I² scaling — **passed=1** ✅

```
edge_flux_scaling_passed=1 complex_scaling_passed=1
imag_only: P(0.5I)/P(1I)=0.25, P(2I)/P(1I)=4.0
complex: P(0.5I)/P(1I)=0.25, P(2I)/P(1I)=4.0
target_P diagnostic only (MMS P_ref≈1.7e-9 W)
```

### 4.2 Complex Picard — **passed=1** ✅

Command:
```bash
test_3d_ophelie_french_self_induction_picard --reload=1 \
 --ophelie-current-form=edge-flux --ophelie-edge-flux-complex=1 --self-induction
```

| Metric | Value |
|------|-----|
| `complex_picard_path` | 1 |
| Log path | `self_induction_complex` + `phi_pcg_real` |
| `self_induction_iters` | 7 |
| `final_J_rel` | 0.0385 (converged=1) |
| `phi_eq_res_vol` | 3.44e-4 (phi_ok=1) |
| `P_joule_W` | 6.39 |
| `max_J_real` / `max_J_imag` | 18.18 / 59.32 |
| `A_ind_over_A_coil` | 0.422 |
| `passed` | 1 |

### 4.3 French reduced complex @50kW — **passed=1** ✅

```bash
test_3d_ophelie_french_reduced --reload=1 \
 --ophelie-current-form=edge-flux --ophelie-edge-flux-complex=1 --target-power=50000
```

| Metric | Value |
|------|-----|
| `P_raw` | 7.279 W |
| `P_scaled` | **50000 W** |
| `phi_eq_res_vol` | 6.16e-4 |
| `phi_solver_passed` | 1 |
| `demo_passed` | 1 |
| real q_antisym | 0 (A_imag≈0, expected) |
| imag q_antisym_rel_l2 | 7.18e-8 |

### 4.4 A_ind one-way regression — **passed=1** ✅

| Metric | coil-only | feedback |
|------|-----------|----------|
| P_complex (W) | 7.279 | 8.644 (+18.8%) |
| `edge_res_red_imag` | — | 0.9998 |
| `edge_res_red_real` | — | **-1** (n/a) |
| `max_J_real` | 0 | 21.86 |

---

## 5. Picard joint convergence + Stage 3.7 power test (2026-06-02)

### 5.1 Picard J_rel + phi_eq_res joint gate ✅

**Implementation**:
- `ophelieSelfInductionPicardConverged`: outer loop break requires `J_rel < j_tol` **and** `phi_eq_res_vol < phi_tol`
- complex: `phi_tol=0.01` (`--self-induction-phi-tol=`); legacy: `phi_tol=phi_eq_res_vol_gate_=0.65`
- Log: `picard_converged=0|1` each outer iter

**Rerun**:
```bash
test_3d_ophelie_french_self_induction_picard --reload=1 \
 --ophelie-current-form=edge-flux --ophelie-edge-flux-complex=1 --self-induction
```

| Metric | Value |
|------|-----|
| `picard_converged` | 1 |
| `j_ok` / `phi_ok` | 1 / 1 |
| `self_induction_iters` | 7 |
| `final_J_rel` | 0.0385 |
| `phi_eq_res_vol` | 3.44e-4 |
| `passed` | 1 |

### 5.2 Complex uniform-field power (Plan §8.2) ✅

```bash
test_3d_ophelie_edge_flux_power_uniform_field
```

| Case | P_recon/P_exact | passed |
|------|-----------------|--------|
| potential_field (H v4) | 1.0 | 1 |
| induction_field (H v4) | 1.0 | 1 |
| complex_potential_real | 1.0 | 1 |
| complex_induction_imag_chain | 1.0 | 1 |
| complex_induction_real_chain | 1.0 | 1 |

---

## 6. Stage 3.7 status note

**§8.2 / Stage 3.7 completed in a previous round** (not “not started”):
- File: `test_3d_ophelie_edge_flux_power_uniform_field.cpp` includes 5 cases (H v4×2 + complex×3)
- Re-verified: `summary passed=1` (2026-06-02)

---

## 7. Picard sweep + thermal one-way (2026-06-02)

### 7.1 Picard sweep ✅

```bash
test_3d_ophelie_french_self_induction_picard_sweep --reload=1
```

| relax \ max_iter | 6 | 8 | 12 |
|------------------|---|---|-----|
| 0.15 | fail | fail | **ok** |
| 0.30 | fail | **ok** (ref) | **ok** |
| 0.45 | **ok** | **ok** | **ok** |

`sweep_passed=1` (reference + sufficient budget grid points all converge)

### 7.2 Complex Joule → heat one-way ✅

```bash
test_3d_ophelie_complex_joule_to_heat_one_way
# P_recon≈103490 W, max_rel_err=0, passed=1
```

Shared interface: `electromagnetic_ophelie_joule_to_heat_one_way.h`

---

## 8. Stage 4 — French reload thermal one-way end-to-end (2026-06-02) ✅

```bash
test_3d_ophelie_french_complex_joule_to_heat_one_way --reload=1 \
 --ophelie-current-form=edge-flux --ophelie-edge-flux-complex=1
```

| Metric | Value |
|------|-----|
| `P_joule_W` | 7.28 |
| `phi_eq_res_vol` | 3.2e-4 |
| `thermal_steps` / `dt` | 3 / 1 s |
| `vol_weighted_rel_err` | 0.124 |
| `energy_vs_power_rel_err` | 0.184 |
| `max_delta_T` | 9.2e-5 K |
| `passed` | 1 |

Device CK: `ApplyOphelieJouleHeatOneWayTemperatureStepCK` 
Shared pipeline: `runFrenchReducedEmThenJouleHeatOneWay`

**Known**: French grid integrated thermal energy ~18% below `P·Δt` (prototype gate 20%); uniform-field MMS passes strictly.

---

## 9. Todo

- [x] French reload complex thermal one-way end-to-end
- [ ] Tighten French energy closure (<5%)
- [ ] Multi-step thermal diffusion / crucible boundary BC (real diffusion, not frozen Q)
