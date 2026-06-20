# OPHELIE φ Operator Consistency Audit

> Internal engineering document. GPT diagnostic stage completed; no further discussion round needed. This document records operator discrepancies and fix directions.

## 1. Issue identification (confirmed)

| Diagnostic | Result | Meaning |
|---|---|---|
| gradient_linear | passed | GradPhi sign correct |
| laplace_rhs_consistency | eq_res ~ 1.9–2.7 | Under MMS, L(φ) ≠ RHS |
| French eq_res | ~ 1e-4 | Real-case solver OK |
| reconstructed_divJ_red | ~ 0.74 (penalty-independent) | Postprocess divJ does not improve |
| lattice vs reload | both ~0.74–0.79 | Boundary is not the main cause |

**Conclusion:** The solver solves `PairwiseLaplace(φ) = RHS`, but the divJ diagnostic follows `GradPhi → E → J → DivKernel(J)` — two inconsistent discrete operator sets.

## 2. Three core kernel comparison

### 2.1 PairwiseLaplace (solver LHS)

File: `electromagnetic_ophelie_laplace.h`

```cpp
pair_weight = harmonicMean(σ_i, σ_j) * (-2 * r * dW_ijV_j / (r² + ε h²))
laplace_i += pair_weight * (φ_i - φ_j)
```

- Uses **distance regularization** `pair_weight_regularization` (default 0.01)
- σ uses **harmonic mean**

### 2.2 ComputeOphelieScalarPhiGradientCK (postprocess)

File: `electromagnetic_ophelie_phi.h`

```cpp
g_ij = -dW_ijV_j * e_ij          // no distance regularization
grad_phi += g_ij * (φ_i - φ_j)
```

### 2.3 PhiRhsFromASrc (solver RHS)

File: `electromagnetic_ophelie_phi.h`

```cpp
g_ij = -dW_ijV_j * e_ij          // same as GradPhi
rhs_i -= ω * harmonicMean(σ) * g_ij · (A_i - A_j)
```

### 2.4 ComputeOphelieVecdDivergenceCK (divJ diagnostic)

File: `electromagnetic_ophelie_diagnostics.h`

```cpp
grad_W_ij = -dW_ijV_j * e_ij     // same as GradPhi
div_i += (J_j - J_i) · grad_W_ij
```

## 3. Inconsistency matrix

| Operator A | Operator B | Weight form | Match? |
|---|---|---|---|
| Laplace(φ) | div(σ grad φ) | distance-weighted vs uncorrected | **No** |
| Laplace(φ) | RHS from A | different flux structure | **No** (MMS proved) |
| GradPhi | Div(J) | same uncorrected g_ij | Yes (mutually consistent) |
| GradPhi + E/J | Laplace | different discretization | **No** (root cause) |

Continuous theory: `div(σ ∇φ) = Laplace(φ)` (under appropriate BCs).  
Discrete implementation: Laplace uses `pairwiseNegativeLaplaceWeight`, Grad/Div use `pairwiseGradientWeightUncorrected` → **even if φ exactly satisfies Laplace=RHS, GradPhi-reconstructed div(σ grad φ) does not equal Laplace(φ)**.

## 4. New test (already run)

```bash
test_3d_ophelie_phi_grad_div_laplace_consistency
```

Results (n=216, dp=0.08, φ=cos MMS):

| Metric | Value | Interpretation |
|---|---|---|
| graddiv_vs_laplace_vol | **0.714** | div(σ∇φ) vs Laplace(φ) relative error ~71%, **operator mismatch** |
| graddiv_vs_laplace_linf | 0.517 | same as above |
| discrete_mms_EImag_norm | **4.4e-8** | E/J chain self-consistent (E→0 under same kernel) |

**Key conclusion:**
- The E/J postprocess chain itself has no sign bug;
- Root cause is **Laplace (distance-weighted) ≠ div(σ grad) (uncorrected)**;
- Even if the solver drives Laplace(φ)=RHS to eq_res~1e-4, GradPhi-reconstructed divJ will not improve.

## 5. Fix implementation (2026-06-02)

**Option A implemented:** Added `OpheliePhiLhsOperatorKind::DivSigmaGrad`

- **LHS:** `GradPhi → σ·GradPhi → Div` (same kernel as postprocess)
- **RHS:** `-ω · Div(σ·A_src)` (same Div kernel)
- **Default:** globally still `LegacyPairwise` (TEAM7/old demo unaffected)
- **literature-mode** automatically enables `DivSigmaGrad`

### Results after fix

| Test | Before fix | After fix |
|---|---|---|
| laplace_rhs (discrete) | eq_res ~2.7, failed | eq_res ~3.5e-7, **passed** |
| manufactured_A (discrete) | failed, divJ_red~0.6 | **passed**, divJ_L2_phi~0.01 |
| french literature divJ_red | ~0.74 | **~1.82** (divJ_continuity passed) |
| french phi_eq_res_vol | ~1e-4 (old operator) | ~0.55 (new operator, BC TBD) |

CLI: `--phi-lhs-operator=div-sigma-grad|legacy-pairwise`

### Shell diagnostic (French reload, r_int=0.85·R)

```
phi_eq_res_shell: interior≈0.56  boundary≈0.54  interior_vol_frac≈0.83
```

Boundary and interior residuals are the same order of magnitude → **eq_res ~0.55 is not simply missing boundary support**; more likely GMRES on 16k particles has not fully converged the div-consistent system.

Post-Jacobi with legacy diagonal + DivSigmaGrad LHS still diverges; disabled (`phi_post_jacobi_refinement_iterations=0`).

### GMRES convergence (French reload, n=16066)

- **max(legacy_diag, div_diag)** preconditioner + 350 outer: eq_res_vol stalls at **0.572** after starting from 0.55
- `rel_res_l2` and `eq_res_vol` stall together → **not insufficient iterations**, but ~57% residual floor from discrete equation on Biot A field
- penalty=0/1e-4/1/1e-2 has almost no effect on eq_res; divJ_red still ~1.75
- MMS with same operator: eq_res ~1e-7 → operator implementation correct; floor from **real-case solvability/discrete compatibility** (not purely boundary)

CLI: `--phi-gmres-max-outer-iter` `--phi-gmres-eq-res-tol` `--phi-eq-res-gate`

### TODO

- Study whether Biot `A_src` satisfies discrete solvability of `div(σ∇φ)=-ω div(σA)` (or needs RHS correction term)
- Surface flux / weak-form boundary

## 6. Client demo scope (unchanged)

```bash
test_3d_ophelie_french_reduced --reload=1 --no-phi --target-power=50000 --state_recording=1
```

Do not externally claim literature-mode passed / divJ verified.
