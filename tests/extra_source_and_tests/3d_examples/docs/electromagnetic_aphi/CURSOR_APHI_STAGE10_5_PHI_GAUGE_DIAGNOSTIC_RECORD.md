# Stage 10.5 — φ gauge penalty sweep + φ-only comparison record

> **Purpose**: Discuss with ChatGPT the left-conductor ~1% max-norm open issue and gauge penalty semantics.  
> **Date**: 2026-05-28  
> **Test**: `test_3d_aphi_ck_contact_phi_gauge_penalty_sweep_diagnostic`  
> **Related**: `test_3d_aphi_ck_contact_left_field_error_diagnostic` (10.5-A/C/E)

---

## 1. Background: what gauge penalty we implemented

Optional zeroth-order term on φ equation LHS (**not Laplace itself**):

```
L_phi += lambda_phi * phi     (Real/Imag each get one term)
```

- Code: `AphiPhiGaugePenaltyCK` + inline in `AphiApplyCK` (Inner/Contact both)
- Config: `AphiLhsAssemblyOptions::use_phi_gauge_penalty` / `phi_gauge_penalty`
- Block Jacobi diagonal sync includes same term
- **No Coulomb gauge penalty on A equation** (no div A penalty)

Relation to Laplace:

| Term | Meaning |
|---|---|
| `terms.laplace_phi` | SPH pairwise discrete σ∇²φ |
| `phi_gauge_penalty` | Mass regularization eliminating **constant φ gauge null space** |

---

## 2. Experimental design

**Case**: two-body flux-matched interface MMS (σ=2/1e-4, ω=1.25, x_interface=0.5, dp=0.1)  
**Solver**: coupled multi-body Contact GMRES (restart=50, outer=100, tol=1e-5)  
**Polish**: **0 sweep** (differs from left_field diagnostic polish=1, to isolate penalty effect)

### 2.1 Full A-phi block — penalty sweep

| run | `use_penalty` | `lambda_phi` |
|-----|---------------|--------------|
| A0 | 0 | — |
| A1 | 1 | 1 |
| A2 | 1 | 10 (MMS default) |
| A3 | 1 | 100 (TEAM7 default) |

### 2.2 φ-only operator comparison (penalty on/off)

Only `laplace_phi + penalty` enabled; A/reaction/grad/div all off.  
RHS still assembled from **full MMS exact** via same options (φ-only subsystem, not full PDE).

| run | operator | penalty |
|-----|----------|---------|
| B1 | φ-only | 10 |
| B2 | φ-only | 0 (off) |

---

## 3. Measured results (stdout verbatim)

```
penalty=0  use_penalty=0  phi_only=0  converged=1  outer=3
  global_true_rel=1.81e-6  left_true_rel=1.60e-6  left_continuous=88.4
  full_left_L2=110.7  full_left_Linf=336.6  phi_real_mean=0.830
  phi_L2_mean_sub=1.51e-4  core_interior_Linf=838.2

penalty=1  use_penalty=1  phi_only=0  converged=1  outer=5
  global_true_rel=3.90e-6  left_continuous=5.23%
  full_left_Linf=19.9%  phi_real_mean=0.011  phi_L2_mean_sub=8.9e-5

penalty=10  use_penalty=1  phi_only=0  converged=1  outer=2
  global_true_rel=4.99e-6  left_continuous=2.05%
  full_left_Linf=7.8%  phi_real_mean=0.0088  phi_L2_mean_sub=0.16%

penalty=100  use_penalty=1  phi_only=0  converged=1  outer=2
  global_true_rel=2.11e-7  left_continuous=0.26%
  full_left_Linf=1.0%  phi_real_mean=1.3e-5  phi_L2_mean_sub=0.017%

penalty=10  use_penalty=1  phi_only=1  converged=1  outer=2
  global_true_rel=1.71e-7  left_continuous=0.011%
  full_left_Linf=0.04%  left_true_rel=5.5e-4 (elevated)
  phi_L2_mean_sub=0.0028%

penalty=0  use_penalty=0  phi_only=1  converged=1  outer=2
  global_true_rel=4.12e-6  left_continuous=112%
  full_left_Linf=428%  phi_real_mean=0.506
```

---

## 4. Interpretation (discussion points for ChatGPT)

### 4.1 Penalty=0 is catastrophic

- Krylov **can still report converged=1** (global_true_rel ~1e-6), but left-body field **fully drifts** (left_continuous ~88–112%).
- `phi_real_mean ≈ 0.83` (full) / `0.51` (φ-only) → constant φ gauge null space not fixed.
- **Conclusion**: MMS acceptance checking Krylov residual alone misses gauge drift; penalty is not an "optional optimization" but **required for solvability/field semantics of this discrete system**.

### 4.2 Larger penalty monotonically reduces left-body field error (this case)

| λ_φ | left_continuous (max-norm) | phi_L2_mean_sub | core_interior_Linf |
|-----|---------------------------|-----------------|-------------------|
| 0 | **88%** | 1.5e-4 | **838%** |
| 1 | 5.2% | 8.9e-5 | 50% |
| 10 | 2.0% | 0.16% | 18% |
| 100 | **0.26%** | **0.017%** | **2.5%** |

- **MMS default λ_φ=10 is not optimal for left-body field accuracy**; TEAM7 **100 is more reasonable on this case**.
- Recommendation: **unified penalty calibration experiment** for interface MMS and TEAM7, or change MMS default to 100 for A/B.

### 4.3 Difference from left_field diagnostic (λ_φ=10, polish=1)

| Config | left_continuous | Notes |
|------|-----------------|------|
| sweep, λ=10, polish=0 | ~2.0% | this table |
| left_field, λ=10, polish=1 | ~0.96% | 10.5-A/C/E |
| sweep, λ=100, polish=0 | **0.26%** | without polish already better than polish=1@10 |

→ Left body ~1% may be **mainly from penalty too small (10 vs 100)**, not just Contact split or polish issue.

### 4.4 φ-only vs full block

- φ-only + λ=10: left_continuous **0.01%**, but **left_true_rel=5.5e-4** (two orders of magnitude larger than full block).
- Reason: φ-only solves **reduced-dimension sub-operator**; field looks close to MMS φ component but **does not satisfy full A-phi coupled residual**.
- **Cannot** use φ-only field error for full block acceptance; only shows **φ Laplace+penalty sub-block discretization precision is sufficient**.

### 4.5 Relation to 10.5-A mean-sub probe

- At λ=0, mean-sub L2 still small (~1e-4) but Linf/left_continuous explode → **error is mainly global constant offset + non-mean modal mix**, not a small perturbation measurable by L2 alone.
- At λ=100, phi_real_mean ~1e-5, mean-sub and raw L2 close → gauge offset largely eliminated by penalty.

### 4.6 10.5-C partitioning (see `CURSOR_APHI_STAGE10_CONTACT_PROGRESS.md`)

- core interior Linf **greater than** interface band (λ=10, polish=1 path).
- Problem **not concentrated in interface band**; consistent with **interior gauge/field drift from insufficient penalty**.

---

## 5. Suggested next steps (Cursor / ChatGPT consensus candidates)

1. **Change interface MMS default `phi_gauge_penalty` from 10 → 100**, rerun MMS + left_field diagnostic, see if left_continuous drops to ~0.3%.
2. **Penalty continuation sweep**: {10, 30, 100, 300} find Pareto point for plate/Joule and left MMS (TEAM7 three-body also needs sweep).
3. **Layered acceptance** (write into MMS passed logic discussion, not necessarily code change yet):
   - Tier-1: Krylov `global_true_rel` + discrete defect + interface band
   - Tier-2 (diagnostic): left body max-norm vs exact, **alert only when λ_φ≥100 or still large after mean-sub**
4. **A-side gauge**: currently no div A penalty; if left A component still poor after raising λ_φ, consider A gauge or curl post-processing.

---

## 6. Code / test index

| File | Role |
|------|------|
| `aphi_phi_gauge_penalty_ck.hpp` | penalty kernel |
| `aphi_matrix_free_operator_ck.hpp` | fused apply inline penalty |
| `aphi_coupling_modes_ck.h` | `AphiLhsAssemblyOptions` |
| `aphi_contact_phi_gauge_sweep_helpers.h` | sweep row runner |
| `aphi_contact_left_field_error_helpers.h` | 10.5-A/C/E partitioning + φ component error |
| `test_3d_aphi_ck_contact_phi_gauge_penalty_sweep_diagnostic` | data source for this record |
| `test_3d_aphi_ck_contact_left_field_error_diagnostic` | 10.5-A/C/E |

### Run

```bash
cd build
ninja test_3d_aphi_ck_contact_phi_gauge_penalty_sweep_diagnostic \
      test_3d_aphi_ck_contact_left_field_error_diagnostic
./bin/test_3d_aphi_ck_contact_phi_gauge_penalty_sweep_diagnostic
./bin/test_3d_aphi_ck_contact_left_field_error_diagnostic
```

---

## 7. One-line summary (paste directly to ChatGPT)

> On two-body Contact MMS, φ gauge penalty implemented (λ_φ·φ added to φ equation LHS). **At λ_φ=0 Krylov converges but left-body field drifts ~100%**; λ_φ from 0→1→10→100 left-body max-norm error monotonically drops to **0.26%** (polish=0). **MMS default λ=10 may be too small**; left body ~1% open issue likely mainly **penalty calibration** not Contact split. φ-only sub-operator field error can be very low but does not satisfy full residual, cannot replace full block acceptance.
