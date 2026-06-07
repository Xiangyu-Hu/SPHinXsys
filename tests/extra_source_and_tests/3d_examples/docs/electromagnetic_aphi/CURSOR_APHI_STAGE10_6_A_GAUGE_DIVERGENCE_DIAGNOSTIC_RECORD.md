# Stage 10.6 — A gauge / div(A) post-solve diagnostic record

> **Date**: 2026-05-28  
> **Basis**: ChatGPT [`stage10_contact_phi_penalty_and_divA_plan_for_cursor.md`](../../../../../../docs/electromagnetic_aphi/stage10_contact_phi_penalty_and_divA_plan_for_cursor.md)  
> **Related**: [`CURSOR_APHI_STAGE10_5_PHI_GAUGE_DIAGNOSTIC_RECORD.md`](CURSOR_APHI_STAGE10_5_PHI_GAUGE_DIAGNOSTIC_RECORD.md)

---

## 1. Clear conclusions (implementation status)

| Mechanism | Status | Notes |
|------|------|------|
| `phi_gauge_penalty` (`L_φ += λ_φ φ`) | ✅ implemented | **Numerical regularization**, locks φ constant gauge; **not** div A = 0 |
| `div(σA)` coupling (φ equation) | ✅ production path | Current continuity, **≠** Coulomb gauge |
| `div A` post-solve diagnostic | ✅ Stage 10.6 | `LinearGradient` + `AphiVectorGradientDivergenceCK` → `DivAReal/Imag` |
| `A ← A - ∇χ` projection | ❌ not implemented | deferred |
| `λ_A ∇(∇·A)` A-divergence penalty | ❌ not implemented | large divA → Stage 10.7 |

**External messaging**: phi penalty **cannot** substitute for div A = 0; current solution **does not** explicitly satisfy Coulomb gauge.

---

## 2. Code

| File | Role |
|------|------|
| `aphi_div_a_diagnostic_helpers.h` | post-solve pipeline + L2/Linf/relative reduction |
| `aphi_contact_phi_gauge_sweep_helpers.h` | sweep row includes divA + three-body sweep |
| `test_3d_aphi_ck_contact_phi_gauge_lambda100_ab_diagnostic` | λ=10 vs 100 A/B |
| `test_3d_aphi_ck_contact_team7_phi_gauge_sweep_diagnostic` | three-body λ∈{10,30,100,300} |
| `test_3d_aphi_ck_contact_phi_gauge_penalty_sweep_diagnostic` | extended divA output |

### divA algorithm (diagnostic-only)

1. Per-body `Inner<>`: `LinearCorrectionMatrix` + `LinearGradient<Vecd>` → `ARealGradient` (Matd)
2. `AphiVectorGradientDivergenceCK`: `div A = trace(∇A)` → register **`DivAReal` / `DivAImag`**
3. Host reduction: `div_A_L2`, `div_A_Linf`, `grad_A_L2`, **`div_A_relative = div_A_L2 / grad_A_L2`**

**Limitations**: per-body Inner gradient (consistent with debug LHS); interface cross-body gradient not in Contact. **Production GMRES does not read DivA fields.**

### Diagnostic levels (initial thresholds)

```text
div_A_relative < 1e-2   → good
[1e-2, 1e-1)            → warn
>= 1e-1                 → high_risk
```

---

## 3. λ=100 A/B + divA (two-body MMS, polish=0)

```
penalty=10:  left_continuous≈2.05%  div_A_relative≈1.113  div_A_level=high_risk
penalty=100: left_continuous≈0.26%  div_A_relative≈1.019  div_A_level=high_risk
             left_continuous_improved=1  ab_ok=1
```

### Interpretation

1. **λ=100 still significantly improves left_continuous** (2% → 0.26%) → candidate default reasonable.
2. **div_A_relative ~1.0 barely decreases with λ_φ** → **phi gauge and A divergence constraint are independent**; Stage 10.5 closed root cause, Stage 10.6 opens A gauge issue.
3. **Three-body TEAM7 λ sweep complete** (see §3b): plate Joule/J insensitive to λ_φ; **λ_φ=100 candidate default OK on three-body**.

---

## 3b. Three-body TEAM7 λ sweep (user run, GPU 1, 2026-05-28, after divA aggregation fix)

**Command**:

```bash
cd /home/yongchuan/sphinxsys/build
CUDA_VISIBLE_DEVICES=1 ./tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_contact_team7_phi_gauge_sweep_diagnostic/bin/test_3d_aphi_ck_contact_team7_phi_gauge_sweep_diagnostic
```

**stdout (`completed_rows=4 passed=1`)**:

| λ_φ | converged | outer | global_true_rel | max_bodywise_true_rel | plate_joule_gap | plate_j_L2_gap | div_A_relative | plate_div_A_relative | div_A_level |
|-----|-----------|-------|-----------------|----------------------|-----------------|----------------|----------------|----------------------|-------------|
| 10 | 1 | 10 | 3.633e-4 | 1.792e-4 | 0.029% | 0.017% | **0.20468** | 0.02312 | high_risk (global) |
| 30 | 1 | 9 | 3.877e-4 | 1.869e-4 | **0.005%** | **0.005%** | **0.20474** | 0.02311 | high_risk (global) |
| 100 | 1 | 9 | 3.106e-4 | 1.510e-4 | 0.026% | 0.010% | **0.20474** | 0.02311 | high_risk (global) |
| 300 | 1 | 9 | 2.746e-4 | 1.321e-4 | 0.034% | 0.015% | **0.20473** | 0.02310 | high_risk (global) |

**Interpretation**:

- **plate Joule / J_L2 gap all < 0.04%** → λ_φ∈[10,300] still insensitive to three-body observables → **λ_φ=100 candidate default OK** (Contact scaffold default already changed to 100).
- **global div_A_relative ≈ 0.205** (high_risk), four λ levels fluctuate only **±0.00006** → **φ penalty has almost no effect on divA**.
- **plate_div_A_relative ≈ 0.023** (warn, not high_risk) → global metric inflated by **air/coil large volume, small |∇A| but nonzero divA** regions; **plate region more physically representative**.
- Krylov: all converged, outer 9–10; increasing λ slightly reduces `global_true_rel` / `max_bodywise_true_rel`, weak dependency like observables.
- **Conclusion**: Stage 10.5 (φ gauge) can close on this case; **Stage 10.7 A-divergence penalty has clear motivation** (global divA does not decrease with λ_φ).

---

## 4. Relation to Stage 10.5

| Question | Stage 10.5 conclusion | Stage 10.6 supplement |
|------|-----------------|-----------------|
| Left body ~1% raw field | Main cause λ_φ too small | λ=100 can suppress to ~0.3% |
| Krylov vs field | Tier-1/2 should be separate | still holds |
| div A = 0? | not checked | two-body global ~1.0; three-body **global ~0.21**, **plate ~0.023** (warn) |
| Next step | calibrate λ_φ | large divA → **Stage 10.7 A-penalty**; projection later |

---

## 5. Execution roadmap (aligned with ChatGPT)

```text
✅ 1. λ=100 A/B (two-body)
✅ 2. divA post-solve diagnostic (reduction + DivAReal/Imag registration)
✅ 3. Three-body TEAM7 λ sweep + divA (user GPU1 rerun confirmed, passed=1)
🔄 4. Stage 10.7 Inner A-penalty prototype (fused vs debug passed)
⏸ 5. Inner/Contact λ_A sweep + whether divA decreases
⏸ 6. Three-body dp/σ/source sweep (after divA closure)
⏸ 7. Cold crucible four-body
```

**λ_φ=100**: **promoted to Contact scaffold default** (two-body MMS / diagnostics + TEAM7 canonical); this sweep data supports that choice.

---

## 6. Run

Long GPU tests (full paths):

```bash
cd /home/yongchuan/sphinxsys/build
CUDA_VISIBLE_DEVICES=1 ./tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_contact_team7_phi_gauge_sweep_diagnostic/bin/test_3d_aphi_ck_contact_team7_phi_gauge_sweep_diagnostic
```

Build + short tests locally; **MMS / full sweep / three-body sweep delegated to user GPU**.

```bash
cd build
ninja test_3d_aphi_ck_contact_phi_gauge_lambda100_ab_diagnostic \
      test_3d_aphi_ck_contact_phi_gauge_penalty_sweep_diagnostic \
      test_3d_aphi_ck_contact_team7_phi_gauge_sweep_diagnostic \
      test_3d_aphi_ck_contact_left_field_error_diagnostic
```

---

## 7. One-liner (for ChatGPT)

> Stage 10.6: two-body global divA~1.0; three-body global~0.21, plate~0.023. λ_φ does not change divA; three-body Joule/J insensitive to λ∈[10,300] → λ_φ=100 candidate default. **phi penalty ≠ div A=0**; Stage 10.7 should prioritize **plate/conductor region** divA, not global L2 alone.
