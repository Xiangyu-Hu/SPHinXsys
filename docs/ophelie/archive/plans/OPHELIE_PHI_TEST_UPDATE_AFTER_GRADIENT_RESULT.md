# OPHELIE French Reduced: Post–Phi-Test Assessment and Next-Step Investigation Plan

## 0. Latest test conclusions (2026-06-02 rerun)

### 0.1 Gradient sign — ruled out

```text
test_3d_ophelie_phi_gradient_linear:
    mean_GradPhi_x = 0.994
    rel_err_x ≈ 0.57%
    passed = 1
```

**Do not change the `(phi_i - phi_j)` sign.**

### 0.2 Laplace/RHS consistency — new test, MMS discrete incompatibility

Added `test_3d_ophelie_phi_laplace_rhs_consistency` (n=216, dp=0.08, penalty=0):

| mms-source | eq_res_vol | eq_res_linf | passed |
|---|---|---|---|
| continuous-grad | 1.88 | 2.07 | 0 |
| discrete-grad | 2.71 | 2.27 | 0 |

**Conclusion: before solving, `L(phi_exact) ≠ RHS`. The issue is not GMRES but that PairwiseLaplace and PhiRhsFromASrc are not the same compatible discretization.**  
discrete-grad did not improve (slightly worse), so this is not merely a continuous-vs-discrete gradient issue.

### 0.3 manufactured_A full-chain — cannot reject solver at penalty=0

`test_3d_ophelie_phi_solve_manufactured_A` (penalty=0):

| mms-source | pre_eq_res_vol | phi_solver_rel_res | post_eq_res_vol | reconstructed_divJ_red | passed |
|---|---|---|---|---|---|
| continuous-grad | 1.88 | 15.67 | 10.92 | 0.60 | 0 |
| discrete-grad | 2.71 | 7.30 | 9.22 | 0.58 | 0 |

`pre_eq_res_vol` matches the consistency test; solver failure to reach an exact solution on MMS is expected behavior.

### 0.4 French literature — small eq_res but divJ does not improve; penalty sweep ineffective

`test_3d_ophelie_french_reduced --reload=1 --literature-mode` (dp=0.02, n_glass=16066):

| phi_gauge_penalty | phi_solver_rel_res | phi_eq_res_vol | reconstructed_divJ_red | phi_solver_passed |
|---|---|---|---|---|
| 0 | 1.86e-3 | 0.0126 | 0.737 | 0 |
| 1e-8 | 1.84e-3 | 0.0126 | 0.737 | 0 |
| 1e-6 | 1.84e-3 | 0.0126 | 0.737 | 0 |
| 1e-4 | 1.87e-3 | 0.0126 | 0.737 | 0 |
| 1e-2 | 2.39e-3 | 0.0120 | 0.737 | 0 |
| **1 (default)** | **4.3e-5** | **1.17e-4** | **0.738** | **1** |

**Conclusion:**

```text
1. Gradient sign is not a bug;
2. Real-case eq_res is small (penalty=1 ~1e-4); solver works normally;
3. reconstructed_divJ_red ≈ 0.74 unchanged across all penalties → gauge penalty is not the main cause of divJ degradation;
4. Bottleneck shifts to: solver-side (Laplace + RHS) and postprocess-side (GradPhi + DivKernel(J)) operator incompatibility;
5. manufactured MMS cannot serve as solver regression; only diagnoses LHS/RHS discrete consistency.
```

---

## 1. Why full-chain manufactured_A failure cannot directly prove the solver is wrong

Current manufactured test sets:

```text
phi_exact = cos(pi*x/Lx) cos(pi*y/Ly) cos(pi*z/Lz)
ASrcReal = -grad(phi_exact) / omega
```

In continuous theory:

```text
EImag = -grad(phi_exact) - omega*ASrcReal = 0
```

But discrete implementation may not be discretely exact for three main reasons.

### 1.1 Gauge penalty changes the exact equation

Default parameters may include:

```text
phi_gauge_penalty = 1.0
```

The solver actually solves something like:

```text
(L + lambda I) phi = RHS
```

But manufactured RHS is built for:

```text
L phi = RHS
```

without the extra term:

```text
lambda * phi_exact
```

Therefore whenever `lambda != 0`, `phi_exact` is not an exact solution of that linear system.

Next steps must include:

```text
--phi-gauge-penalty=0
--phi-gauge-penalty=1e-8
--phi-gauge-penalty=1e-6
...
```

or explicitly add the penalty term in manufactured RHS.

### 1.2 ASrc uses continuous gradient, not necessarily compatible with discrete RHS/Laplace

In the test:

```text
ASrcReal = -continuous_grad(phi_exact)/omega
```

But RHS is discretized by pairwise operator:

```text
RHS = -div_discrete(omega * sigma * ASrcReal)
```

LHS is another pairwise Laplace:

```text
LHS = Laplace_pairwise(phi)
```

Even if satisfied continuously, discretely we may not have:

```text
Laplace_pairwise(phi_exact) == RHS_from_continuous_grad(phi_exact)
```

A stricter test should build a discrete-compatible source:

```text
1. Set PhiExact;
2. Use current ComputeOphelieScalarPhiGradientCK to get GradPhiDiscrete;
3. Set ASrcReal = -GradPhiDiscrete / omega;
4. Assemble RHS;
5. Compare RHS and Laplace(PhiExact).
```

If still inconsistent, then PhiRhsFromASrc and PairwiseLaplace themselves are not a compatible pair.

### 1.3 reconstructed DivJImag and solver residual are not the same operator

Small solver residual means:

```text
PairwiseLaplace(phi) - RHS ≈ 0
```

But reconstructed divJ is:

```text
J = sigma(-GradPhi - omega A)
DivJ = DivKernel(J)
```

a different composite operator. Even with small solver-consistent continuity residual, `DivKernel(J)` may remain large.

Therefore literature-mode should split:

```text
phi_solver_passed:
    ||Lphi - RHS|| / ||RHS|| small

divJ_warning:
    whether reconstructed DivJImag improves

demo_passed:
    A/B/E/J/Q finite, Q>=0, P target OK
```

---

## 2. What not to do now

```text
1. Do not flip ComputeOphelieScalarPhiGradientCK sign;
2. Do not add self-induction yet;
3. Do not expand geometry;
4. Do not fail demo based on reconstructed DivJImag alone;
5. Do not interpret full-chain manufactured_A failure as Biot or French geometry failure.
```

---

## 3. Recommended next-step order

### Step A: update documentation

Change the old review judgment "GradPhi sign highly suspicious" to:

```text
Gradient linear test confirmed current ComputeOphelieScalarPhiGradientCK sign is correct;
current issue shifts to operator compatibility / gauge penalty / boundary support.
```

### Step B: fix manufactured_A test design

Add at least three modes to `test_3d_ophelie_phi_solve_manufactured_A`:

```text
--mms-source=continuous-grad
--mms-source=discrete-grad
--phi-gauge-penalty=...
```

Where:

```text
continuous-grad:
    ASrc = -grad_exact(phi_exact)/omega

discrete-grad:
    PhiExact -> ComputeOphelieScalarPhiGradientCK -> ASrc = -GradPhiDiscrete/omega
```

And output:

```text
||Laplace(phi_exact) - RHS|| / ||RHS||
||Lphi_solved - RHS|| / ||RHS||
E_phi / E_L0
DivJ_phi / DivJ_L0
```

### Step C: add Laplace/RHS consistency test

Add:

```text
test_3d_ophelie_phi_laplace_rhs_consistency
```

Flow:

```text
1. Set PhiExact;
2. Compute LaplacePairwise(PhiExact);
3. Build ASrc with continuous-grad or discrete-grad;
4. Compute PhiRhsFromASrc;
5. Compare Lphi_exact and RHS.
```

If consistency is poor, the φ equation LHS/RHS discretization is incompatible; fix operator consistency before tuning the solver.

### Step D: solver-consistent continuity residual

Output in French reduced and manufactured tests:

```text
eq_res = ||Lphi - RHS|| / ||RHS||
```

Also keep:

```text
reconstructed_divJ_red = ||Div(J_phi)|| / ||Div(J_L0)||
```

Do not output only `literature_passed`.

### Step E: penalty sweep

After test design is sound, sweep:

```text
phi_gauge_penalty = 0, 1e-8, 1e-6, 1e-4, 1e-2, 1
```

Record:

```text
eq_res
divJ_L2_red
E_phi/E_L0
P_raw
maxPhi
maxGradPhi
```

### Step F: lattice vs reload

Finally run:

```text
--skip-relax --literature-mode
--reload=1 --literature-mode
```

to assess whether boundary/support amplifies reconstructed DivJ.

---

## 4. Specific reply recommendations for Cursor

You can reply to Cursor directly:

```text
A. Agreed. gradient_linear already proved the current sign is correct; do not change (phi_i - phi_j).

B. full-chain manufactured_A failure does not directly mean gradient is wrong. More likely:
   1) gauge penalty=1 makes phi_exact no longer an exact solution;
   2) ASrc uses continuous grad, incompatible with pairwise Laplace/RHS discretization;
   3) solver residual and reconstructed DivJImag are not the same operator.

Next steps: do not add self-induction or expand geometry. Continue in this order:

1. Update docs; remove old "gradient sign suspicious" judgment;
2. Add phi_gauge_penalty=0/1e-8 modes to manufactured_A;
3. Add discrete-grad source mode:
   PhiExact -> ComputeOphelieScalarPhiGradientCK -> ASrc=-GradPhiDiscrete/omega;
4. Add Laplace/RHS consistency test:
   compare PairwiseLaplace(PhiExact) and PhiRhsFromASrc;
5. In French reduced and manufactured tests, output separately:
   eq_res=||Lphi-RHS||/||RHS||
   reconstructed_divJ_red=||DivJ_phi||/||DivJ_L0||;
6. Then run penalty sweep and lattice/reload comparison.

Client demo continues with --no-phi + 50kW scaling.
```

---

## 5. Client demo scope

Still valid for demo:

```text
--reload=1 --no-phi --target-power=50000 --state_recording=1
```

Show:

```text
analytic cylinder glass
relaxed/reloaded particles
multiloop coil source
Biot-Savart A/B
E/J
JouleHeat scaled to 50 kW
```

Do not show:

```text
literature-mode passed
divJ continuity verified
full OPHELIE reproduced
```

---

## 5.1 DivSigmaGrad operator unification (2026-06-02 implementation)

- LHS/RHS changed to `div(σ∇φ)` / `-ω·div(σA)`, consistent with GradPhi+DivJ postprocess
- MMS: `laplace_rhs` eq_res ~3e-7 ✅; `manufactured_A` passed, divJ_L2_phi ~0.01
- French literature: `reconstructed_divJ_red` **0.74 → 1.82**, `literature_passed=1`
- `phi_eq_res_vol ~ 0.57` (**GMRES stall floor**; 350 outer still ~0.57; not boundary-dominated: interior≈boundary≈0.56)
- GMRES preconditioner `max(legacy,div)` + CLI tuning wired; more iterations cannot lower eq_res
- Post-Jacobi refinement (mismatched diagonal) diverges; disabled
- CLI: `--phi-lhs-operator=div-sigma-grad|legacy-pairwise`

See [`OPHELIE_PHI_OPERATOR_AUDIT.md`](OPHELIE_PHI_OPERATOR_AUDIT.md).

---

## 6. Latest conclusion

```text
Gradient sign already ruled out (gradient_linear passed).

Laplace/RHS consistency confirmed MMS discrete incompatibility (eq_res_vol ~ 1.9–2.7);
manufactured_A cannot reject the solver.

French real-case: phi_eq_res_vol ~ 1e-4 (solver OK), but reconstructed_divJ_red ~ 0.74 unchanged;
penalty sweep (0 → 1) has almost no effect on divJ → gauge penalty is not the main divJ cause.

Current bottleneck: solver-side PairwiseLaplace+RHS and postprocess-side GradPhi+Div(J) are not the same compatible discretization.

Next steps (by priority):
1. ~~lattice vs reload literature-mode comparison~~ ✅ completed; divJ both ~0.74–0.79; boundary not main cause;
2. Study whether GradPhi and Laplace/RHS flux forms should unify (same pairwise kernel);
3. Do not flip gradient, do not add self-induction, do not expand geometry.
```

---

## 7. Completed Cursor implementation manifest

| Step | Content | Status |
|---|---|---|
| A | Update docs; remove gradient sign suspicion | ✅ |
| B | manufactured_A add `--mms-source` + `--phi-gauge-penalty` | ✅ |
| C | Add `test_3d_ophelie_phi_laplace_rhs_consistency` | ✅ |
| D | French/manufactured separate output for `eq_res` and `reconstructed_divJ_red` | ✅ |
| E | penalty sweep (0, 1e-8, …, 1) | ✅ run; divJ no improvement |
| F | lattice vs reload comparison | ✅ run; both divJ ~0.74–0.79 |

### Step F results (literature-mode, dp=0.02, penalty=1 default)

| particles | n_glass | phi_eq_res_vol | reconstructed_divJ_red | divJ_continuity |
|---|---|---|---|---|
| lattice | 16106 | 1.01e-4 | 0.788 | failed |
| reload | 16066 | 1.18e-4 | 0.738 | failed |

Boundary/support (lattice vs reload) has limited effect on divJ; both have small solver eq_res but reconstructed divJ does not improve.

---

## 8. Biot RHS solvability diagnostic (2026-06-02)

### Operator decoupling

- `OpheliePhiRhsOperatorKind`: `DivSigmaA` (`-ω·div(σA)`, consistent with DivSigmaGrad LHS) vs `LegacyFlux` (pairwise flux).
- `phi_rhs_operator_kind_` and `phi_lhs_operator_kind_` independent; literature-mode default LHS=`DivSigmaGrad`, RHS=`DivSigmaA`.
- CLI: `--phi-rhs-operator=div-sigma-a|legacy-flux`
- Test: `test_3d_ophelie_phi_biot_rhs_solvability` (French cylinder + multiloop Biot, `dp=0.06`, GMRES 80 outer)

### Results (n=562, no reload)

| Metric | div-sigma-a | legacy-flux |
|---|---|---|
| `||RHS_div - RHS_legacy||_vol / ||RHS_div||` | — | **~2.0** (large RHS difference) |
| `phi_eq_res_vol` | 0.205 | 0.205 (almost same, GMRES stall) |
| `divJ_L2_reduction` | **~4.87** | ~0.51 |

**Interpretation:**

1. Two RHS discretizations differ significantly on Biot field, but **eq_res floor almost identical** → real-case ~0.57 (reload) is more like "Biot A not in DivSigmaGrad range space" solvability floor, not a simple RHS formula typo.
2. **div-sigma-a RHS far better for divJ**; literature must keep `DivSigmaA`, not legacy-flux.
3. reload literature: `phi_eq_res_vol≈0.574`, `reconstructed_divJ_red≈1.74`, `literature_passed=1` (consistent with Option A).

### Cross-RHS residual (same φ, swap b)

| particles | n | `rhs_cosine` | div `eq_res` | legacy `eq_res` | div solution `cross_eq_res_legacy` | div `divJ_red` | legacy `divJ_red` |
|---|---|---|---|---|---|---|---|
| lattice dp=0.06 | 566 | **-1** | 0.294 | 0.294 | **1.93** | 3.39 | 0.52 |
| reload dp=0.02 | 16066 | **-1** | 0.574 | 0.574 | **1.79** | 1.74 | 0.56 |

**Interpretation:**

1. `rhs_cosine_div_legacy = -1`: two RHS nearly opposite in volume-weighted inner product, not a small discrete difference.
2. Each RHS after GMRES has almost same `eq_res` → stall floor independent of RHS choice.
3. φ from div solve on legacy RHS gives `cross_eq_res ≈ 1.79` (much larger than 0.57) → **single φ cannot satisfy both b**; on Biot field the equation under DivSigmaGrad discretization is closer to "b not in Ran(L)" than a typo.
4. div RHS still greatly improves divJ (reload: 1.74 vs 0.56); literature-mode fixed to `DivSigmaA` (CLI `--phi-rhs-operator` does not override literature default).

If both RHS give similar `eq_res` but large `cross_eq_res` → Biot A's b is not in L's range space, not merely an RHS typo.

**Run** (working directory must be `~/SPHinXsysSYCL/build`; do not use `...` for paths):

```bash
cd ~/SPHinXsysSYCL/build
cmake ..
cmake --build . --target test_3d_ophelie_phi_rhs_flux_sign_audit \
  test_3d_ophelie_phi_biot_rhs_solvability -j$(nproc)

# Sign audit (linear A, ~1s)
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_phi_rhs_flux_sign_audit/bin/test_3d_ophelie_phi_rhs_flux_sign_audit

# Biot RHS comparison — lattice
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_phi_biot_rhs_solvability/bin/test_3d_ophelie_phi_biot_rhs_solvability

# Biot RHS comparison — reload (needs build/reload/Reload.xml)
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_phi_biot_rhs_solvability/bin/test_3d_ophelie_phi_biot_rhs_solvability --reload=1

# If no Reload.xml, generate under build:
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_reduced/bin/test_3d_ophelie_french_reduced --relax=1
```

You can also run `--reload=1` from `.../test_3d_ophelie_phi_biot_rhs_solvability/bin/` (auto-finds `../../../../../reload`).

CSV appends to `./output/ophelie_phi_rhs_solvability.csv` (for matplotlib bar comparison).

---

## 9. Legacy flux vs div(σA) sign audit (2026-06-02)

### Continuous/discrete correspondence (constant σ)

| Form | Discrete implementation |
|------|----------|
| **DivSigmaA** | `div_i = Σ_j (σA)_j−(σA)_i · ∇W_ij`, RHS `= −ω·div` (`ComputeOphelieVecdDivergenceCK`) |
| **LegacyFlux** | `rhs_i −= ω σ_ij ∇W_ij·(A_i−A_j)` ≈ `+ω σ div(A)` (same uncorrected ∇W, σ harmonic mean) |

Hence **b_legacy ≈ −b_div** (not random): `cosine(div,legacy)=−1`, `||b_div−b_legacy||/||b_div||≈2`.

Measured (2026-06-02): `cosine(div,−legacy)=1`, `||b_div+b_legacy||_vol/||b_div||≈3e−7` (linear A and Biot reload) → **"opposite" on Biot is exact sign relation under same discretization, not approximate**.

### Tests

| Test | Purpose |
|------|------|
| `test_3d_ophelie_phi_rhs_flux_sign_audit` | Linear A(x); require `cosine(div,−legacy)>0.95`, `rel_diff(div,−legacy)<0.35` |
| `test_3d_ophelie_phi_biot_rhs_solvability` | Real Biot field; also prints `rhs_cosine_div_neg_legacy`, `sign_flip_ok` (diagnostic, not pass condition) |

### reload measured (user n=16066)

| Metric | Value |
|------|-----|
| `eq_res` div / legacy | 0.574 / 0.574 |
| `cross_eq_res` | ~1.79 |
| `divJ_red` div / legacy | **1.74** / 0.56 |
| `cosine(div,legacy)` | −1 |

**Conclusion:** legacy is not "another interchangeable RHS" but **old flux path sign-opposite to div**; `eq_res≈0.57` remains Biot solvability floor under L; delivery keeps **DivSigmaGrad + DivSigmaA**.
