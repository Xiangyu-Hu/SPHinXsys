# Stage 10.8 Metric Interpretation and Follow-up Validation Plan: Correct Understanding of true_rel, divA, E/J/Joule Δ

This document guides Cursor in continuing Stage 10 A--phi / A-divergence penalty mainline work. The core question is: how should we interpret the current Inner solenoidal gate numbers, for example

```text
eta_A = 0:
    true_rel ≈ 1.0e-5
    div_A_rel ≈ 0.54

eta_A = 0.1:
    true_rel ≈ 2.7e-5
    div_A_rel ≈ 0.20
    Joule Δ ≈ 13%
    E_L2 Δ ≈ 6.5%

eta_A = 0.3:
    true_rel ≈ 5.2e-6
    div_A_rel ≈ 0.10
    Joule Δ ≈ 14%
    E_L2 Δ ≈ 7.0%
```

What do these values mean? What should their "theoretically correct" values be? Are they merely "lower" but still far from the correct answer?

Conclusion first:

> The current results only show that the pairwise A-divergence penalty has numerical effect on the Inner solenoidal gate and that GMRES can converge; they do not yet show that `div A = 0` has reached theoretical correctness, nor that E/J/Joule are quantitatively correct.

Subsequently we must add manufactured divergence-free benchmark, dp convergence, absolute/reference error, and observable reference comparison to move this from research gate to validation gate.

---

## 1. What Each Metric Means

### 1.1 What is `true_rel`

`true_rel` is the relative value of the linear system residual, roughly:

```text
true_rel = ||K(x) - b|| / ||b||
```

It measures:

```text
Whether GMRES has solved the discrete linear system well.
```

Theoretically, if the linear system is solved exactly:

```text
true_rel = 0
```

Numerically it cannot equal 0, so we set tolerances, e.g.:

```text
strict tolerance: 1e-4 or 1e-5
engineering tolerance: on the order of 1e-3
```

Therefore:

```text
true_rel ≈ 1e-5 ~ 3e-5
```

means the linear system itself is solved very well.

But:

```text
small true_rel ≠ physics must be correct;
small true_rel ≠ divA must be small;
small true_rel ≠ E/J/Joule must be close to the theoretical solution.
```

It is only an algebraic metric that "the discrete equations were solved."

---

### 1.2 What is `div_A_rel`

In the current Stage 10.8 pairwise path, `div_A_rel` is defined similarly to:

```text
div_A_rel = ||div A|| / ||A||
```

or some normalized relative divergence metric.

It measures:

```text
How large the longitudinal / divergence component in the A field is.
```

If Coulomb gauge is strictly satisfied:

```text
div A = 0
```

then the theoretical target is:

```text
div_A_rel = 0
```

Numerically it should converge to a small value with resolution, and at least decrease with `dp` in a divergence-free manufactured solution.

Therefore the current:

```text
eta_A = 0:   div_A_rel ≈ 0.54
eta_A = 0.3: div_A_rel ≈ 0.10
```

only shows:

```text
The A-divergence penalty reduced the relative divergence under the current definition by about 5×.
```

But this does not show:

```text
divA has reached theoretical correctness.
```

Because `0.10` is still not a very small quantity. It is improvement, not final validation.

---

### 1.3 What are `Joule Δ` and `E_L2 Δ`

The current:

```text
Joule Δ ≈ 13% ~ 14%
E_L2 Δ ≈ 6.5% ~ 7.0%
```

without an external theoretical solution or FEM reference, usually mean:

```text
Change relative to the eta_A = 0 baseline.
```

That is:

```text
How much E/J/Joule observables changed after enabling the A-divergence penalty.
```

This is not "error relative to the theoretical solution" unless the test explicitly provides a theoretical/analytic reference.

Therefore the correct interpretation of these two values is:

```text
The A-penalty did not destroy E/J/Joule to order-of-magnitude error;
but it did change observables, and 13%–14% Joule change is not very small.
```

So "E/J/Joule not obviously destroyed" should more accurately be:

```text
E/J/Joule change is within research-gate acceptable range, but not validation-grade accuracy.
```

For formal validation, 13% Joule change is too large to claim "physical quantities are accurate."

---

## 2. What These Values Should Be If Computed Correctly

### 2.1 For `true_rel`

Theoretical correct value:

```text
0
```

Numerical target:

```text
<= solver tolerance
```

For example:

```text
strict: true_rel < 1e-4 or 1e-5
engineering: true_rel < 1e-3
```

Current:

```text
true_rel ≈ 1e-5 ~ 3e-5
```

is very good.

---

### 2.2 For `div_A_rel`

If the target is Coulomb gauge:

```text
div A = 0
```

then the theoretical correct value:

```text
0
```

Numerical targets should not be fixed by experience alone; they should be determined by:

```text
1. divergence-free manufactured A_exact;
2. dp refinement;
3. relaxed vs regular particles;
4. divA diagnostic with the same discretization;
5. observe whether div_A_rel converges with dp.
```

For the current solenoidal gate:

```text
eta_A = 0.3: div_A_rel ≈ 0.10
```

only shows "it decreased," not "close to theoretical correctness."

Ideally, in a divergence-free MMS, we should see:

```text
div_A_rel(dp) -> 0 as dp -> 0
```

and at the same resolution, after penalty:

```text
div_A_rel_with_penalty << div_A_rel_without_penalty
```

But whether it is sufficient ultimately needs calibration against exact/reference.

---

### 2.3 For `Joule Δ` / `E_L2 Δ`

If there is a theoretical solution or FEM reference, ideal targets are:

```text
E_error -> 0
J_error -> 0
Joule_error -> 0
```

If the current comparison is only against eta_A=0 baseline, the ideal value is not 0, because enabling penalty will change the solution; but we hope it does not excessively change physically relevant observables.

Two types of metrics are recommended:

```text
observable_change_vs_eta0:
    perturbation of penalty on baseline.

observable_error_vs_exact_or_reference:
    true physical error.
```

Current `Joule Δ≈13%–14%` should be written as:

```text
penalty-induced observable change
```

not:

```text
theoretical error
```

Reference must be added subsequently to judge whether penalty improved physical results or deviated from them.

---

## 3. Correct Evaluation of the Current Results

### 3.1 What we can say

```text
Pairwise A-divergence penalty in the Inner solenoidal gate can be solved to GMRES convergence;
eta_A=0.1~0.3 significantly reduces div_A_rel under the current definition;
E/J/Joule show no order-of-magnitude destruction;
PC and apply are now clearly more consistent than Stage 10.7.
```

### 3.2 What we cannot say

```text
divA has reached theoretical correctness;
Coulomb gauge is strictly satisfied;
Joule is physically correct;
eta_A=0.3 is production default;
13% Joule change can be ignored;
Inner gate is equivalent to real TEAM7 validation.
```

### 3.3 More accurate description

```text
Stage 10.8 demonstrates a convergent pairwise A-divergence penalty prototype.
It reduces the discrete divA measure substantially, but the absolute divA accuracy and physical observable accuracy still require manufactured/reference validation.
```

---

## 4. Current Biggest Gap

The current gap is not whether GMRES converges, but:

```text
Missing theoretical/reference target values.
```

Specifically missing:

```text
1. A manufactured A-field benchmark with divA_exact = 0;
2. divA_rel convergence curve vs dp;
3. E/J/Joule error relative to exact/reference;
4. Whether A-penalty improves or perturbs physical observables;
5. Whether convergence still holds under relaxed particles.
```

Without these, `divA reduction` can only serve as research evidence, not final validation.

---

# Part II. Validation Cursor Must Do Next

## 5. New Case A: Divergence-free A MMS benchmark

### 5.1 Goal

Construct a manufactured solution that analytically satisfies:

```text
div A_exact = 0
```

and verify:

```text
divA diagnostic gives results close to 0;
A-divergence penalty does not destroy exact divergence-free solution;
error converges with dp.
```

### 5.2 Recommended construction

Use:

```text
A_exact = curl Psi
```

because any curl field satisfies:

```text
div(curl Psi) = 0
```

For example optionally:

```text
Psi = (0, 0, psi)
psi = sin(pi x) sin(pi y) sin(pi z)
```

Then:

```text
A_x =  ∂psi/∂y
A_y = -∂psi/∂x
A_z =  0
```

Analytically:

```text
div A = 0
```

### 5.3 Test content

Add:

```text
test_3d_aphi_ck_inner_divergence_free_a_mms
```

or:

```text
test_3d_aphi_ck_inner_solenoidal_a_mms
```

### 5.4 Two RHS modes

Two versions are recommended.

#### Mode 1: consistency RHS

```text
b_eta = K_eta(A_exact, phi_exact)
```

i.e. RHS includes the current penalty operator.

Purpose:

```text
Verify solver can solve the discrete system with penalty.
```

#### Mode 2: physical invariance RHS

```text
b_0 = K_0(A_exact, phi_exact)
```

then solve:

```text
K_eta x = b_0
```

If `A_exact` is also nearly `divA=0` discretely, the penalty term should be near zero, so the solution should not change noticeably.

Purpose:

```text
Verify penalty does not noticeably perturb divergence-free physical solution.
```

### 5.5 Output metrics

```text
eta_A
true_rel
A_L2_error
A_Linf_error
phi_L2_error
div_A_L2_abs
div_A_rel
div_A_exact_error
E_L2_error
J_L2_error
Joule_error
observable_change_vs_eta0
```

### 5.6 Acceptance targets

For divergence-free MMS:

```text
div_A_rel should be significantly smaller than the current source-driven case's 0.10;
and decrease with dp refinement.
```

This is the benchmark that truly shows divA computation is correct.

---

## 6. New Case B: dp refinement for divA

### 6.1 Goal

Determine the theoretical numerical level of current divA_rel.

### 6.2 Resolution

At least three levels:

```text
dp = coarse
dp = medium
dp = fine
```

If current box dp=0.1, can use:

```text
dp = 0.15, 0.1, 0.075
```

or choose based on existing case cost.

### 6.3 Output

```text
dp
particle_count
eta_A
true_rel
A_L2_error
div_A_L2_abs
div_A_rel
div_A_reduction_factor
E_L2_error
Joule_error
```

### 6.4 Judgment

If correct, should see:

```text
div_A_rel decreases with dp
A/E/J/Joule error decreases or remains bounded
```

If div_A_rel does not decrease with dp, then:

```text
Current divA diagnostic or penalty discretization still has systematic error.
```

---

## 7. New Case C: reference observable benchmark

### 7.1 Goal

Judge whether `E/J/Joule` change is "improvement" or "perturbation."

Current:

```text
Joule Δ≈13%–14%
E_L2 Δ≈7%
```

is only change relative to eta=0 baseline; we do not know which is closer to the correct solution.

### 7.2 Optional references

By priority:

```text
1. analytic MMS exact E/J/Joule;
2. high-resolution eta=0 or eta=small self-reference;
3. external FEM/MFEM reference;
4. monolithic inner reference.
```

### 7.3 Output

Do not only output:

```text
Joule Δ vs eta=0
```

Also output:

```text
Joule_error_vs_exact
E_error_vs_exact
J_error_vs_exact
```

or:

```text
Joule_error_vs_reference
```

### 7.4 Acceptance

If eta_A=0.1/0.3 makes:

```text
divA decrease;
E/J/Joule closer to exact/reference;
GMRES converge;
```

then we can say A-penalty improved physical results.

If:

```text
divA decreases;
but E/J/Joule far from reference;
```

then penalty is too strong or wrong form.

---

## 8. New Case D: observable sensitivity gate

A research gate can be retained:

```text
A-penalty should not excessively perturb observables.
```

But thresholds should be clearer.

Recommendation:

```text
eta_A = 0.1:
    Joule change vs eta0 < 10% preferred
    E_L2 change vs eta0 < 5% preferred

eta_A = 0.3:
    Joule change vs eta0 < 15% research acceptable
    E_L2 change vs eta0 < 10% research acceptable
```

Current:

```text
eta_A=0.3:
    Joule Δ≈14%
    E_L2 Δ≈7%
```

So it can count as:

```text
research acceptable upper bound
```

but not production default.

Recommend setting default candidate first to:

```text
eta_A = 0.1 or 0.2
```

rather than directly `0.3`.

---

# Part III. Corrections to Existing Conclusions

## 9. How "Inner gate passed" should be rewritten

Original statement:

```text
Inner gate passed.
```

Recommended more rigorous version:

```text
Inner pairwise A-divergence penalty research gate passed for solver convergence and divA-reduction trend.
However, absolute divA correctness and physical observable accuracy still require manufactured/reference validation.
```

English summary:

```text
Inner pairwise A-divergence penalty has passed the research gate:
it converges and can reduce divA.
But absolute divA correctness and E/J/Joule physical accuracy still require analytic or reference validation.
```

---

## 10. How "E/J/Joule not obviously destroyed" should be rewritten

The original statement is easily misunderstood.

Recommended rewrite:

```text
E/J/Joule show no order-of-magnitude destruction, but 6%–14% change relative to eta_A=0 baseline.
This shows penalty has substantial impact on physical observables; exact/reference must judge whether this impact is improvement or perturbation.
```

This is more accurate than "not obviously destroyed."

---

## 11. How "divA reduction" should be rewritten

Recommendation:

```text
divA reduction is a necessary but not sufficient condition.
The final target is not just a lower divA than baseline, but a discretely convergent divA that approaches the theoretical value, which is zero for Coulomb gauge.
```

English summary:

```text
divA reduction is only a necessary condition, not sufficient. The final target is not merely lower than baseline, but approaching theoretical value 0 as resolution increases.
```

---

# Part IV. Cursor Execution Checklist

## 12. Task 1: Revise Stage 10.8 record wording

Update:

```text
CURSOR_APHI_STAGE10_8_PAIRWISE_DIVA_RECORD.md
```

Add section:

```text
Interpretation of divA and observable metrics
```

Clarify:

```text
true_rel:
    solver residual; target is 0; numerically judged by tolerance.

div_A_rel:
    Coulomb-gauge diagnostic; theoretical target is 0.
    Current 0.54 -> 0.10 is only reduction, not final accuracy.

Joule/E_L2 Δ:
    currently penalty-induced change vs eta0 baseline, not theoretical error.
```

---

## 13. Task 2: Add divergence-free A MMS

New test:

```text
test_3d_aphi_ck_inner_divergence_free_a_mms
```

Functionality:

```text
A_exact = curl Psi
div A_exact = 0
phi_exact can be 0 or smooth field
generate RHS
solve eta_A = 0, 0.1, 0.2, 0.3
output A/E/J/Joule/divA error
```

---

## 14. Task 3: Add dp refinement

Based on divergence-free A MMS:

```text
dp sweep
```

Output convergence table:

```text
dp
N
eta_A
true_rel
A_L2_error
div_A_rel
E_L2_error
Joule_error
```

---

## 15. Task 4: Rewrite observable gate

Split current gate into two types:

### 15.1 Research sensitivity gate

```text
observable_change_vs_eta0
```

To judge whether penalty perturbs too much.

### 15.2 Validation error gate

```text
observable_error_vs_exact/reference
```

To judge whether physical quantities are correct.

Currently only the first type exists; the second is not yet available.

---

## 16. Task 5: Reposition eta_A default candidates

Recommendation:

```text
eta_A = 0.1:
    primary candidate

eta_A = 0.2:
    optional middle value

eta_A = 0.3:
    upper research bound

eta_A = 1.0:
    over-penalty diagnostic only
```

Do not use `eta_A=0.3` directly as default.

---

## 17. Task 6: Do not enable A-penalty on Contact yet

Until completing:

```text
divergence-free MMS
dp convergence
observable reference
```

At least complete these on Inner before entering Contact A-penalty.

Otherwise even if divA decreases in Contact, we cannot judge whether it is correct.

---

# Part V. Current External Messaging

## 18. What we can say

```text
The pairwise A-divergence penalty can be solved robustly in the Inner setting and reduces the discrete divA measure by several times.
```

## 19. What we cannot say

```text
The Coulomb gauge is now accurately enforced.
The physical E/J/Joule solution is validated.
The A-penalty is production ready.
```

## 20. Most accurate statement

```text
Stage 10.8 establishes a convergent research prototype for pairwise A-divergence penalty. It demonstrates divA-reduction capability, but absolute gauge accuracy and physical observable accuracy require divergence-free MMS and reference-based validation.
```

---

# 21. Final Recommendations

Current route should be adjusted to:

```text
1. Do not use divA reduction alone as success criterion.
2. Do not treat E/J/Joule Δ vs eta0 as theoretical error.
3. Do divergence-free A MMS first.
4. Then dp convergence.
5. Then E/J/Joule exact/reference validation.
6. For eta_A, recommend 0.1/0.2 first; 0.3 as upper research bound.
7. Contact A-penalty remains deferred until Inner absolute correctness validation is complete.
```

In one sentence:

> The current pairwise A-divergence penalty has proven "can converge, can reduce divA," but not yet "reduced to the correct value, physical quantities more accurate." Next step must introduce divergence-free manufactured solution and reference observables, rather than only continuing source-driven reduction sweeps.
