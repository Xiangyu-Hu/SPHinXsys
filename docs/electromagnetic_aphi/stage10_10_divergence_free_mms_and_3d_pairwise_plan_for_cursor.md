# Stage 10.10 Execution Plan: Divergence-Free MMS, 3D Pairwise divA Validation, and Contact A-Penalty Unblocking Conditions

This document summarizes judgments from the latest Stage 10.9 results and provides detailed execution plan for Stage 10.10 for Cursor. Goal: avoid blindly pushing Contact A-penalty or cold-crucible four-body; first solve two foundational issues:

1. Current `div-free MMS` degenerates, i.e. `K(A_exact) ≈ 0`, cannot serve as effective manufactured benchmark;
2. `3D pairwise divA` has large error on true 3D divergence-free fields, e.g. `div_A_rel≈0.48`; cannot yet prove pairwise divA is reliable Coulomb-gauge diagnostic.

Recommendations here should serve as next development mainline.

---

## 0. One-Sentence Conclusion

Stage 10.9 value is proving:

```text
Inner pairwise A-divergence penalty route "may work".
```

But not yet proving:

```text
It accurately satisfies div A = 0;
It improves E/J/Joule physical accuracy;
It can connect to Contact;
It can be production default.
```

Therefore next step is not Contact A-penalty directly, but first complete:

```text
Stage 10.10-A: non-degenerate divergence-free MMS;
Stage 10.10-B: 3D pairwise divA fix/explanation;
Stage 10.10-C: E/J/Joule exact/reference validation;
Stage 10.10-D: eta_A policy freeze;
Stage 10.10-E: then decide whether to unblock Contact A-penalty.
```

---

# Part I. Current Progress Positioning

## 1. Contact Baseline Still Stable

Contact mainline has not regressed; baseline should remain:

```text
lambda_phi = 100
lambda_A = off
coupled multi-body GMRES
two-body Contact MMS passed
three-body TEAM7-like scaffold passed
E/J/Joule observable available
```

This baseline can continue as engineering scaffold.

Note:

```text
lambda_phi = 100
```

Only controls `phi` gauge drift, not `divA`.

Still should not enable:

```text
Contact A-divergence penalty
```

Also should not enter:

```text
cold-crucible four-body
10D thermal coupling
```

---

## 2. Inner Pairwise A-Penalty Research Prototype Established

Stage 10.8/10.9 completed:

```text
pairwise divA
pairwise grad(divA)
pairwise 3×3 PC
eta_A relative scale
PC gate
Inner solenoidal gate
```

Key progress over Stage 10.7.

Stage 10.7 main issue:

```text
A-penalty apply: B-corrected pipeline
PC: pairwise 3×3 block
=> apply/PC inconsistent, GMRES residual floor
```

Stage 10.8/10.9 changed to:

```text
A-penalty apply: pairwise divA + pairwise grad(divA)
PC: pairwise 3×3 block
=> apply/PC same-family consistent
```

This shows:

```text
Inner pairwise A-divergence penalty research prototype can continue.
```

---

## 3. But Coulomb-Gauge Validation Not Complete

Cannot currently say:

```text
div A = 0 is accurately satisfied.
```

Reason:

```text
Linear2D divergence-free field:
    pairwise div_A_rel ≈ 1e-9, very good

Sinusoidal3DCurlPsi:
    B-corrected div_A_rel ≈ 0.011
    pairwise div_A_rel ≈ 0.48, very poor
```

This shows:

```text
pairwise divA works well on simple 2D / 2.5D-like fields;
but not yet validated on true 3D curl-generated fields.
```

Therefore current biggest technical risk:

```text
Whether 3D pairwise divA can serve as reliable Coulomb-gauge diagnostic.
```

If unresolved, Contact A-penalty even if running cannot explain physical meaning.

---

## 4. Case C Is Observable Sensitivity, Not Reference Validation

Current Case C shows:

```text
After enabling eta_A, E/J/Joule do not collapse by orders of magnitude;
eta_A=0.1/0.3 can lower current pairwise divA metric;
observable change within research-acceptable range.
```

But cannot show:

```text
E/J/Joule closer to true solution.
```

Reason:

```text
No analytic E/J/Joule reference currently;
No trustworthy fine-dp reference;
Existing dp=0.1 vs dp=0.075 self-reference gap very large, e.g. Joule gap up to O(90%).
```

So Case C should be written as:

```text
observable sensitivity gate
```

Not:

```text
physical accuracy validation
```

---

# Part II. Q1–Q6 Decisions for Cursor

## Q1. div-free MMS Degeneration: How to Construct `b != 0` with `divA = 0`?

### Conclusion

Do not continue using degenerate fields as main benchmark, e.g.:

```text
A = (-y/2, x/2, 0)
A = curl(0,0,psi) but K(A)≈0 construction
```

These may be harmonic / near-null modes causing:

```text
rhs_norm ≈ 0
```

MMS loses meaning.

### Recommended Field 1: 2.5D Divergence-Free Nonzero-RHS Field

Construct:

```text
A_x = 0
A_y = 0
A_z = sin(pi x) sin(pi y)
phi = 0
```

Analytically:

```text
div A = ∂A_z/∂z = 0
```

But:

```text
laplace A_z = -2 pi^2 sin(pi x) sin(pi y) != 0
```

Therefore:

```text
K(A_exact) should not be 0;
E = -j omega A nonzero;
Joule nonzero;
B = curl A nonzero.
```

Most recommended non-degenerate divergence-free MMS next.

### Recommended Field 2: True 3D Divergence-Free Field

Construct:

```text
A_x = sin(pi y) sin(pi z)
A_y = sin(pi z) sin(pi x)
A_z = sin(pi x) sin(pi y)
phi = 0
```

Analytically:

```text
∂A_x/∂x = 0
∂A_y/∂y = 0
∂A_z/∂z = 0
=> div A = 0
```

Also:

```text
laplace A_x = -2 pi^2 A_x
laplace A_y = -2 pi^2 A_y
laplace A_z = -2 pi^2 A_z
```

True 3D field; use as 3D validation case, but run Field 1 first.

### RHS Modes

For each field, two RHS modes.

#### Mode A: Consistency RHS

```text
b_eta = K_eta(A_exact, phi_exact)
```

Purpose:

```text
Verify penalized discrete system self-consistently solvable.
```

#### Mode B: Physical Invariance RHS

```text
b_0 = K_0(A_exact, phi_exact)
solve K_eta x = b_0
```

Purpose:

```text
If A_exact already div-free, penalty should not significantly change solution.
```

---

## Q2. 3D pairwise div_A_rel≈0.48: Must Fix, or Is B-Corrected 3D Enough?

### Conclusion

Must fix, or at least explain and downgrade.

If pairwise route is A-penalty mainline, giving on true 3D divergence-free field:

```text
div_A_rel ≈ 0.48
```

Cannot ignore.

Reason:

```text
If continuous theory div A = 0,
but pairwise diagnostic says divA large,
penalty will forcibly modify originally correct A field.
```

This makes A-penalty physical meaning unstable.

### Required Stage 10.10-B: 3D Pairwise divA Diagnostic

Add or extend test:

```text
test_3d_aphi_ck_pairwise_diva_3d_diagnostic
```

Test fields:

```text
1. Linear2D field
2. 2.5D field: A=(0,0,sin(pi x)sin(pi y))
3. 3D cross-sine field:
   A_x=sin(pi y)sin(pi z)
   A_y=sin(pi z)sin(pi x)
   A_z=sin(pi x)sin(pi y)
4. Existing Sinusoidal3DCurlPsi
```

Per field output:

```text
pairwise div_A_abs_L2
pairwise div_A_rel_A_norm
pairwise div_A_rel_gradA_norm
B_corrected div_A_abs_L2
B_corrected div_A_rel_A_norm
B_corrected div_A_rel_gradA_norm
core_shell sensitivity
dp sensitivity
boundary vs core contribution
```

### Judgment

If:

```text
pairwise 3D divA does not decrease with dp
```

Then pairwise divA cannot serve as full 3D Coulomb-gauge validation; only penalty regularization route.

If:

```text
pairwise improves significantly in core / larger core_shell
```

Then main issue is boundary support / normalization.

If:

```text
B-corrected stable good, pairwise stable bad
```

Then consider:

```text
pairwise penalty route as solver regularization;
B-corrected divA as diagnostic reference;
```

But must explain the difference.

---

## Q3. Case C: Can eta=0.1 Be Research Candidate?

### Conclusion

Yes.

But only say:

```text
eta_A = 0.1 is a research candidate for source-driven Inner experiments.
```

Cannot say:

```text
eta_A = 0.1 makes E/J/Joule more accurate.
```

### Current Case C Meaning

Case C shows:

```text
eta_A = 0.1 lowers current divA metric;
GMRES converges;
E/J/Joule vs eta=0 baseline no order-of-magnitude breakdown;
Milder than eta=0.2/0.3.
```

Therefore suggest:

```text
primary research candidate: eta_A = 0.1
optional: eta_A = 0.2
upper research bound: eta_A = 0.3
over-penalty diagnostic: eta_A = 1.0
```

### Still Missing

To prove "closer to true solution", need:

```text
analytic MMS exact E/J/Joule
external FEM/MFEM reference
trustworthy high-resolution convergence reference
```

Current fine-dp self-reference gap too large; not reliable reference.

---

## Q4. Should eta_A Default Be Frozen?

### Conclusion

Freeze as research defaults, not production defaults.

Suggest:

```text
primary_eta_a = 0.1
optional_eta_a = 0.2
upper_eta_a = 0.3
over_penalty_eta_a = 1.0
```

Code level can go in:

```text
AphiADivergencePenaltyResearchDefaults
```

But production/default still:

```text
use_a_divergence_penalty = false
```

i.e.:

```text
lambda_phi = 100
lambda_A = off
```

---

## Q5. Next Priority Ordering

### First Priority: Non-Degenerate Divergence-Free MMS

Do first:

```text
A = (0,0,sin(pi x)sin(pi y))
```

Target:

```text
rhs_norm != 0
divA_exact = 0
E/J/Joule exact computable
```

### Second Priority: 3D Pairwise divA Fix/Explanation

Compare:

```text
2.5D field
3D cross-sine field
Sinusoidal3DCurlPsi
pairwise vs B-corrected
dp / core_shell / denominator sensitivity
```

### Third Priority: Documentation Freeze

Update record:

```text
Inner pairwise A-penalty research prototype convergent;
absolute gauge validation still open;
Contact A-penalty paused;
eta_A=0.1 primary research candidate.
```

### Fourth Priority: Contact A-Penalty

Only after first two resolved.

Do not do Contact A-penalty now.

---

## Q6. Should External Wording Change?

Make more rigorous.

### Can Say

```text
Stage 10.9 establishes an Inner pairwise A-divergence penalty research prototype.
The operator, PC, and eta_A scaling are now algebraically consistent.
The prototype converges in the Inner setting and reduces the discrete pairwise divA measure.
eta_A=0.1 is the current primary research candidate.
```

Summary:

```text
Stage 10.9 established Inner pairwise A-divergence penalty research prototype.
Current operator, PC, and eta_A calibration are same-family consistent;
Converges in Inner setting and lowers pairwise divA metric;
eta_A=0.1 is current primary research candidate.
```

### Cannot Say

```text
Coulomb gauge has been accurately enforced.
E/J/Joule are validated against the true solution.
Pairwise divA is fully validated in 3D.
A-penalty is production ready.
Contact A-penalty can be enabled.
```

### Most Accurate Wording

```text
Stage 10.9 confirms numerical feasibility of the pairwise A-divergence penalty in Inner research tests, but full 3D divergence-free validation and reference-based observable validation remain open. Contact A-penalty should remain paused.
```

---

# Part III. Stage 10.10 Detailed Execution Plan

## Stage 10.10-A: Non-Degenerate Divergence-Free MMS

### 10.10-A.1 New Test

```text
test_3d_aphi_ck_inner_divergence_free_a_nonzero_rhs_mms
```

### 10.10-A.2 Exact Field 1: 2.5D

```text
A_x = 0
A_y = 0
A_z = sin(pi x) sin(pi y)
phi = 0
```

Analytic quantities:

```text
divA_exact = 0

B_x = ∂A_z/∂y = pi sin(pi x) cos(pi y)
B_y = -∂A_z/∂x = -pi cos(pi x) sin(pi y)
B_z = 0

E = -j omega A - grad(phi) = -j omega A
J = sigma E
q = 0.5 sigma |E|^2
```

### 10.10-A.3 Exact Field 2: 3D Cross-Sine

```text
A_x = sin(pi y) sin(pi z)
A_y = sin(pi z) sin(pi x)
A_z = sin(pi x) sin(pi y)
phi = 0
```

Analytic:

```text
divA_exact = 0

laplace A_x = -2 pi^2 A_x
laplace A_y = -2 pi^2 A_y
laplace A_z = -2 pi^2 A_z
```

### 10.10-A.4 RHS Modes

Mode A:

```text
b_eta = K_eta(A_exact, phi_exact)
```

Mode B:

```text
b_0 = K_0(A_exact, phi_exact)
solve K_eta x = b_0
```

### 10.10-A.5 Output

```text
rhs_norm
eta_A
true_rel
A_L2_error
A_Linf_error
phi_L2_error
E_L2_error
J_L2_error
Joule_error
div_A_pairwise_rel
div_A_Bcorrected_rel
observable_change_vs_eta0
```

### 10.10-A.6 Pass Criteria

Initial criteria:

```text
rhs_norm > 1e-10 or nonzero threshold per case scale
GMRES converged
A_L2_error does not noticeably worsen with eta
div_A_rel clearly smaller than source-driven case
E/J/Joule error explainable
```

---

## Stage 10.10-B: Pairwise 3D divA Diagnostic

### 10.10-B.1 New/Extend Test

```text
test_3d_aphi_ck_pairwise_diva_3d_diagnostic
```

### 10.10-B.2 Field List

```text
1. Linear2D
2. Az2D = (0,0,sin(pi x)sin(pi y))
3. CrossSine3D
4. Sinusoidal3DCurlPsi
```

### 10.10-B.3 Output

```text
field_name
dp
core_shell
pairwise_div_A_rel_A_norm
pairwise_div_A_rel_gradA_norm
Bcorrected_div_A_rel_A_norm
Bcorrected_div_A_rel_gradA_norm
pairwise_div_A_abs_L2
Bcorrected_div_A_abs_L2
boundary_fraction
core_fraction
```

### 10.10-B.4 Analysis Goals

Judge:

```text
whether pairwise 3D divA converges with dp;
whether error concentrates at boundary;
whether reasonable after denominator change;
whether B-corrected stably better than pairwise;
whether pairwise only penalty route not diagnostic reference.
```

---

## Stage 10.10-C: Observable Exact/Reference Validation

From 10.10-A exact fields, output analytic:

```text
E_exact
J_exact
Joule_exact
B_exact
```

Compare:

```text
eta_A = 0, 0.1, 0.2, 0.3
```

Output:

```text
E_error_vs_exact
J_error_vs_exact
Joule_error_vs_exact
B_error_vs_exact
observable_change_vs_eta0
```

Goal:

```text
Distinguish "penalty-induced change" from "physical error".
```

If:

```text
eta_A=0.1 lowers divA while E/J/Joule closer to exact
```

Penalty has physical benefit.

If:

```text
eta_A=0.1 lowers divA but E/J/Joule farther from exact
```

Penalty too strong or discretization issue.

---

## Stage 10.10-D: eta_A Policy Freeze

From 10.10-A/B/C results, freeze research policy:

```text
eta_A = 0.1 primary
eta_A = 0.2 optional
eta_A = 0.3 upper
eta_A = 1.0 over-penalty diagnostic
```

If 10.10-A/C shows `eta=0.1` still noticeably destroys exact observable, downgrade:

```text
eta_A = 0.03 or 0.05
```

As new primary candidate.

---

## Stage 10.10-E: Contact Unblocking

Only start Contact A-penalty when all conditions met:

```text
1. non-degenerate div-free MMS passed;
2. 3D pairwise divA either fixed or clearly scoped;
3. eta_A=0.1 exact observable error acceptable;
4. Inner documentation frozen;
5. production default still off.
```

If not met, continue pausing Contact A-penalty.

---

# Part IV. What Not to Do Now

Do not:

```text
1. Directly enable Contact A-penalty;
2. Directly sweep eta_A on three-body TEAM7;
3. Set eta_A=0.1 as production default;
4. Claim Coulomb gauge satisfied;
5. Use Case C to claim E/J/Joule more accurate;
6. Continue using old degenerate div-free MMS as main benchmark;
7. Keep tuning eta_A blindly without fixing 3D pairwise divA.
```

---

# Part V. Final Cursor Execution Checklist

## Task 1: Add Non-Degenerate Divergence-Free MMS

Suggested files:

```text
test_3d_aphi_ck_inner_divergence_free_a_nonzero_rhs_mms.cpp
diagnostics/aphi_divergence_free_nonzero_rhs_mms_helpers.h
```

Must support:

```text
field_type = Az2D
field_type = CrossSine3D
rhs_mode = consistency
rhs_mode = physical_invariance
eta_A = 0, 0.1, 0.2, 0.3
```

---

## Task 2: Add 3D Pairwise divA Diagnostic

Suggested file:

```text
test_3d_aphi_ck_pairwise_diva_3d_diagnostic.cpp
```

Must support:

```text
Linear2D
Az2D
CrossSine3D
Sinusoidal3DCurlPsi
```

Output pairwise and B-corrected divA metrics.

---

## Task 3: Add Observable Exact Reference

In MMS helper add:

```text
E_exact
J_exact
Joule_exact
B_exact
```

Output error vs exact, not only vs eta=0.

---

## Task 4: Update Stage 10.9 Documentation

In:

```text
CURSOR_APHI_STAGE10_8_PAIRWISE_DIVA_RECORD.md
```

Or add:

```text
CURSOR_APHI_STAGE10_10_DIVERGENCE_FREE_VALIDATION_RECORD.md
```

Record:

```text
Stage 10.9 proves numerical feasibility, not final Coulomb gauge validation.
Stage 10.10 targets non-degenerate div-free MMS and 3D pairwise divA credibility.
```

---

## Task 5: Continue Freezing Contact A-Penalty

Keep:

```text
use_a_divergence_penalty = false
```

Do not enable in Contact tests until Stage 10.10-E conditions met.

---

# Part VI. Final Decision

Current recommendation:

```text
1. Inner pairwise A-penalty research prototype can continue;
2. eta_A=0.1 is primary research candidate;
3. Contact A-penalty continues paused;
4. 3D pairwise divA≈0.48 is current biggest technical risk;
5. Next must do non-degenerate divergence-free MMS;
6. Need exact E/J/Joule reference; cannot only look at vs eta=0 change;
7. After 10.10-A/B/C decide whether Contact unblocked.
```

One sentence:

> Stage 10.9 proves pairwise A-penalty numerical feasibility, but not correct Coulomb-gauge enforcement. Stage 10.10 targets closing this foundation with non-degenerate div-free MMS and 3D pairwise divA diagnostic.
