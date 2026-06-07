# Stage 10.12 Execution Plan: Boundary divA, B=curlA, Contact A-Penalty Unblocking, and Follow-Up Mainline

This document is for Cursor execution. It summarizes judgments from the latest Stage 10.11 results and provides a detailed development roadmap for Stage 10.12. Focus areas:

1. Whether current 10.10 / 10.11 are complete;
2. How to interpret the boundary divA issue;
3. Whether `B=curlA` error is blocking;
4. Whether ghost/buffer enters production;
5. Whether Contact A-penalty can be unblocked;
6. Whether projection remains deferred;
7. Whether `eta_A=0.1` is frozen;
8. How the next batch of Stage 10.12 code should be written.

---

## 0. One-Sentence Conclusion

The current direction is correct. Stage 10.10 is largely closed, and Stage 10.11 further adds boundary divA specialization, B/E/J/Joule/H exact observables, ghost/buffer experiments, projection design, and `B=curlA` dual-track diagnostics.

The most important judgment right now is:

```text
The divA issue has narrowed from
"possible overall discretization error"
to
"core is correct, boundary is open, B post-processing needs dp refinement".
```

Therefore the next step can cautiously enter:

```text
Contact A-penalty two-body research gate
```

But must not directly enter:

```text
three-body TEAM7 A-penalty;
cold-crucible four-body;
10D thermal coupling;
production default A-penalty enabled;
projection implementation.
```

Current production / stable baseline remains:

```text
lambda_phi = 100
use_a_divergence_penalty = false
coupled multi-body Contact GMRES
two-body Contact MMS passed
three-body TEAM7-like scaffold passed
E/J/Joule observable available
```

---

# Part I. Current Stage Assessment

## 1. Is Stage 10.10 Complete?

Largely complete.

Stage 10.10 was originally meant to address:

```text
1. div-free MMS must not degenerate, i.e. rhs_norm must not be near 0;
2. whether pairwise divA is actually wrong in 3D divergence-free fields;
3. whether E/J/Joule can have exact/reference validation;
4. whether eta_A can serve as a research candidate;
5. whether Contact A-penalty can be unblocked.
```

According to current records, key 10.10 conclusions are:

```text
Az2D / CrossSine3D:
    rhs_norm is nonzero;
    consistency MMS passed;
    E_combined / J_combined / Joule exact can be validated;
    core-region pairwise divA can reach O(1e-7);
    large global divA mainly from boundary and denominator artifact.
```

Therefore we can say:

```text
10.10 non-degenerate MMS and core pairwise divA validation are largely complete.
```

But we cannot say:

```text
global Coulomb gauge is fully satisfied;
boundary divA is resolved;
Contact A-penalty is ready for production enablement.
```

---

## 2. Is Stage 10.11 Complete?

Main goals of Stage 10.11 were:

```text
1. boundary divA specialization;
2. B/E/J/Joule/H exact observables;
3. ghost/buffer experiment;
4. projection design;
5. B=curlA dual-track diagnostic.
```

Currently we can consider:

```text
10.11 is complete as a diagnostic stage;
but boundary divA and B=curlA refinement remain open issues.
```

The most important new results are:

```text
E_combined / J_combined / Joule vs exact ≈ 1e-6~1e-7;
B=curlA vs exact ≈ 3.5%;
boundary divA dominates global energy;
core divA remains very small.
```

This indicates:

```text
The A/phi -> E/J/Joule electrothermal source chain is already very accurate;
B=curlA as first-derivative post-processing still needs dp refinement;
boundary divA should still be managed as an open diagnostic.
```

---

# Part II. Correct Interpretation of Boundary divA

## 3. What Does Large Boundary divA Mean?

Current results show:

```text
core pairwise divA ≈ 1e-7
boundary divA dominates global divA energy
large global divA rel mainly from boundary and denominator artifact
```

Therefore large boundary divA must no longer be interpreted as:

```text
pairwise divA 3D discretization failed overall;
the entire A field is wrong;
A-phi solver is untrustworthy.
```

The more accurate interpretation is:

```text
core-region pairwise divA behaves correctly for non-degenerate divergence-free MMS;
the boundary layer contributes most divA energy due to incomplete kernel support, MMS/boundary-condition mismatch, unsuitable normalization denominator, etc.
```

This is consistent with common SPH boundary issues: near boundaries, kernel support is truncated and applying interior formulas directly loses consistency.

---

## 4. Can boundary divA ~0.03–0.04 Be Permanently Accepted?

No, not permanently.

But at the current stage it can be treated as:

```text
known SPH boundary/support diagnostic issue
```

Temporarily not blocking:

```text
Inner research prototype;
Contact two-body research gate;
E/J/Joule exact validation.
```

It must remain an open issue because:

```text
1. real Contact interfaces are also a kind of "boundary/interface";
2. cold-crucible wall/melt/crucible interfaces amplify boundary issues;
3. if target physical quantities concentrate near boundaries, boundary divA may affect results;
4. ghost/buffer is not universally effective currently.
```

Recommended external wording:

```text
Boundary divA is localized and currently treated as a support-deficiency diagnostic issue, not a core-operator failure. It is acceptable for the Inner research prototype, but not considered fully solved.
```

Summary:

```text
boundary divA is now classified as a boundary support/diagnostic issue, not a core-operator failure; it need not block the current Inner research prototype, but cannot be regarded as fully resolved.
```

---

## 5. Should Global `||divA|| / ||A||` Remain a Hard Gate?

It should not remain a hard gate.

Formally deprecate:

```text
global_div_A_rel_Anorm = ||divA|| / ||A||
```

as a pass/fail metric.

It may still be output, but only named:

```text
informational_global_Anorm_metric
```

New hard/diagnostic metrics should use:

```text
core_pairwise_divA_gradDen_rel
boundary_pairwise_divA_gradDen_rel
core_Bcorrected_divA_gradDen_rel
boundary_Bcorrected_divA_gradDen_rel
boundary_energy_fraction
```

Where:

```text
core metrics judge interior operator;
boundary metrics diagnose support/boundary;
B-corrected serves as a posteriori reference;
gradDen denominator is more suitable than ||A|| for measuring divA scale.
```

---

# Part III. Validation Status of B=curlA and E/J/Joule/H

## 6. Status of E/J/Joule Exact Validation

Current Az2D / CrossSine MMS already output:

```text
E_combined_error_vs_exact
J_combined_error_vs_exact
Joule_error_vs_exact
```

Representative results reach:

```text
~1e-6 to ~1e-7
```

This indicates:

```text
The A/phi solution recovery to E/J/Joule electrothermal source chain is very accurate.
```

Therefore we can currently say:

```text
E/J/Joule observable chain is validated for the current manufactured fields.
```

But note:

```text
this is still MMS manufactured-field validation;
not real TEAM7 or cold-crucible reference validation.
```

---

## 7. Is B=curlA ~3.5% Blocking?

It depends on context.

### 7.1 For the E/J/Joule Mainline

Not blocking.

Because E/J/Joule already have exact validation with very small error. For electromagnetic heating, Joule is the core heat source, and this chain is accurate.

### 7.2 For B-Field Validation

Blocking.

If claiming:

```text
B=curlA post-processing is validated
```

then `~3.5%` is not enough; dp refinement or improved gradient/curl discretization is needed.

Therefore the current wording should be:

```text
B=curlA is available as a diagnostic observable, but its current oscillatory-field error (~3.5%) requires dp-refinement before claiming high-accuracy B-field validation.
```

Summary:

```text
B=curlA is already available as diagnostic output, but current oscillatory-field error is about 3.5%; resolution convergence is needed before claiming high-accuracy B-field validation is complete.
```

---

## 8. Is B dp Refinement Needed?

Yes.

Recommended for Stage 10.12:

```text
dp = 0.15, 0.10, 0.075, maybe 0.05
```

At least for these fields:

```text
Az2D
CrossSine3D
Linear2D
```

Output:

```text
B_error_vs_exact
E_combined_error_vs_exact
J_combined_error_vs_exact
Joule_error_vs_exact
core_divA_gradDen
boundary_divA_gradDen
```

Judgment:

```text
If B_error decreases with dp:
    B=curlA post-processing can be said to converge, but current dp has ~3.5% error.

If B_error does not decrease:
    check B-corrected curl implementation, boundary/curl operator, field frequency, or kernel correction.
```

---

## 9. Are H and curlH Needed?

Yes, but not the current top priority.

Basic relations:

```text
B = curl A
H = nu B
curl H ≈ J_total
```

Do not write:

```text
curl B = H
```

That is incorrect.

For constant `nu`:

```text
H = nu B
```

is just scaling; once B is validated, H follows naturally.

For variable materials:

```text
H = nu(x) B
```

pay attention near interfaces.

`curlH` is a second-order derived quantity: curl of A first, then curl of H. It may be more sensitive than B, so wait for B dp refinement conclusions first.

Recommended Stage 10.12 ordering:

```text
B dp refinement first;
then H=nuB;
finally curlH≈J_total.
```

---

# Part IV. ghost/buffer and SPH Boundary Treatment

## 10. Should ghost/buffer Enter Production Diagnostics?

Not as production default for now.

Current ghost/buffer experiment results:

```text
Az2D / CrossSine3D:
    boundary rel worsens after ghost

Sinusoidal3DCurlPsi:
    boundary rel improves after ghost
```

This indicates:

```text
ghost/buffer is not a universal fix;
it helps some support-deficient fields;
it may over-correct fields with sufficient support or boundary-value mismatch.
```

Therefore current classification should be:

```text
conditional boundary diagnostic tool
```

Not:

```text
default production postprocessor
```

---

## 11. Should Contact Wall Dummy Be Done Now?

Not yet.

Real Contact / wall dummy introduces:

```text
Contact relation;
dummy state assignment;
multi-body variable sync;
interface support;
wall boundary condition;
operator consistency;
PC consistency.
```

Stage 10.12 step one should be Contact two-body A-penalty operator/PC gate, not wall dummy directly.

wall dummy / Contact boundary support can follow as:

```text
Stage 10.13 or later
```

---

# Part V. Should Projection Remain Deferred?

## 12. Projection Remains Deferred for Now

Recommendation:

```text
projection remains at design-document stage, not entering code mainline.
```

Reasons:

```text
1. core divA is already very good;
2. current main issue is boundary artifact;
3. projection may over-correct boundary artifact;
4. projection requires extra Poisson solve;
5. Contact projection chi boundary and interface conditions are complex;
6. phi synchronized update is a blocking issue;
7. current pairwise A-penalty route is still progressing and has not failed.
```

---

## 13. Is phi Synchronization a Blocking Issue for Projection?

Yes.

If doing projection:

```text
A_new = A - grad chi
```

To keep the electric field:

```text
E = -j omega A - grad phi
```

unchanged, usually need synchronized update:

```text
phi_new = phi + j omega chi
```

or the corresponding sign per current phasor convention.

Without phi update, projection changes E/J/Joule. That is not "pure gauge correction".

Therefore projection design must first answer:

```text
1. chi equation;
2. chi boundary conditions;
3. A update;
4. phi synchronized update;
5. whether E remains unchanged;
6. how to handle Contact multi-body.
```

Do not implement projection mainline until these are closed.

---

# Part VI. Current Positioning of eta_A

## 14. Can eta_A=0.1 Be Frozen?

Can be frozen as:

```text
primary research default
```

But not as production default.

Formally write:

```text
primary_eta_A = 0.1
optional_eta_A = 0.2
upper_eta_A = 0.3
over_penalty_eta_A = 1.0
```

production/default still keeps:

```text
use_a_divergence_penalty = false
```

i.e.:

```text
lambda_phi = 100
lambda_A = off
```

---

## 15. Is eta Sweep × Observable Matrix Still Needed?

Yes, but not to decide whether `eta_A=0.1` can be research default—that is already set—but to build a complete sensitivity table.

Recommended for Stage 10.12:

```text
Az2D / CrossSine3D
eta_A = 0, 0.1, 0.2, 0.3
```

Output:

```text
A_error
B_error
E_combined_error
J_combined_error
Joule_error
core_divA_gradDen
boundary_divA_gradDen
```

Purpose:

```text
confirm eta=0.1 is a mild choice;
confirm perturbation upper bound for eta=0.2/0.3;
determine future Contact initial eta.
```

---

# Part VII. Positioning of Invariance MMS

## 16. Should Invariance MMS Be Permanently Downgraded to Informational?

Yes. At minimum currently downgrade to:

```text
penalty sensitivity diagnostic
```

Not hard fail.

Reason:

```text
under b0 mode, continuous exact field is div-free;
but discrete boundary divA is nonzero;
penalty changes boundary;
solution change is expected.
```

So its role is:

```text
observe how much the solution changes relative to eta=0 after enabling eta_A.
```

Not to require:

```text
solution remains completely unchanged after enabling penalty.
```

---

# Part VIII. Can Contact A-Penalty Be Unblocked?

## 17. Can Contact A-Penalty Begin?

Can begin as **research prototype**, but strictly limited to two-body gate.

Must not directly enter:

```text
three-body TEAM7 A-penalty;
cold-crucible;
production default;
thermal coupling.
```

Why can it begin?

Minimum conditions are now met:

```text
core divA validated;
boundary issue localized;
E/J/Joule exact chain accurate;
eta_A=0.1 has research default;
projection not needed yet;
production baseline stable.
```

But Contact A-penalty still needs validation:

```text
1. Contact pairwise divA consistent with monolithic;
2. Contact gradDivA correct;
3. Contact PC includes cross-body neighbor contribution;
4. interface band divA controllable;
5. Contact A-penalty does not break two-body MMS / E/J/Joule.
```

---

## 18. What Is the First Contact A-Penalty Test?

Do not run three-body TEAM7 directly. First step:

```text
test_3d_aphi_ck_contact_a_divergence_penalty_pairwise_operator_diagnostic
```

or equivalent test.

### C1: Contact pairwise divA / gradDivA apply equivalence

Setup:

```text
monolithic body
split two-body Contact
same A field
```

Compare:

```text
divA
gradDivA
penalty contribution
```

Output:

```text
core_diff
interface_band_diff
boundary_diff
bodywise_diff
passed
```

### C2: Contact PC consistency gate

Similar to Inner PC gate, but must include Contact neighbor contribution.

Compare:

```text
finite-difference penalty apply diagonal
vs
Contact JacobiGradDivABlock
```

Output:

```text
core_pc_diff
interface_pc_diff
boundary_pc_diff
passed
```

### C3: two-body Contact A-penalty MMS

Settings:

```text
lambda_phi = 100
eta_A = 0, 0.1, 0.2
```

Output:

```text
converged
global_true_rel
max_bodywise_true_rel
interface_band_rel
E_combined_error
J_combined_error
Joule_error
core_divA
interface_divA
boundary_divA
```

Only after C1–C3 pass consider three-body.

---

# Part IX. Stage 10.12 Priority Ordering

## Priority 1: Document / External Wording Freeze

Freeze 10.10/10.11 rules first to avoid repeated misreading.

Must record:

```text
global Anorm divA not a hard gate;
core-only divA validated;
boundary divA open diagnostic;
B=curlA 3.5% informational, needs dp convergence;
E/J/Joule exact passed;
eta_A=0.1 research primary;
projection defer;
Contact A-penalty research-only unblocked.
```

---

## Priority 2: Contact A-Penalty Two-Body Research Prototype

This is the Stage 10.12 mainline push.

In order:

```text
C1: Contact pairwise divA / gradDivA apply equivalence
C2: Contact PC consistency gate
C3: two-body Contact A-penalty MMS
```

Do not go directly to three-body.

---

## Priority 3: B=curlA dp Refinement

To judge whether `~3.5%` converges.

Tests:

```text
Az2D
CrossSine3D
Linear2D
dp = 0.15, 0.10, 0.075, optional 0.05
```

Output:

```text
B_error_vs_exact
E_combined_error_vs_exact
Joule_error_vs_exact
core_divA
boundary_divA
```

---

## Priority 4: eta Sweep Observable Matrix

Build full table:

```text
Az2D / CrossSine3D
eta_A = 0, 0.1, 0.2, 0.3
```

Output:

```text
A_error
B_error
E_combined_error
J_combined_error
Joule_error
core_divA
boundary_divA
```

---

## Priority 5: curlH vs J Exact Validation

Recommended but not immediate blocking.

Start from simple fields:

```text
Linear2D or Az2D
```

Note:

```text
curlH = J_total
```

Clarify what `J_total` includes in the current case.

---

## Priority 6: Sinusoidal / Impressed / Cold Crucible

Continue as informational, not current gate.

---

## Priority 7: Projection Post-Solve Prototype

Continue defer; design document only.

---

## Priority 8: 10D Thermal Coupling

Continue defer.

---

# Part X. How Far From "Done"?

## 19. Inner A-Penalty Research

Current completion:

```text
~85%
```

Completed:

```text
pairwise route;
PC consistency;
eta policy;
non-degenerate MMS;
core divA validation;
E/J/Joule exact;
boundary diagnosis;
projection design.
```

Still missing:

```text
B dp convergence;
eta observable matrix;
boundary policy final documentation.
```

---

## 20. Contact A-Penalty

Current completion:

```text
~30%
```

Reason:

```text
Contact baseline mature;
but Contact A-penalty not yet truly integrated.
```

Next steps:

```text
Inner+Contact divA;
Inner+Contact gradDivA;
Contact PC;
two-body Contact gate;
then three-body TEAM7 eta sweep.
```

---

## 21. Production / True TEAM7 / Cold-Crucible Readiness

Current completion:

```text
~35%~45%
```

Still needed:

```text
Contact A-penalty;
B/curlH validation;
real TEAM7 reference;
real coil/source geometry;
cold-crucible four-body;
thermal coupling.
```

So the most accurate judgment is:

```text
We are on the right track and closer to closure than 10.10;
but not yet production / real TEAM7 / cold-crucible ready.
```

---

# Part XI. Cursor Execution Checklist

## Task 1: Update Stage 10.11 Record

Add summary:

```text
10.10/10.11 closure summary
```

Record:

```text
core divA validated;
boundary divA open;
global Anorm metric deprecated;
E/J/Joule exact passed;
B requires dp refinement;
ghost conditional only;
projection defer;
eta_A=0.1 research primary;
Contact research-only can start from two-body gate.
```

---

## Task 2: Add Contact A-Penalty Operator Diagnostic

Suggested file:

```text
test_3d_aphi_ck_contact_a_divergence_penalty_pairwise_operator_diagnostic.cpp
```

Suggested helper:

```text
aphi_contact_a_divergence_penalty_diagnostic_helpers.h
```

Function:

```text
monolithic vs split Contact
compare divA / gradDivA / penalty contribution
```

---

## Task 3: Add Contact PC Consistency Gate

Suggested file:

```text
test_3d_aphi_ck_contact_graddiv_pc_consistency_diagnostic.cpp
```

Function:

```text
finite-difference penalty apply diagonal
vs
Contact JacobiGradDivABlock
```

Must include:

```text
Inner neighbor contribution
Contact neighbor contribution
interface band particles
```

---

## Task 4: Add Two-Body Contact A-Penalty MMS

Suggested file:

```text
test_3d_aphi_ck_contact_a_divergence_penalty_two_body_mms.cpp
```

Parameters:

```text
lambda_phi = 100
eta_A = 0, 0.1, 0.2
```

Output:

```text
converged
global_true_rel
max_bodywise_true_rel
interface_band_rel
A_error
E_combined_error
J_combined_error
Joule_error
core_divA
interface_divA
boundary_divA
```

---

## Task 5: B=curlA dp Refinement

Extend existing:

```text
test_3d_aphi_ck_curl_b_dual_track_diagnostic
```

Add:

```text
dp sweep
```

Output:

```text
dp
N
B_error
E_error
Joule_error
core_divA
boundary_divA
```

---

## Task 6: eta Observable Matrix

Extend:

```text
test_3d_aphi_ck_inner_divergence_free_a_nonzero_rhs_mms
```

Add full matrix:

```text
field = Az2D, CrossSine3D
eta_A = 0, 0.1, 0.2, 0.3
```

Output:

```text
A/B/E/J/Joule errors
core/boundary divA
```

---

# Part XII. Final Decision

Current route recommendation:

```text
1. Formally freeze dual-track architecture:
   pairwise for penalty/PC;
   B-corrected + core/boundary for diagnostic.

2. Deprecate global ||divA||/||A|| hard gate.

3. boundary divA as open diagnostic; does not block Contact two-body research.

4. B=curlA 3.5% does not block E/J/Joule mainline, but blocks "B validated" claim.

5. Contact A-penalty can be unblocked as research prototype:
   phase one only two-body operator/PC/MMS gate.

6. Projection continues defer.

7. eta_A=0.1 fixed as research primary;
   production still off.

8. cold-crucible four-body and 10D thermal continue defer.
```

One sentence:

> You are on the right track. 10.11 has narrowed the divA issue from "possible overall discretization error" to "core correct, boundary open, B post-processing needs convergence validation". Stage 10.12 should begin Contact A-penalty two-body research gate, while in parallel doing B dp refinement and eta observable matrix.
