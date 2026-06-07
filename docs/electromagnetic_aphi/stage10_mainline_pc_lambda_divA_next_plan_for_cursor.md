# Stage 10 Mainline Progress Notes: PC, λ_phi, λ_A, divA Discretization Choices and Follow-up Execution Plan

This document is for Cursor execution, clarifying several easily confused concepts in the current Stage 10 A--phi / Contact mainline, and giving the next mainline execution plan.

This round is not about teacher demo. Teacher demo is temporarily set aside; continue mainline:

```text
Keep Contact baseline;
Keep phi gauge penalty;
Continue A-side divA / Coulomb gauge work;
Shift focus from "continue sweeping λ_A" to "unify and validate div/grad/penalty/PC discretization forms first."
```

---

## 0. Conclusion First

Do not continue focusing on teacher demo. Next mainline should be:

```text
1. Keep Contact mainline: λ_phi = 100, λ_A = off.
2. Clarify λ_phi and λ_A are two completely different numerical regularization parameters.
3. Continue A-divergence penalty, but do not connect Contact yet.
4. Do not force B-corrected gradient; if B-corrected pipeline makes PC hard to close, can switch to more consistent, stable pairwise/direct divA discretization.
5. Key is not "must use which divergence formula," but:
   divA diagnostic, A-penalty apply, PC should use the same discretization definition as much as possible.
6. Do first on Inner: direct/pairwise divA + grad(divA) penalty + matching PC.
7. After Inner converges, then consider Contact.
```

---

# Part I. Concept Clarification

## 1. What is PC?

PC is **preconditioner**.

In current GMRES solve, we do not explicitly assemble the matrix, but compute via matrix-free apply:

```text
y = K(x)
```

GMRES needs to repeatedly solve linear system:

```text
K x = b
```

If `K` has poor condition number, GMRES is slow or residual stalls. Preconditioner `M` approximates inverse of `K`, or approximates main diagonal/block structure of `K`, so GMRES actually solves:

```text
M^{-1} K x = M^{-1} b
```

This concentrates the spectrum and makes GMRES easier to converge.

### 1.1 What is PC in current code

Currently mainly block Jacobi form:

```text
Construct local small block diagonal approximation per particle;
Locally approximate inverse for 8 unknowns of A/phi;
Used for GMRES preconditioning.
```

For A--phi block, 8 unknowns roughly are:

```text
A_real: 3 components
A_imag: 3 components
phi_real: 1
phi_imag: 1
```

So often called `8×8 PC` or block Jacobi.

### 1.2 Why A-divergence penalty needs new PC

Original operator mainly includes:

```text
A Laplace
phi Laplace
j omega sigma A reaction
sigma grad(phi)
j omega div(sigma A)
phi gauge penalty
```

Corresponding PC can already handle these terms.

But now added:

```text
A equation += λ_A * penalty_operator(A)
```

where penalty_operator is similar to:

```text
-grad(div A)
```

This operator is not ordinary scalar Laplace; it is longitudinal / component-coupled operator. In continuous Fourier space similar to:

```text
k k^T
```

It couples `A_x, A_y, A_z` three components. So if PC is still only:

```text
d_a += λ_A * laplace_diag
```

it cannot well describe new operator's spectral structure. Current GMRES stopping at 4%--8% residual floor is likely because new penalty operator and PC do not match.

---

## 2. What is λ_phi?

`λ_phi` is the coefficient of **phi gauge penalty**.

Current form is:

```text
LhsPhi += λ_phi * phi
```

i.e. add zero-order mass-like regularization term in `phi` equation.

### 2.1 What problem does λ_phi solve

It solves:

```text
phi constant gauge / null mode / raw phi drift
```

That is, `phi` can add a constant globally without noticeably affecting `grad(phi)`, causing:

```text
GMRES residual very good, but raw phi / raw field drifts greatly vs exact.
```

Previous sweep already proved:

```text
λ_phi = 0:
    Krylov residual can be very good;
    but left_continuous drifts to very large.

λ_phi = 100:
    raw field drift significantly reduced.
```

Therefore current Contact baseline using:

```text
λ_phi = 100
```

is reasonable engineering default.

### 2.2 What λ_phi does not solve

`λ_phi` does not solve:

```text
div A = 0
```

It is not Coulomb gauge. It only locks `phi` gauge, does not automatically reduce `divA`.

---

## 3. What is λ_A?

`λ_A` is the coefficient of **A-divergence penalty**, controlling divergence of A:

```text
div A ≈ 0
```

Target form is similar to:

```text
A equation += λ_A * [-grad(div A)]
```

If sign and discrete pairing are correct, this term penalizes `divA`.

### 3.1 λ_A and λ_phi have no direct relation

Very important:

```text
λ_phi and λ_A are not the same thing;
λ_phi cannot replace λ_A;
λ_A cannot replace λ_phi.
```

Differences:

| Parameter | Where added | Form | Effect | Order |
|---|---|---|---|---|
| `λ_phi` | phi equation | `λ_phi * phi` | Lock phi constant gauge | zero-order |
| `λ_A` | A equation | `λ_A * [-grad(div A)]` | Penalize A divergence | second-derivative |

Therefore cannot infer from:

```text
λ_phi = 100 reasonable
```

that:

```text
λ_A = 100 also reasonable
```

These parameters have different dimensions, scales, numerical meaning.

### 3.2 What does λ_A=100 mean

Current sweep's:

```text
λ_A = 100
```

is only an absolute numerical coefficient, weight of penalty term in A equation.

But this value itself lacks sufficient physical/discrete scale explanation. More reasonable calibration should use relative scale, e.g.:

```text
η_A = ||λ_A * graddiv_operator|| / ||LaplaceA_operator||
```

or:

```text
η_A = median(|λ_A * D_graddiv|) / median(|D_laplace|)
```

Therefore next step should not blindly sweep:

```text
λ_A = 10, 50, 100, 300, 1000
```

but first clarify:

```text
What is current λ_A's penalty_to_laplace_ratio.
```

---

# Part II. B-corrected Gradient Is Not Mandatory Route

## 4. If B-corrected gradient has issues, can skip it

Agreed.

Current `B-corrected` version advantages:

```text
Uses corrected gradient, theoretically better for irregular particles and boundary support.
```

But disadvantages:

```text
Composite pipeline:
    grad A -> trace -> grad(divA)
and PC hard to construct fully matching.
```

If current B-corrected pipeline causes:

```text
apply can run;
divA can decrease;
but GMRES residual floor very high;
PC hard to close;
```

then can temporarily switch to simpler, more consistent discretization form.

### 4.1 Current priority is not "must B-corrected"

Current priority should be:

```text
First construct a convergent, discretely consistent, energy sign correct A-divergence penalty prototype.
```

Not force:

```text
Must use B-corrected grad/div.
```

### 4.2 Allow trying direct / pairwise divA

Can try more direct SPH divergence:

```text
divA_i = Σ_j V_j (A_j - A_i) · ∇W_ij
```

or form consistent with current pairwise operator.

Then penalty:

```text
penalty_A_i = - grad(divA)_i
```

where grad(divA) also uses pairwise gradient matching divA:

```text
grad(divA)_i = Σ_j V_j (divA_j - divA_i) ∇W_ij
```

Key point:

```text
div operator D
grad operator G
penalty operator P = -G D
```

These three must be a set; cannot use corrected for divA diagnostic, another set for penalty apply, third set for PC.

---

## 5. trace divergence vs direct divergence

Current trace way is:

```text
gradA = gradient tensor of A
divA = trace(gradA)
```

Mathematically:

```text
div A = ∂A_x/∂x + ∂A_y/∂y + ∂A_z/∂z
```

So trace(gradA) is correct continuous formula.

But discretely there are two routes:

### 5.1 Route A: compute gradA first then trace

```text
A -> gradA -> trace(gradA)
```

Advantages:

```text
Can also use for curlA;
Unified with gradient tensor diagnostic;
Convenient to output gradA / divA / curlA.
```

Disadvantages:

```text
Composite operator complex;
After B-corrected gradient, PC hard to match.
```

### 5.2 Route B: compute divA directly

```text
A -> divA
```

e.g. pairwise:

```text
divA_i = Σ_j V_j (A_j - A_i) · ∇W_ij
```

Advantages:

```text
More direct;
Easier to pair with -grad(divA);
Easier to construct PC;
Easier to do energy test.
```

Disadvantages:

```text
Different from curlA gradient tensor route;
Need separate accuracy validation.
```

### 5.3 Recommendation

Next step should allow parallel comparison of two routes:

```text
Route 1:
    B-corrected gradA -> trace -> B-corrected grad(divA)

Route 2:
    pairwise/direct divA -> pairwise/direct grad(divA)
```

Compare:

```text
accuracy
energy sign
GMRES convergence
PC consistency
divA reduction
E/J/Joule observable impact
```

If Route 2 more stable, can use Route 2 first for A-divergence penalty prototype.

---

# Part III. How to Progress Current Mainline

## 6. Teacher demo paused, mainline continues

Agreed. Teacher demo can be set aside.

Current mainline priority should change to:

```text
1. Fix Contact baseline;
2. Continue A-divergence penalty Inner prototype;
3. Compare B-corrected vs direct/pairwise divA two discretizations;
4. Change to relative scale λ_A calibration;
5. Then consider Contact integration.
```

---

## 7. Freeze current baseline

First fix a stable version:

```text
Contact baseline:
    λ_phi = 100
    λ_A = off
    coupled multi-body GMRES
    two-body MMS passed
    three-body TEAM7 scaffold passed
```

Clarify in documentation:

```text
A-divergence penalty is research-only and disabled by default.
```

This protects current Contact results from A-penalty research path impact.

---

# Part IV. Cursor Next Execution Plan

## 8. Step 1: Concept and naming cleanup

### 8.1 Parameter naming

Keep:

```cpp
phi_gauge_penalty
```

For:

```text
LhsPhi += λ_phi * phi
```

Recommend naming A penalty clearly as:

```cpp
a_divergence_penalty
```

Or more accurately:

```cpp
a_coulomb_gauge_penalty
```

But documentation must state:

```text
λ_phi:
    phi gauge regularization coefficient.

λ_A:
    A divergence / Coulomb-gauge penalty coefficient.

They are independent numerical regularization parameters.
```

### 8.2 Documentation update

Update:

```text
CURSOR_APHI_STAGE10_7_A_DIVERGENCE_PENALTY_RECORD.md
```

Add section:

```text
Parameter clarification:
    λ_phi vs λ_A
    PC meaning
    B-corrected route vs direct pairwise route
```

---

## 9. Step 2: PC explanation and residual decomposition

Do not continue outputting only one `true_rel`.

Add residual block decomposition:

```text
res_A_real_norm
res_A_imag_norm
res_phi_real_norm
res_phi_imag_norm
res_A_total_norm
res_phi_total_norm
res_A_fraction
res_phi_fraction
```

If finer breakdown possible:

```text
res_A_longitudinal_proxy
res_A_transverse_proxy
```

Goal:

```text
Judge whether residual floor mainly from A penalty block, or phi coupling block.
```

---

## 10. Step 3: Compare two divA discretization routes

New test:

```text
test_3d_aphi_ck_div_a_discretization_comparison_diagnostic
```

Compare two divA:

### Route 1: B-corrected trace

```text
gradA_Bcorrected
divA = trace(gradA)
gradDivA_Bcorrected
```

### Route 2: direct pairwise

```text
divA_pairwise = Σ_j V_j (A_j - A_i) · ∇W_ij
gradDivA_pairwise = Σ_j V_j (divA_j - divA_i) ∇W_ij
```

Output:

```text
divA_B_L2
divA_pairwise_L2
divA_B_vs_pairwise_L2_diff
gradDivA_B_L2
gradDivA_pairwise_L2
gradDivA_B_vs_pairwise_L2_diff
energy_sign_B
energy_sign_pairwise
```

Purpose:

```text
Judge whether B-corrected pipeline necessary;
Judge whether pairwise/direct version more stable.
```

---

## 11. Step 4: energy sign diagnostic for both routes

For arbitrary A field, compute separately:

```text
P_B = -grad(divA)_Bcorrected
P_P = -grad(divA)_pairwise
```

Check:

```text
<A, P_B> / ||divA_B||²
<A, P_P> / ||divA_pairwise||²
```

Expect:

```text
Should be positive, magnitude close to 1 or at least stable positive value.
```

If one route's energy sign unstable or locally negative, that route unsuitable as penalty.

---

## 12. Step 5: First make convergent prototype with direct/pairwise route

If Route 2 energy sign and PC more consistent, recommend prioritizing:

```text
pairwise/direct A-divergence penalty prototype
```

Not because necessarily higher order, but because current stage needs one that is:

```text
Convergent;
Explainable;
Preconditionable;
Can reduce divA;
Does not destroy E/J/Joule;
Can migrate to Contact.
```

Tests:

```text
Inner impressed-current box
Inner solenoidal source
Inner manufactured div-free A
```

Output:

```text
converged
outer
true_rel
div_A_relative
res_A/res_phi decomposition
penalty_to_laplace_ratio
E/J/Joule if available
```

---

## 13. Step 6: Change λ_A to relative scale sweep

Do not sweep absolute values only.

Add output:

```text
penalty_to_laplace_diag_ratio
penalty_to_reaction_ratio
```

Recommend sweep:

```text
eta_A = 0, 0.01, 0.03, 0.1, 0.3, 1.0
```

If code still uses `λ_A`, back-calculate from local diagonal scale:

```text
λ_A = eta_A * median(|D_laplace|) / median(|D_graddiv|)
```

Or first roughly record ratio.

---

## 14. Step 7: PC route

### 14.1 If using B-corrected route

Need B-corrected/local FD-based PC:

```text
local finite-difference diagonal / 3x3 block of actual B-corrected penalty apply
```

### 14.2 If using pairwise route

Can construct more consistent pairwise 3x3 PC:

```text
D = direct div
G = direct grad
P = -G D
PC diagonal/block = local diagonal of P
```

Priority:

```text
apply and PC from same discrete formula.
```

Rather than pursuing highest-order corrected gradient from the start.

---

## 15. Step 8: Contact integration conditions

Only connect Contact when meeting:

```text
1. Inner penalty route converged;
2. true_rel meets target;
3. divA decreases on converged solution;
4. residual decomposition not abnormal;
5. E/J/Joule not destroyed;
6. λ_A has stable interval;
7. PC and apply consistent.
```

Contact version should do:

```text
A field
    -> Inner + Contact divA
    -> Inner + Contact grad(divA)
    -> A penalty contribution
```

Do not only do per-body Inner-only penalty, unless explicitly labeled:

```text
body-interior diagnostic only
```

---

# Part V. Answers to User's Current Questions

## 16. "What is PC?"

PC is preconditioner, helping GMRES solve linear system faster and more stably.

After adding current `λ_A grad(divA)`, system spectrum becomes harder, original PC does not match, so GMRES residual stops at 4%--8%. This is not solved by simply increasing iterations.

---

## 17. "What is λ_A=100, what is it for?"

`λ_A=100` is A-divergence penalty coefficient, adding divergence control term to A equation, making:

```text
div A -> 0
```

But current `100` is only absolute numerical sweep point, not proven reasonable scale. It has no direct relation to `λ_phi=100`.

---

## 18. "B-corrected gradient has issues, can skip it?"

Yes. Current priority is get it right, stable, convergent. B-corrected is not mandatory. If direct/pairwise divA + gradDivA more consistent, better preconditionable, can prioritize direct/pairwise route.

---

## 19. "trace divergence, can compute directly too?"

Yes. trace(gradA) is one way, direct pairwise divergence another. Key is:

```text
divA diagnostic;
A penalty apply;
PC;
```

Best use same discretization definition, avoid inconsistency.

---

## 20. "What's the relation between λ_phi and λ_A?"

No direct relation.

```text
λ_phi:
    zero-order gauge penalty in phi equation;
    controls phi constant drift.

λ_A:
    second-derivative divergence penalty in A equation;
    controls divA.
```

They are independent numerical regularization parameters, cannot replace each other.

---

# Part VI. Current Final Route Decision

## 21. Current immediate route

```text
1. Teacher demo paused.
2. Contact baseline frozen: λ_phi=100, λ_A=off.
3. A-divergence penalty continues Inner research.
4. Compare B-corrected route vs direct/pairwise route.
5. Prioritize discretely consistent, convergent route.
6. λ_A use relative scale calibration.
7. residual decomposition + PC consistency diagnostic must be done.
8. Connect Contact after Inner closes.
```

---

## 22. Do not do

```text
Do not connect Contact A-penalty now;
Do not use divA on unconverged solution as success evidence;
Do not treat λ_A=100 as calibrated default;
Do not conflate λ_phi=100 and λ_A=100;
Do not force B-corrected route;
Do not immediately switch to projection;
Do not enter cold crucible four-body / 10D.
```

---

## 23. Success criteria

For Inner A-divergence penalty prototype, success criteria should be:

```text
GMRES converged or true_rel reaches 1e-4~1e-3 engineering tolerance;
div_A_relative decreases on converged solution;
E/J/Joule not obviously destroyed;
residual decomposition no abnormal block;
λ_A relative scale sweep has stable interval;
PC and apply discretely consistent.
```

When not converged:

```text
divA reduction only shows directional trend, not success.
```

---

## 24. Final one-liner

Current issue is not Contact mainline stuck, but A-side Coulomb-gauge penalty discretization and preconditioning not yet closed.  
Next step: do not blindly sweep λ_A, do not force B-corrected gradient; compare corrected vs direct/pairwise two divA discretizations, choose route with correct energy sign, matching PC, GMRES convergent.
