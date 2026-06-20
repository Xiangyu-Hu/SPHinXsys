# French Reduced OPHELIE 04.06-01 Code and Test Review: Blockers, Risks, and Next-Step Roadmap

## 0. Review Conclusion

The current `test_3d_ophelie_french_reduced` already completes a fairly full engineering chain:

```text
analytic cylinder GlassBody
    -> SYCL-CK relax / reload
        -> multiloop filament CoilSource
            -> host Biot-Savart A/B
                -> Level0 E/J/Q
                    -> optional PhiImag solve
                        -> power calibration / scaling
                            -> VTP output + divJ diagnostics
```

This means: **the customer-demo chain “external coil → A/B/E/J/Q in glass → JouleHeat” is basically presentable**, but it must be labeled **demo mode / source-to-JouleHeat chain**, not “OPHELIE literature current-continuity validation is already complete.”

The current real blocker is:

```text
phi solver residual is already small, but post-phi divJ is worse than Level0.
```

Logs show:

```text
phi rel residual ~8.3e-5
divJ_L2_L0 = 123.102
divJ_L2_phi = 166.907
divJ_L2_red = 0.737547
literature_passed = 0
```

This is not mainly a Biot-Savart or geometry-generation issue; it is a **PhiImag correction discrete consistency / sign / gradient-divergence diagnostic / boundary support** issue.

---

## 1. Current Progress: What Is Already Correct

### 1.1 Geometry generation direction is correct

The current French reduced case already avoids STL and uses:

```text
GlassBody:
    analytic cylinder

CoilSource:
    multiloop circular filament line source

Optional visual bodies:
    CoilVisualBody / CrucibleWallVisualBody
```

This is the correct direction. French literature does not provide a full CAD dimension table; the current reduced geometry should not chase exact CAD, but should validate the parametric chain.

### 1.2 relax / reload is already wired

The current test supports:

```text
--relax=1
--reload=1
--skip-relax
--reload-dir
```

In code:

```cpp
glass_body.generateParticles<BaseParticles, Lattice>();
relaxFrenchReducedBodies(...);
writeFrenchReducedReload(...);
```

And reload:

```cpp
glass_body.generateParticles<BaseParticles, Reload>(glass_body.Name());
```

This means that even with an analytic cylinder, body-fitted relax can be done as requested. The engineering flow here is correct.

### 1.3 Source -> Biot-Savart is already fairly clean

`electromagnetic_ophelie_current_moment_source.h` already abstracts the source as:

```text
CurrentMoment = A·m
```

And supports:

```text
volume:   J * dV
filament: I * dl
```

In Biot-Savart, the B direction is:

```cpp
moment.cross(r)
```

This is the correct direction. `multiloop_source` also has on-axis Bz analytic comparison.

### 1.4 Demo mode can run end-to-end

Logs show demo mode:

```text
phi_correction=0
P_raw=4.77409 W
P_scaled=50000 W
field_scale=102.339
passed=1
```

This means that for customer presentation only:

```text
coil excitation
    -> A/B
        -> E/J
            -> JouleHeat
                -> scaled 50 kW heat source
```

The basic conditions are already in place.

---

## 2. What Cannot Be Claimed Yet

Currently you cannot say:

```text
1. French OPHELIE is fully reproduced;
2. phi correction guarantees divJ≈0;
3. literature-mode has passed;
4. current results can serve as quantitative literature validation.
```

Because literature-mode explicitly failed:

```text
literature_passed=0
divJ_L2_red=0.737547 < gate 1.25
```

That is, current phi correction does not improve current continuity; it makes the `DivJImag` metric worse.

---

## 3. Main Blocker: Small Phi Residual but Worse divJ

### 3.1 Phenomenon

Current logs:

```text
phi residual:
    rel_res_l2_vol ≈ 8.3e-5
    rel_res_linf ≈ 2.3e-4 ~ 3.0e-4

but divJ:
    Level0 divJ_L2 = 123.102
    Phi divJ_L2 = 166.907
    divJ_L2_red = 0.737547
```

This means the linear system is solved under the current definition, but the computed `PhiImag` does not make post-processed `JImag` more divergence-free.

### 3.2 What this usually means

This kind of behavior is generally not “GMRES did not converge,” but one of:

```text
A. RHS/LHS sign of the phi equation is inconsistent with E = -grad(phi) - omega A;
B. ComputeOphelieScalarPhiGradientCK gradient sign may be flipped;
C. Laplace operator and post-process gradient/divergence are not the same compatible discretization;
D. divJ diagnostic operator and phi equation residual are inconsistent;
E. relaxed cylinder boundary support deficiency makes boundary divJ dominate;
F. gauge penalty changes the original Neumann problem and prevents physical divJ from dropping;
G. only A_src is used, not A_ind; self-consistent OPHELIE is not done yet, so physics is still incomplete.
```

Highest priority is A/B/C/D, i.e. **sign and discrete consistency**.

---

## 4. Code-Level Main Suspects

## 4.1 `ComputeOphelieScalarPhiGradientCK` sign — **ruled out (2025-06 update)**

`test_3d_ophelie_phi_gradient_linear` results:

```text
PhiImag = x
mean_GradPhi_x = 0.994
rel_err_x ≈ 0.57%
passed = 1
```

**Conclusion: current kernel outputs +∇φ, consistent with the (f_j - f_i)·g_ij convention used in divJ diagnostics. Do not flip (phi_i - phi_j).**

Current issue shifts to:

```text
- PairwiseLaplace(phi) and GradPhi + Div(J) are discretely incompatible;
- gauge penalty changes the exact MMS system;
- continuous vs discrete gradient construction compatibility for ASrc;
- solver eq_res and reconstructed DivJImag should be reported separately.
```

~~Original “GradPhi sign highly suspicious” judgment is obsolete.~~

---

## 4.1-legacy (obsolete) `ComputeOphelieScalarPhiGradientCK` sign needs priority check

Current code:

```cpp
grad_phi += g_ij * (phi_i - phi_j);
```

But the more common SPH gradient form is:

```cpp
grad_phi += (phi_j - phi_i) * gradW_ij * V_j;
```

Current `pairwiseGradientWeightUncorrected()` is already used in divergence diagnostics as:

```cpp
div_i += (vec_j - vec_i).dot(grad_W_ij);
```

If the `grad_W_ij` convention here is consistent with divergence, the compatible gradient form is more likely:

```cpp
grad_phi += (phi_j - phi_i) * g_ij;
```

rather than:

```cpp
grad_phi += (phi_i - phi_j) * g_ij;
```

If the current gradient is actually `-grad(phi)`, then:

```cpp
EImag = -GradPhiImag - omega * ASrcReal;
```

becomes:

```text
E = +true_grad(phi) - omega A
```

which directly flips phi correction direction and makes divJ worse.

### Required test

Add:

```text
test_3d_ophelie_phi_gradient_linear
```

Setup:

```text
PhiImag = x
```

Check:

```text
mean GradPhiImag_x should be +1
```

If output is -1, the current gradient sign is flipped.

This test is more important than continuing dp / penalty sweeps.

---

## 4.2 `OpheliePairwiseLaplaceCK` and gradient-divergence incompatibility

Current phi solve uses pairwise Laplace operator:

```cpp
laplace_i += pair_weight * (phi_i - phi_j);
```

But divJ diagnostic is:

```cpp
div(J)_i = sum_j (J_j - J_i) dot g_ij
```

And `J` is assembled after another `ComputeOphelieScalarPhiGradientCK` computes `gradPhi`.

This means:

```text
small phi residual
```

only means:

```text
PairwiseLaplace(phi) + penalty phi ≈ RHS
```

but does not guarantee:

```text
div( sigma * (-gradPhi - omega A) ) ≈ 0
```

because this div comes from “gradient operator + divergence operator,” not the same pairwise flux operator used in the phi solve.

### Two divJ metrics should be added

Recommended split:

#### Metric 1: equation residual continuity

Use phi equation residual directly:

```text
R_eq = L_phi - RHS
```

or by sign definition:

```text
R_eq = PairwiseLaplace(phi) + penalty*phi - PhiRhs
```

This is solver-consistent discrete continuity.

#### Metric 2: postprocess divJ

Current `DivJImag`, i.e.:

```text
div( reconstructed J )
```

This is a physics/visualization diagnostic, but should not alone decide whether the solver is correct.

Current literature-mode uses `DivJImag` as a hard gate, which may be too strict and may misjudge when operators are incompatible.

---

## 4.3 Phi RHS sign has not been truly validated

Already present:

```text
test_3d_ophelie_phi_rhs_constant_A
```

But constant A has theoretical RHS 0; it only proves “no garbage,” not correct sign.

Already present:

```text
test_3d_ophelie_phi_reduces_divJ
```

But it does not solve phi; it manually sets:

```cpp
grad_phi_imag = -omega * a_src_real
```

So it does not validate:

```text
RHS -> solve Phi -> compute GradPhi -> reduce divJ
```

as a full chain.

### Required test

Add:

```text
test_3d_ophelie_phi_solve_linear_A
```

Setup:

```text
AReal = (alpha*x, 0, 0)
constant sigma
box or cylinder interior
```

Flow:

```text
1. setup RHS
2. solve PhiImag
3. compute GradPhiImag
4. compute E/J
5. compute divJ
```

Acceptance:

```text
divJ_phi < divJ_level0
```

If this test fails, it is not a French geometry issue, but a sign/discrete issue in the phi solve chain itself.

---

## 4.4 gauge penalty should be optional, not default physical gate

Current:

```text
phi_gauge_penalty = 1.0
```

It makes the system easier to solve, but actually changes the Neumann null-space problem. For reduced cylinder, if penalty is too strong, it may suppress correct phi amplitude and prevent divJ from dropping.

Recommended sweep:

```text
--phi-gauge-penalty=0
--phi-gauge-penalty=1e-8
--phi-gauge-penalty=1e-6
--phi-gauge-penalty=1e-4
--phi-gauge-penalty=1e-2
--phi-gauge-penalty=1
```

Long term, more reasonable:

```text
pin one particle phi=0
```

or:

```text
subtract volume-weighted mean after each iteration/solve
```

rather than a large penalty.

---

## 5. Literature-Level Key Reminder: Still Not Full OPHELIE

French OPHELIE is not only:

```text
A = A_src
```

Literature model writes:

```text
V and j are unknowns;
A is given by Biot-Savart law over current density j;
matrix is full, non-symmetric, complex.
```

So strict OPHELIE is a self-consistent integral problem:

```math
A = A_src + K[j]
j = -sigma(grad V + i omega A)
div j = 0
```

Current French reduced case is mainly:

```text
A ≈ A_src from multiloop filament coil
```

That is, `self_induction` is not yet wired into the French case. Literature OPHELIE also treats metal conductor elements with surface/skin depth, but molten glass is a volume mesh. Ignoring metal-wall EM in the current simplification is acceptable, but if claiming “close to OPHELIE,” be aware:

```text
Stage 4 PhiImag correction is not full OPHELIE;
Stage 5 self-induction is closer to paper integral equation.
```

Therefore current literature-mode acceptance should not prematurely require literature agreement; call it:

```text
OPHELIE-inspired reduced A_src + phi mode
```

---

## 6. Can Customer Demo Happen Now?

### 6.1 Presentable version

Can present:

```text
Demo mode:
    --reload=1 --no-phi
    power scaling to 50 kW
```

Presentation content:

```text
1. code-generated cylindrical glass geometry;
2. body-fitted relax/reload;
3. external multiloop coil source;
4. Biot-Savart A/B;
5. induced E/J;
6. JouleHeat distribution;
7. total JouleHeat power normalized to 50 kW;
8. VTP/ParaView visualization.
```

This already answers the customer’s main chain:

```text
external coil excitation -> conductive glass EM response -> Joule heat source
```

But presentation must say:

```text
This is a reduced demo, not full OPHELIE literature reproduction.
Current phi correction / divJ continuity is still under validation.
```

### 6.2 Version not recommended for presentation

Do not present:

```text
--literature-mode
```

Because it currently has `literature_passed=0`, and post-phi divJ is worse. Do not show this to the customer.

### 6.3 If customer only cares about pipeline

If the customer only wants a “pipeline works” video/plot:

```text
can present.
```

### 6.4 If customer wants “literature alignment / OPHELIE validation”

```text
cannot present as completed.
```

---

## 7. Next-Step Roadmap

## P0: Fix phi chain first; do not expand geometry

Current geometry, relax, Biot, Level0 are sufficient. Priority now:

```text
1. test_3d_ophelie_phi_gradient_linear
2. test_3d_ophelie_phi_solve_linear_A
3. solver-consistent divJ/equation residual diagnostic
4. phi_gauge_penalty sweep
5. lattice vs relax comparison
```

### Recommended commands

Run first:

```bash
# lattice no relax
test_3d_ophelie_french_reduced --skip-relax --literature-mode --state_recording=1

# relaxed reload
test_3d_ophelie_french_reduced --reload=1 --literature-mode --state_recording=1
```

Compare:

```text
divJ_L2_red
phi_res
max PhiImag
max GradPhiImag
```

If lattice is good and reload bad, relax/boundary support is the main cause.  
If both are bad, prioritize sign/discrete consistency.

## P1: Split literature acceptance into two tiers

Current `literature_passed` is too strong; recommended split:

```text
demo_passed:
    A/B/E/J/Q finite
    Q>=0
    P_scaled=target

ophelie_phi_solve_passed:
    phi residual OK
    equation residual continuity OK

ophelie_physical_divJ_warning:
    reconstructed divJ improved or not
```

So customer demo is not failed by a single `DivJImag` gate.

## P2: Consider self-induction, but not immediately

Only after:

```text
Phi solve chain sign / gradient / residual diagnostic is fixed
```

then add:

```text
A_ind = K[J]
```

Otherwise self-induction amplifies current phi issues and prevents localization.

---

## 8. Explicit Tasks for Cursor

Can send directly to Cursor:

```text
Current French reduced geometry / relax / multiloop Biot / Level0 E/J/Q can serve as a demo chain, but literature-mode failed because PhiImag correction did not improve divJ. Next steps: do not expand geometry, and do not add self-induction first. Please investigate in this order:

1. Add test_3d_ophelie_phi_gradient_linear:
   Set PhiImag=x, check whether ComputeOphelieScalarPhiGradientCK outputs +ex.
   If -ex, fix gradient (phi_i - phi_j) / (phi_j - phi_i) sign.

2. Add test_3d_ophelie_phi_solve_linear_A:
   Set ASrcReal=(alpha*x,0,0), run full RHS -> solve Phi -> GradPhi -> E/J -> divJ.
   This test must prove the phi solve chain itself lowers divJ.

3. Add solver-consistent continuity residual:
   Do not only look at reconstructed DivJImag.
   Also output equation_residual = ||PhiLhsImag - PhiRhsImag|| / ||PhiRhsImag||,
   and continuity residual based on pairwise flux operator.

4. For French reduced case, compare lattice vs relax:
   --skip-relax --literature-mode
   --reload=1 --literature-mode
   If reload is clearly worse, body-fitted relax / boundary support is a main cause.

5. Run phi_gauge_penalty sweep:
   0, 1e-8, 1e-6, 1e-4, 1e-2, 1
   Record phi_res, divJ_L2_red, P_raw.

6. Temporarily split literature_passed into:
   demo_passed
   phi_solver_passed
   divJ_warning
   Do not let reconstructed divJ block the demo chain.

7. Customer demo uses demo mode only:
   --reload=1 --no-phi or phi as experimental overlay
   50 kW scaling
   Output A/B/E/J/JouleHeat.
```

---

## 9. Final Judgment

### Current blocker

```text
Not geometry;
not relax itself;
not Biot-Savart;
not power calibration;
but PhiImag correction discrete consistency / sign / divJ acceptance method.
```

### What can be shown now

```text
Can show demo chain:
coil -> A/B -> E/J -> JouleHeat -> 50 kW heat source
```

### What cannot be shown now

```text
Cannot show as:
OPHELIE literature mode already passed;
phi correction satisfies divJ continuity;
full French literature reproduction.
```

### Next-step core

```text
First fix PhiImag correction with manufactured solve tests so it truly lowers divJ;
then consider self-induction;
then thermal coupling or more complex geometry.
```
