# OPHELIE / TEAM7 High-σ Complex Edge-Flux Bottleneck Explanation and Next-Steps Plan

**Purpose**: For Cursor / Codex as basis for next-round code remediation, test design, and discussion.  
**Topic**: Why TEAM7 currently does not pass; whether complex edge-flux route is correct; which issues to address first; how to test and proceed.  
**Conclusion summary**: Current TEAM7 issues do not mean integral-form solver or edge-flux route failed, but that current implementation has not yet passed high-conductivity standard eddy-current benchmark. Next stage should prioritize cleaning normalization / restore mechanism, building high-σ analytic/semi-analytic benchmark, refactoring physical-scale Picard, then return to TEAM7 L3 validation.

---

## 0. Answer Core Questions First

User questions:

> What does "exposing complex edge-flux induced-current amplitude closure issue on high-conductivity TEAM7 aluminum plate" mean?  
> Does it mean the integral-form solver itself is defective and cannot solve TEAM7?  
> Or the solver is unsuitable for TEAM7?  
> Or only that current implementation is not yet resolved?

Clear answer:

**Not that integral-form solver inherently cannot solve TEAM7.**  
**Not that complex edge-flux / Biot–Savart integral form is inherently unsuitable for TEAM7.**  
More reasonable current judgment:

> TEAM7 exposed issues not yet closed in current complex edge-flux implementation, especially induced-current amplitude calibration in high-conductivity aluminum plate, normalization / restore, self-induction feedback, boundary/geometry discretization consistency, and TEAM7 validation definition / scope.

That is, **route is not wrong; implementation has not yet reached strict validation level for TEAM7 standard benchmark**.

---

## 1. Is Integral Form / Complex Edge-Flux Theoretically Suitable for TEAM7?

Yes.

TEAM7 is essentially a typical eddy-current problem:

\[
A_{\mathrm{coil}} \rightarrow E \rightarrow J \rightarrow A_{\mathrm{ind}}, B_{\mathrm{ind}}
\]

This matches current OPHELIE-like complex edge-flux route:

\[
A_{\mathrm{total}} = A_{\mathrm{coil}} + A_{\mathrm{ind}}
\]

\[
E = -\nabla \phi - i\omega A_{\mathrm{total}}
\]

\[
J = \sigma E
\]

\[
A_{\mathrm{ind}} = K[J]
\]

Therefore, integral form / Biot–Savart is naturally suitable for unbounded air domain, coil source field, conductor interior eddy currents, self-inductance feedback, and probe-line magnetic field comparison. French OPHELIE route is also integral form, so TEAM7 current non-pass cannot be directly interpreted as "integral form unsuitable".

More accurate statement:

> TEAM7 is a very suitable high-conductivity standard benchmark to test current complex edge-flux OPHELIE-like solver; current non-pass means implementation and validation definition / scope still need closure, not that solver route itself is rejected.

---

## 2. Why French Reduced Can Pass but TEAM7 Is Hard?

French reduced glass case conductivity is approximately:

```text
sigma_glass ~ O(10) S/m
```

TEAM7 aluminum plate conductivity:

```text
sigma_aluminum ~ 3.5e7 S/m
```

About 6 orders of magnitude difference.

This amplifies many issues not obvious in low-conductivity glass:

1. whether edge-flux RHS / pair conductance `C_ij` physical scale is fully consistent;
2. whether `phi` correction is sufficient to offset non-conservative part of induced electric field;
3. whether `EEdgeRecon` is strictly V/m;
4. whether `JEdgeRecon = sigma * EEdgeRecon` is strictly A/m²;
5. whether `A_ind` feedback must be self-consistent Picard;
6. whether normalization / restore pollutes feedback and amplitude judgment;
7. whether SPH discretization near boundary is robust enough at high `sigma`;
8. whether coil source, probe position, reference phase convention are strictly aligned.

Therefore, French reduced passing means main chain is effective in low-conductivity cold-crucible reduced case; but TEAM7 is stricter high-conductivity eddy-current stress test and cannot automatically pass.

---

## 3. Specific Bottlenecks Exposed by Current TEAM7

Current key phenomenon:

```text
J_sim / J_ref already too large;
B_sim / B_ref also too large;
(B/J)_sim not outrageous enough to explain all B error alone.
```

Therefore current main issue is not simple Biot probe amplification, but:

> absolute current amplitude closure issue on chain `A_coil -> phi -> E_edge -> J_edge`.

Must check:

```text
1. edge_drop physical unit strictly V;
2. phi unit is V;
3. EEdgeRecon strictly V/m;
4. JEdgeRecon = sigma * EEdgeRecon strictly A/m²;
5. whether pair conductance C_ij only used for residual or implicitly affects current amplitude;
6. whether normalization / restore distorts feedback test;
7. whether high-sigma Picard must use physical-scale state iteration;
8. whether TEAM7 geometry, coil source, probe and reference data fully aligned.
```

---

## 4. Wrong Conclusions Not to Draw

### 4.1 Cannot Say "Integral-Form Solver Cannot Solve TEAM7"

No evidence supports this. Integral form/Biot–Savart is common for eddy-current problems; TEAM7 is suitable to test it.

### 4.2 Cannot Say "Edge-Flux Method Unsuitable for TEAM7"

More accurate:

```text
current L2 one-way edge-flux insufficient to represent TEAM7 standard eddy-current coupled solution.
```

TEAM7 standard solution is not simply:

```text
A_coil -> solve once for J -> postprocess B_ind
```

Should at least:

```text
A_total = A_coil + A_ind
A_total -> solve phi/J
J -> A_ind
repeat until self-consistent
```

i.e. need L3 Picard or stronger coupled solve.

### 4.3 Cannot Immediately Apply Empirical Scaling to J or C_ij

For example:

```text
J *= 0.25
C_ij *= 0.25
RHS *= 0.25
```

This may make TEAM7 phase90 look better superficially, but breaks French reduced, uniform-field test or other analytic benchmarks.

Before changing operator scale, must first do:

```text
normalization / restore invariance
high-sigma analytic/semi-analytic benchmark
physical-scale Picard
J/Biot/probe unit audit
```

---

## 5. Most Reasonable Current Judgment

For three possibilities user raised, current ranking:

### Possibility 1: Integral-form solver defective, cannot solve TEAM7

**Not supported currently.**

### Possibility 2: Solver unsuitable for TEAM7

**Also not supported.** TEAM7 is instead a very suitable validation benchmark.

### Possibility 3: Current implementation has unresolved issues

**This is the most reasonable current judgment.**

Specific issues:

```text
complex edge-flux absolute J amplitude in high sigma aluminum plate, normalization/restore, self-induction feedback and TEAM7 validation definition / scope not yet closed.
```

---

## 6. Impact on Project Main Line

This is not bad news, but good news.

In French reduced case:

```text
sigma lower;
geometry simpler;
validation focus on total power, Q distribution and Picard main chain;
```

Many high-conductivity issues may be masked. TEAM7 amplifies issues and helps discover earlier:

```text
1. whether current amplitude reliable in high sigma limit;
2. whether A_ind feedback is real;
3. whether complex phasor closes;
4. whether J/B reference can align;
5. whether SPH edge-flux credible under standard eddy-current benchmark.
```

If TEAM7 eventually passes, EM solver credibility increases greatly.

---

# 7. Next-Steps Overall Route

Do not immediately "tune parameters to make TEAM7 look like it passes". Proceed as:

```text
1. clean normalization / restore, ensure physical output not polluted by numerical scaling;
2. build high-sigma analytic / semi-analytic benchmark;
3. refactor physical-scale Picard so feedback iterates at physical scale;
4. return to TEAM7 L3 Picard;
5. add Jey + Bz dual validation;
6. only then consider adjusting edge-flux operator / RHS / C_ij.
```

---

# 8. P0 Issues Cursor Must Address Immediately

## P0.1 normalization / restore invariance must come first

Current most important issue is whether normalization / restore affects physical results and feedback judgment.

### Current risk

If flow is:

```text
1. scale physical A/E/J fields down by scale;
2. solve phi/J;
3. finally restore physical fields;
4. use restored or pre-restore results for feedback judgment;
```

then feedback / Picard may be polluted by scale mechanism.

Especially already seen:

```text
feedback pre-restore changes greatly;
post-restore returns to one-way;
```

This means current feedback test cannot directly serve as physical solution.

### CLI to add

```bash
--ophelie-edge-flux-normalization-mode=off
--ophelie-edge-flux-normalization-mode=field-scale-restore
--ophelie-edge-flux-normalization-mode=solver-local
```

Where:

```text
field-scale-restore = current approximate implementation;
solver-local = target implementation, physical fields keep true dimensions, only normalize RHS / phi / operator temporarily inside solver.
```

### normalization sweep to add

Add test or TEAM7 CLI:

```bash
--team7-normalization-sweep=1
```

Sweep:

```text
safe_rhs_l2 = 1e4, 1e5, 1e6, 1e12
safe_rhs_max_abs = default / enlarged
```

Output CSV:

```text
safe_rhs_l2
safe_rhs_max_abs
normalization_mode
input_scale
J_imag_L2_vol_pre_restore
J_imag_L2_vol_post_restore
Bind_over_Bcoil_pre_restore
Bind_over_Bcoil_post_restore
phase90_RMS_pre_restore
phase90_RMS_post_restore
P_recon_pre_restore
P_recon_post_restore
restore_invariance_error_J
restore_invariance_error_B
restore_invariance_error_P
```

### Acceptance requirements

In one-way linear problem:

```text
post-restore J/B/P should be insensitive to safe_rhs_l2;
restore_invariance_error should be near numerical error;
```

If post-restore results change clearly across safe_rhs_l2, normalization still affects physical output and must be fixed first.

---

## P0.2 TEAM7 validation logic must stay strict

L2 one-way cannot masquerade as TEAM7 validation. Recommendation:

```text
L2 one-way:
  smoke_passed=1 OK;
  diagnostic_only=1;
  team7_validation_passed=0;

L3 Picard:
  can be validation candidate.
```

`team7_validation_passed` should at least include:

```cpp
report.team7_validation_passed =
    strict_candidate &&
    report.team7_phase0_source_passed &&
    report.team7_phase90_probe_passed &&
    picard_ok &&
    em_ok &&
    coil_path_audit_ok &&
    !report.diagnostic_only;
```

If `team7_validation_passed` logic does not include `phase90_probe_passed`, future misjudgment likely; must fix.

---

## P0.3 freeze phasor convention, no longer use a_sign to tune phase90

Continue freezing:

```text
imag_a_sign = +1
phase90 convention = minus-imag at probe layer
```

If user sets:

```text
imag_a_sign != +1
```

must:

```text
diagnostic_only=1
team7_validation_passed=0
```

Output warning:

```text
[team7][warning] imag_a_sign != +1 changes the solver equation.
[team7][warning] This is debug-only and cannot be used as TEAM7 validation.
```

phase90 sign should only be handled at probe comparison layer:

```bash
--team7-phase90-convention=minus-imag|plus-imag
```

Do not use solver equation `a_sign` to fit reference phase.

---

## P0.4 source_scale / j_post_scale / diagnostic_only definition / scope

Strict validation requires:

```text
source_scale = 1.0
j_post_scale = off or 1.0
imag_a_sign = +1
normalization audited
level = picard
```

If:

```text
source_scale != 1.0
j_post_scale = auto or != 1.0
imag_a_sign != +1
normalization_mode = field-scale-restore and invariance not passed
```

then must:

```text
diagnostic_only=1
team7_validation_passed=0
```

`j_post_scale=auto` only as debug, not validation.

---

# 9. P1: high-sigma analytic / semi-analytic benchmark

Before changing edge-flux operator, must use simpler high-sigma cases to judge whether issue is systematic scale error.

## P1.1 uniform E high-sigma box

Construct conductor box / slab:

\[
A = 0,\quad \phi = -E_0\cdot x
\]

Exact:

\[
E = E_0,\quad J=\sigma E_0,\quad Q=\frac{1}{2}\sigma |E_0|^2
\]

Test different:

```text
sigma = 16, 1e3, 1e5, 3.5e7 S/m
```

Check:

```text
E_recon / E_exact
J_recon / (sigma E_exact)
P_recon / P_exact
q_antisymmetry
edge_residual
```

If high-sigma uniform E passes, `J = sigma E` and Joule heat basically OK.

## P1.2 harmonic A high-sigma box

Set:

\[
\phi = 0,\quad A_{\mathrm{real}} = A_0
\]

Exact imaginary chain:

\[
E_{\mathrm{imag}} = -\omega A_0,\quad J_{\mathrm{imag}}=\sigma E_{\mathrm{imag}}
\]

Check:

```text
E_imag_recon / (-omega A0)
J_imag_recon / (sigma * -omega A0)
P_recon / P_exact
omega scaling
sigma scaling
```

This directly checks core TEAM7 chain `A_real -> J_imag`.

## P1.3 conducting cylinder skin-depth benchmark

Priority recommend adding this case because closest to induction heating:

```text
cylindrical conductor;
axial or circumferential harmonic magnetic/vector-potential excitation;
validate skin depth delta;
compare J(r), Q(r), total P;
frequency scaling at f = 50, 200, 300kHz etc.
```

Use:

```text
1. judge whether edge-flux amplitude correct under high-sigma + skin effect;
2. judge whether J depth distribution reasonable at 50 Hz / 200 Hz;
3. provide analytic/semi-analytic validation anchor beyond TEAM7.
```

---

# 10. P2: physical-scale Picard

Current Picard cannot only sweep relaxation; must first ensure iteration at physical scale.

## P2.1 Picard state must be physical scale

Each Picard iteration should be:

```text
A_total_physical^k = A_coil_physical + A_ind_physical^k
solve edge-flux using physical-scale active A
reconstruct J_physical^k
A_ind_physical^{k+1} = K[J_physical^k]
under-relax A_ind or A_total in physical scale
repeat
```

Should not use pre-restore intermediate fields for physical feedback judgment.

## P2.2 Prefer relaxing A_ind / A_total, not only J

Current fixed point:

```text
A_total -> solve J -> K[J] -> A_ind -> A_total
```

So more natural under-relaxation:

\[
A_{\mathrm{ind}}^{k+1} \leftarrow
(1-\alpha) A_{\mathrm{ind}}^{k}
+
\alpha K[J^k]
\]

or:

\[
A_{\mathrm{total}}^{k+1} =
A_{\mathrm{coil}} + A_{\mathrm{ind}}^{k+1}
\]

Not only relax `J`. J-relax can be kept as diagnostic, but production Picard should prioritize A-field relaxation.

## P2.3 Picard CSV must output

Each iteration output:

```text
iter
normalization_mode
normalization_scale
J_pre_restore
J_post_restore
J_rel_raw_physical
J_rel_relaxed_physical
Aind_rel_physical
P_recon
P_rel_physical
Aind_over_Acoil
Bind_over_Bcoil
phase90_RMS
source_model
edge_res_red_real
edge_res_red_imag
q_antisym_real
q_antisym_imag
max_J_real
max_J_imag
```

Acceptance cannot look only at `J_rel`; also:

```text
P_rel
Aind_rel
phase90 trend
Bind/Bcoil trend
finite fields
```

---

# 11. P3: TEAM7 L3 Picard validation candidate

After normalization definition / scope and physical-scale Picard are cleaned, run strict L3:

```bash
test_3d_ophelie_team7_complex_edge_flux \
  --team7-level=picard \
  --team7-validation-mode=strict \
  --team7-source-model=filament-racetrack \
  --team7-coil-source-scale=1.0 \
  --team7-phase90-convention=minus-imag \
  --ophelie-edge-flux-imag-a-sign=1 \
  --team7-frequency=50
```

First goal is not immediate pass, but clear failure reasons:

```text
team7_validation_passed=0
fail_reason:
  phase90 failed
  picard not converged
  Jey failed
  source failed
  normalization not strict
```

If after L3:

```text
Bind/Bcoil drops clearly from one-way ~4;
phase90 RMS drops;
Picard converges;
Jey close;
```

then missing feedback in one-way is main cause.

If after L3 still:

```text
Bind/Bcoil ≈ 4;
Jey still clearly too large;
phase90 still large deviation;
```

then enter P4 operator scale audit / boundary audit.

---

# 12. P4: edge-flux operator scale audit

Only if P0–P3 still show systematic current amplitude too large, consider modifying edge-flux operator scale.

Audit content:

## P4.1 `C_ij` dimension and consistency

Current form roughly:

```cpp
C_ij = harmonicMean(sigma_i, sigma_j)
       * pairwiseNegativeLaplaceWeight(dW_ij * Vol_j, ...)
```

Need output statistics:

```text
C_ij_min/max/mean
sum_j C_ij * |r_ij|^2
sigma_i * Vol_i
ratio = sum_j C_ij * |r_ij|^2 / (sigma_i * Vol_i)
```

Goal is to judge scale relation between pairwise Laplace conductance and physical admittance.

## P4.2 constant/linear potential consistency

Test:

\[
\phi=-E_0\cdot x,\quad A=0
\]

Check:

```text
edge_drop / |r|
E_recon
J_recon
P_recon
residual
```

## P4.3 harmonic A consistency

Test:

\[
\phi=0,\quad A=A_0
\]

Check:

```text
edge_drop = omega A · r
E_recon = -omega A
J_recon = sigma E
```

## P4.4 Boundary effects

TEAM7 is finite plate + hole geometry; boundary conditions very important. Must distinguish:

```text
interior particles
near external boundary
near hole boundary
top/bottom surface
```

Output partitioned J/E/Q/edge residual.  
If error concentrates near boundary, may not be global scale but boundary edge-flux / Neumann handling issue.

---

# 13. P5: Jey reference and 200 Hz

## P5.1 Jey reference usage rules

Reference value zero points cannot use relative error directly or ratio blows up. Recommendation:

```text
Jey_RMS_abs_all
Jey_RMS_abs_nonzero_ref
Jey_rel_RMS_nonzero_ref
Jey_sign_agreement
Jey_peak_value_sim/ref
Jey_peak_location_sim/ref
```

`reference=0` points separately:

```text
Jey_zero_ref_abs_error
```

## P5.2 Meaning of 200 Hz

200 Hz currently mainly means error scales linearly/quadratically with `f`, no obvious `omega^2` equation error. Later 200 Hz for:

```text
skin-depth trend
frequency response
phase90 scaling
Jey scaling
```

But priority lower than 50 Hz baseline scale and normalization/Picard closure.

---

# 14. Code Structure Remediation Recommendations

## 14.1 `team7_validation_passed` logic

Must ensure strict validation includes phase90:

```cpp
report.team7_validation_passed =
    strict_candidate &&
    report.team7_phase0_source_passed &&
    report.team7_phase90_probe_passed &&
    picard_ok &&
    em_ok &&
    coil_path_audit_ok &&
    !report.diagnostic_only;
```

## 14.2 `ComputeOphelieEdgeFluxPhiRhsFromASrcCK` naming

Class now reads active A from component, not necessarily `ASrc`. Recommend rename:

```text
ComputeOphelieEdgeFluxPhiRhsFromAComponentCK
```

If not renaming yet, at least comment:

```text
Historical name. The kernel actually uses component.active_a_field, not necessarily ASrc.
```

## 14.3 `measureOphelieEdgeFluxRhsNormalizationScale()` side effects

If function internally setup/finalizes RHS, not pure measure. Recommendation:

```text
1. rename to setupAndMeasureOphelieEdgeFluxRhsNormalizationScale;
or
2. truly make side-effect-free measure with temporary field or restore original RHS.
```

Otherwise call order easily pollutes RHS.

## 14.4 complex power metrics

In complex mode, power metrics should read:

```text
joule_heat_edge_recon_complex
```

Not imag-only field:

```text
joule_heat_edge_recon_imag
```

Log should also distinguish:

```text
P_recon_complex
P_recon_imag_only
P_recon_real_only
```

## 14.5 coil_source_scale default/strict logic

Can temporarily keep legacy smoke 0.754, but strict validation must:

```text
source_scale = 1.0
```

And any non-1.0 scale automatically:

```text
diagnostic_only=1
team7_validation_passed=0
```

## 14.6 output tag

All CSV/log must carry tag to avoid overwrite:

```text
team7_probe_<level>_<tag>.csv
team7_picard_<tag>.csv
team7_bz_<level>_<tag>.csv
team7_jey_<level>_<tag>.csv
team7_summary_<level>_<tag>.txt
```

---

# 15. Do Not Do

Next round do not:

```text
1. Do not directly multiply J by 0.25;
2. Do not directly multiply C_ij or RHS by empirical factors;
3. Do not write TEAM7-specific A-phi equation;
4. Do not use a_sign=-1 to fix phase90;
5. Do not let j_post_scale=auto participate in validation;
6. Do not do 50 kW power scaling for TEAM7;
7. Do not mark L2 one-way as team7_validation_passed=1;
8. Do not continue interpreting feedback/Picard physical trends before normalization/restore audited;
9. Do not pivot back to div-grad residual floor.
```

---

# 16. Direct Execution Summary for Cursor

Can send directly to Cursor:

```text
Current judgment:

TEAM7 non-pass does not mean integral / complex edge-flux route failed.
It means current implementation has not passed high-conductivity standard eddy-current benchmark.
Issues mainly exposed in absolute amplitude closure of A_coil -> phi -> EEdgeRecon -> JEdgeRecon,
and physical-scale consistency of normalization/restore and self-induction feedback.

Next steps: do not use empirical factors on J or C_ij. Proceed in order:

P0 normalization / restore:
1. Add --ophelie-edge-flux-normalization-mode=off|field-scale-restore|solver-local.
2. Add --team7-normalization-sweep=1, sweep safe_rhs_l2.
3. Output pre/post restore J/B/P/phase90, compute restore_invariance_error.
4. feedback/Picard must use physical-scale state; field-scale-restore un-audited only diagnostic_only.

P0 validation:
5. Fix team7_validation_passed, must include phase90_probe_passed.
6. L2 one-way always diagnostic_only, team7_validation_passed=0.
7. imag_a_sign=+1 frozen; phase90 only minus-imag at probe layer.
8. source_scale!=1, j_post_scale!=1/auto, a_sign!=+1 all must diagnostic_only=1.

P1 high-sigma benchmark:
9. Add uniform E high-sigma box.
10. Add harmonic A high-sigma box.
11. Add conducting cylinder skin-depth benchmark.
12. Check E/J/P absolute scale vs sigma, omega, A.

P2 physical-scale Picard:
13. Picard iteration uses physical A_total/J/A_ind.
14. Prioritize A_ind or A_total under-relaxation, not only J.
15. Each iteration CSV output J_rel_raw, J_rel_relaxed, Aind_rel, P_rel, Bind/Bcoil, phase90_RMS.

P3 TEAM7 L3:
16. Run L3 Picard with filament source_scale=1, phase90=minus-imag, imag_a_sign=+1.
17. Output clear failure reasons, do not demand immediate pass.

P4 operator audit:
18. If high-sigma benchmark and L3 still show J systematically too large, audit C_ij/RHS/edge reconstruction.
19. Do not multiply empirical factors just to get Bind/Bcoil≈1.
```

---

# 17. Summary

Current TEAM7 issue should be understood as:

> TEAM7 as high-conductivity standard eddy-current benchmark exposes that complex edge-flux implementation J amplitude, normalization/restore, physical-scale feedback and high `sigma` discretization consistency are not yet closed.

This is not route failure, but entering stricter validation stage.

Final goals remain:

```text
French reduced:
  validate cold-crucible power scale, Q distribution, Picard EM main line.

TEAM7:
  validate high sigma standard eddy-current benchmark B/J/frequency response.

high-sigma analytic benchmark:
  validate edge-flux operator and J/E/Q absolute scale.
```

All three together support a credible OPHELIE-like complex edge-flux electromagnetic solver.
