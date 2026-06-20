# OPHELIE / TEAM7 Round 3 Analysis and Next-Steps Cursor Execution Plan

Date: 2026-06-11 (execution progress update: 2026-06-12)  
Topic: TEAM7 L2/L3 complex edge-flux validation, high-σ aluminum plate bottleneck, P1a test definition / scope correction, follow-on development tasks

### Execution progress (2026-06-12)

| Item | Status | Record |
|----|------|------|
| P0 P1a gauge cancellation | done | `HIGH_SIGMA_EDGE_FLUX_P1A_LOG.md` v4 |
| P1 rotational-A | done | J_phi/J_edge≈0.91, no 4× |
| P2 hole_lateral + Jn/Jt | done | `TEAM7_VALIDATION_LOG.md` §P2 |
| P3 L1 strict definition / scope doc | done | phase0 total ≠ pure source reference |
| P4 A_ind/Lenz audit | done | box corr≈−0.80; TEAM7 phasor≈−0.59, lenz_pass=1 |
| F probe CSV skip | done | missing reference does not write empty file |
| P4.1 moment complete | done | aniso~0.06; boundary~0.13; ratio/inv(h³)≈4 |
| P4.5 hole edge-recon | done | e_edge_em 0.35 vs 0.05; fallback=0 |
| GPT Round 3 closure pack | done | `OPHELIE_GPT_DISCUSSION_TEAM7_ROUND3_CLOSURE.md` |
| P5 Picard relax sweep | frozen | |

---

## 0. This Round Conclusion Summary

This round Round 3 data very valuable, but must avoid misjudgment. Current most important judgment:

1. **normalization / restore already largely ruled out as accidental source of one-way amplitude deviation.**  
   After safe_rhs_l2 sweep from `1e4` to `1e12`, post-restore `J`, `Bind/Bcoil`, `phase90_RMS` basically unchanged; one-way physical output insensitive to current normalization threshold.

2. **solver-local normalization and field-scale-restore equivalent under one-way.**  
   Means `solver-local` can be cleaner default route later; `field-scale-restore` kept as legacy diagnostic.

3. **Biot probe is not separate amplification source.**  
   Current `J_sim/J_ref` already too large, `B_sim/B_ref` also too large; `B/J` mapping not outrageous enough to explain all error alone. Core still:
   ```text
   A_coil -> φ -> E_edge -> J_edge
   ```
   absolute amplitude closure on this chain.

4. **Round 3 P1a `harmonic_A_imag_phi_solve` conclusion needs correction.**  
   Current test uses spatially constant vector potential `A0`. Spatially constant `A0` is pure gauge, `curl(A0)=0`; isolated conductor φ equation should cancel it, edge-drop and J near zero. Therefore after full φ-solve `J≈0` is not failure, but reasonable gauge cancellation.  
   Cannot use this test to prove "φ re-solve path broken" or cancel φ solve accordingly.

5. **Cannot directly decide to modify `C_ij` scale.**  
   `sum(C*r^2)/(σ*Vol) ≈ 1.5e8` looks large but not a clearly defined dimensionless error. For `dp≈3 mm`, kernel moment scale similar to `1/h^3` can reach `1e7–1e8`. Cannot directly conclude `C_ij` wrong from this ratio, cannot multiply `C_ij` or `J` by empirical factor.

6. **TEAM7 still has not validation passed.**  
   Should treat TEAM7 as high-σ standard eddy-current benchmark exposing issues, not use empirical scaling to pass phase90. Next steps must separate validation of gauge, non-conservative A, boundary normal current, A_ind/Lenz-law, Picard and source/probe.

---

## 1. Already Ruled Out or Basically Clear

### 1.1 normalization / restore not accidental source of one-way deviation

Round 3 P0 sweep main results:

```text
safe_rhs_l2: 1e4 -> 1e12
input_scale: 0.00413 -> 1
J_post:      63993 basically unchanged
Bind/B:      3.987 basically unchanged
phase90:     11.81 basically unchanged
```

Means:

```text
one-way post-restore physical output is invariant to safe_rhs_l2.
```

Therefore `Bind/B≈4` not accidentally caused by specific normalization threshold.

But does not mean feedback / Picard physical interpretation already clean. Picard still must ensure:

```text
A_ind, A_total, J, E, Q all participate in iteration in physical scale.
```

Later TEAM7 strict validation should prioritize:

```text
--ophelie-edge-flux-normalization-mode=solver-local
```

`field-scale-restore` can be kept legacy diagnostic.

### 1.2 solver-local and field-scale-restore equivalent under one-way

Round 3 P2β shows:

```text
field-scale-restore: Bind/B=3.987, J=63993, phase90=11.81
solver-local:        Bind/B=3.987, J=63993, phase90=11.81
```

One-way path same physical output from both normalization modes. Later recommend:

```text
TEAM7 / complex edge-flux strict validation:
    use solver-local normalization
```

Keep:

```text
field-scale-restore:
    legacy / diagnostic only
```

### 1.3 Biot probe not independent amplification source

Earlier J/B separation:

```text
J_sim/J_ref already too large
B_sim/B_ref also too large
(B/J)_sim roughly comparable to (B/J)_ref
```

Cannot attribute phase90 / Bz deviation to:

```text
Biot–Savart probe map alone.
```

Main issue still:

```text
A_coil -> φ -> E_edge -> J_edge
```

recovered current amplitude.

### 1.4 frequency scaling no obvious ω² error

200 Hz results roughly:

```text
J ~ 4 × J_50Hz
P ~ 16 × P_50Hz
```

Matches frequency-domain linear expectation. Not low-level error like extra ω factor.

---

## 2. Round 3 Key Definition / Scope Correction: constant A is gauge-cancellation test

### 2.1 Issue with current P1a test

Current P1a `harmonic_A_imag_phi_solve` setup roughly:

```cpp
A_coil_real = A0;  // spatially constant
phi_imag = 0;
expected E = -omega * A0;
```

Then saw:

```text
edge-only: passes
phi_solve + restore: J nearly zero, marked failed
```

Recorded as:

```text
edge-only correct, phi_solve path failed.
```

This judgment needs correction.

### 2.2 Why is J≈0 after constant A full φ-solve reasonable?

Spatially constant vector potential satisfies:

\[
\nabla \times A_0 = 0.
\]

Pure gauge / conservative drive. For isolated conductor, φ equation can generate:

\[
\phi \approx -\omega A_0 \cdot x,
\]

so pairwise edge-drop:

\[
e_{ij}
=
(\phi_j-\phi_i)
+
\omega A_0 \cdot (x_j-x_i)
\approx 0.
\]

Thus:

\[
J \approx 0.
\]

After full φ-solve `J≈0` is not φ solve wrong, likely **gauge cancellation correct**.

### 2.3 How to modify P1a test name and acceptance?

Rename:

```text
harmonic_A_imag_phi_solve
```

to:

```text
constant_A_gauge_cancellation
```

New expected:

```text
edge-only:
    E = -ω A0
    J = σ E
    edge reconstruction diagnostic only, not isolated-conductor physical solve.

full φ-solve:
    φ cancels constant A
    edge_drop_after_phi ≈ 0
    J_after_phi ≈ 0
    gauge_cancellation_passed = true
```

Add output fields:

```text
gauge_cancellation_passed
phi_linear_fit_error
edge_drop_after_phi_l2
edge_drop_after_phi_linf
J_after_phi_norm
J_after_phi_over_edge_only
```

Recommend acceptance:

```text
J_after_phi / J_edge_only << 1
edge_drop_after_phi_l2 small
```

### 2.4 Cannot therefore cancel φ solve

Must clarify:

```text
Do not switch production to edge-only.
```

TEAM7 coil field not spatially constant A; has nonuniform spatial structure and nonzero curl. φ solve still key for conductor interior current conservation and normal-current boundary.

---

## 3. Truly Unclosed Issues

### 3.1 TEAM7 L1 coil source still not strictly closed

Although filament source connected and runs with `source_scale=1.0`, L1 still clear differences:

```text
phase0 source RMS ≈ 0.527
peak_sim/ref ≈ 9.85 / 7.81 mT
peak_x_sim/ref ≈ 198 / 126 mm
```

Two issues:

1. source field spatial shape not fully aligned;
2. using phase0 total reference for source validation may not be strict.

TEAM7 table phase0 is total field at one time of standard problem, not necessarily pure coil-only or f=0 source-only field. Strict L1 should use:

```text
f=0 source-only reference
or no-plate / σ=0 coil-only trusted reference
or trusted FEM coil-only reference
```

Therefore before L1 source fully closed, should not use phase90 as hard validation gate.

### 3.2 high σ boundary / normal-current treatment may be major error source

P4 partition results:

```text
interior:       Bind/B ≈ 3.79, e_edge_em_mis ≈ 0.068
top_surface:    Bind/B ≈ 3.87, e_edge_em_mis ≈ 0.198
bottom_surface: Bind/B ≈ 4.99, e_edge_em_mis ≈ 0.239
lateral:        Bind/B ≈ 4.19, e_edge_em_mis ≈ 0.202
```

Means:

1. interior also too large, not pure boundary issue;
2. surface mismatch ~3× interior, boundary treatment clearly has issues.

TEAM7 perforated thick plate: top/bottom, outer side, hole wall boundaries. Current:

```text
hole_lateral_proxy n=0
```

Hole wall diagnostic not captured. Must-fix diagnostic gap.

Next must add:

```text
n · J audit
J_normal / J_tangential
top / bottom / outer_lateral / hole_lateral / interior partition
edge residual by partition
edge-drop / EMF mismatch by partition
```

If significant normal current at conductor boundary, φ correction and boundary condition not truly closed.

### 3.3 `C_ij` ratio cannot directly prove conductance scale error

P4.1:

```text
sum(C*r^2)/(σ*Vol) median ≈ 1.5e8
```

Cannot directly treat as error; not clearly defined dimensionless quantity. For `h≈dp≈0.003 m`:

\[
1/h^3 \approx 3.7\times 10^7.
\]

With kernel moment, neighbor count and coefficients, `1e8` scale not necessarily abnormal.

More important: if multiply all pair `C_ij` by same constant, φ equation LHS and RHS may scale together, φ solution may not change; edge reconstruction LS weight overall scaling may cancel. Therefore:

```text
Do not patch Cij merely based on conductance_ratio≈1.5e8.
```

Next need moment-consistency audit:

```text
M_i = Σ C_ij r_ij ⊗ r_ij
trace(M_i)
eigenvalues(M_i)
anisotropy
interior/boundary split
comparison to σ, volume, kernel moment
```

Not directly multiply `C_ij` by empirical factor.

### 3.4 Picard currently cannot be main hope to fix Bind/B

Current L3 sweep:

```text
one-way:   Bind/B≈3.99, phase90≈11.81
relax=0.05 Bind/B≈4.32, phase90≈7.11
relax=0.15 Bind/B≈6.52, phase90≈14.66
relax=0.35 Bind/B≈42.6, phase90≈112
```

Current Picard does not naturally push Bind/B down, may worsen. phase90 RMS drops at `relax=0.05` but Bind/B still ~4, not stable validation.

Therefore wording:

```text
L3 Picard is required for final TEAM7 validation,
but current L3 is only a diagnostic.
Do not assume Picard will fix Bind/B≈4.
```

Later Picard should rerun after:

```text
source/probe convention clearer
boundary audit complete
A_ind / Lenz-law sign audit complete
rotational-A benchmark passed
solver-local physical-scale state used
```

---

## 4. Decisions on Round 3 Key Questions

### Q1: Cij ratio≈1.5e8 means must change Cij?

**Insufficient to conclude. Cannot change directly.**

Must first define correct moment / dimensional audit, then judge whether generic pair conductance needs correction.

### Q2: P1a means TEAM7 should not re-solve φ, only edge-only?

**No. Absolutely not.**

P1a constant-A full solve should be interpreted as gauge cancellation. TEAM7 production must keep φ solve.

### Q3: solver-local should replace field-scale-restore?

**Yes.**

Recommend later TEAM7 strict and complex edge-flux default:

```text
--ophelie-edge-flux-normalization-mode=solver-local
```

field-scale-restore kept legacy diagnostic.

### Q4: Priority boundary Neumann or global Cij scale?

**Priority boundary audit and non-conservative A benchmark, not direct Cij fix.**

Order:

```text
1. fix hole wall/boundary partitions;
2. add Jn/Jt normal-current audit;
3. add rotational-A / uniform-B benchmark;
4. do A_ind / Lenz-law audit;
5. then judge whether generic operator needs modification.
```

### Q5: e_edge_em_mismatch surface ~3× interior means boundary not closed?

**Yes, at least boundary edge reconstruction / continuum EMF consistency poorer.**

Output by region:

```text
interior
top
bottom
outer side
hole wall
```

And dp refinement.

### Q6: hole_lateral_proxy n=0 what to do?

**Must fix.**

TEAM7 hole is key geometry feature; cannot lack hole wall diagnostic. Recommend from TEAM7 explicit geometry parameters or level-set / STL distance:

```text
distance_to_hole_wall
distance_to_outer_side
distance_to_top
distance_to_bottom
```

Then partition:

```text
hole_lateral_skin
outer_lateral_skin
top_surface
bottom_surface
interior
```

### Q7: Continue Picard relaxation sweep?

**Not recommended large relax sweep now.**

Next Picard should change to:

```text
A_ind / A_total under-relaxation
```

Not only J relaxation. Fixed point essentially:

```text
A_total -> solve J -> K[J] -> A_ind -> A_total
```

Therefore relaxing A_ind or A_total more natural.

### Q8: When can team7_validation_passed=1?

Still far. At least need:

```text
source_scale=1
filament or validated source
imag_a_sign=+1
phase90 convention fixed
level=picard
solver-local normalization
no diagnostic_only
coil source phase0 passed
phase90 passed
Jey passed or optional passed
Picard converged
q_antisym / edge_res / phi_eq_res passed
boundary normal-current acceptable
```

L2 one-way never `team7_validation_passed=1`.

### Q9: TEAM7 acceptance use physical-scale solver-local?

**Yes.**

Recommend definition:

```text
TEAM7 strict validation uses solver-local physical-scale output.
field-scale-restore is diagnostic only.
```

### Q10: Round 3 data enough to lock Cij/RHS change?

**Insufficient.**

Can only lock:

```text
TEAM7 J amplitude is too large;
normalization threshold is not the accidental source;
source/edge/boundary/Picard still not fully closed.
```

Cannot lock:

```text
Cij/RHS must be modified.
```

### Q11: Jey ~10× and Bind/B≈4 same root cause?

Likely same source, both J-driven. Different ratios mean probe, spatial distribution, reference points, boundary factors superimposed. Jey needs finer statistics:

```text
absolute RMS
relative RMS over nonzero reference points
sign agreement
peak location
median ratio over significant points
```

### Q12: Next code priority?

Recommend:

```text
1. fix P1a benchmark physical expected;
2. add rotational-A / uniform-B nonconservative benchmark;
3. complete TEAM7 boundary and hole wall Jn/Jt audit;
4. do A_ind / Lenz-law sign audit;
5. return to physical-scale Picard;
6. last judge whether generic edge-flux operator needs change.
```

---

## 5. Concrete Next-Step Development Plan

### P0: fix P1a test definition / scope and documentation

#### P0.1 rename constant-A test

Change:

```text
harmonic_A_imag_phi_solve
```

to:

```text
constant_A_gauge_cancellation
```

#### P0.2 modify expected result

```text
edge-only:
    expected E = -ω A0
    diagnostic only

full φ solve:
    expected edge_drop ≈ 0
    expected J ≈ 0
    gauge_cancellation_passed = true
```

#### P0.3 add output

```text
gauge_cancellation_passed
phi_linear_fit_error
edge_drop_after_phi_l2
edge_drop_after_phi_linf
J_after_phi_norm
J_after_phi_over_edge_only
```

#### P0.4 modify log conclusion

No longer write:

```text
phi_solve + restore failed
```

Should write:

```text
constant-A full solve performs gauge cancellation.
```

### P1: add non-conservative A benchmark

Currently missing A field test not fully gauge-cancelable by φ. Recommend add:

```text
rotational_A_uniform_B
```

Set:

\[
A = \frac{1}{2} B_0 \times r
\]

For example:

\[
A=(-\frac12 B_0 y,\frac12 B_0 x,0),
\]

then:

\[
\nabla \times A = B_0 \hat{z}.
\]

Non pure-gauge A field; φ cannot cancel it entirely. For cylinder or regular conductor, one-way should produce circumferential current.

Recommend new or extend:

```text
test_3d_ophelie_high_sigma_edge_flux_scaling
case_name=rotational_A_uniform_B
```

Output:

```text
E_theta_sim / E_theta_exact
J_theta_sim / (σ E_theta_exact)
radial_E_leak
normal_J_boundary
P_recon / P_exact
e_edge_em_mismatch
q_antisym
edge_res_red
```

### P2: complete TEAM7 boundary audit

#### P2.1 fix hole_lateral_proxy

Current:

```text
hole_lateral_proxy n=0
```

Must-fix diagnostic gap.

From TEAM7 geometry parameters or level-set / STL distance:

```text
distance_to_hole_wall
distance_to_outer_side
distance_to_top
distance_to_bottom
```

Partition:

```text
interior
top_surface
bottom_surface
outer_lateral
hole_lateral
corner / edge optional
```

#### P2.2 partitioned output

Each partition output:

```text
n_particles
J_normal_L2
J_tangential_L2
Jn_over_Jt
E_normal_L2
E_tangential_L2
edge_residual_L2
e_edge_em_mismatch
Bind/Bcoil
J_magnitude_mean
J_magnitude_max
```

Goal: judge whether high-σ error boundary-dominated, interior-dominated, or both.

### P3: L1 source / probe strict audit

Filament source not fully closed; do not look only at peak fit.

Output:

```text
phase0_RMS_abs
phase0_RMS_rel
peak_sim/ref
peak_x_sim/ref
best_fit_scale_L2
best_fit_scale_peak
correlation
```

If no f=0 source-only reference, log must clarify:

```text
phase0 total reference is not a pure source-only reference.
```

Later preferably prepare:

```text
f=0 source-only reference
or no-plate / σ=0 trusted FEM reference
```

### P4: A_ind / Lenz-law sign audit

Current Picard did not lower Bind/B; possible reasons:

```text
Picard not converged
A_ind feedback sign/phase wrong
J amplitude too large
boundary normal-current not closed
```

Recommend simple Lenz-law audit.

In `rotational_A_uniform_B` or cylinder benchmark:

```text
imposed A -> solve J -> compute A_ind/B_ind
check whether induced field opposes imposed time-varying magnetic drive
```

TEAM7 at least output:

```text
corr(B_ind, B_coil)
corr(B_total, B_coil)
A_ind_real/imag sign relative to coil
B_ind_real/imag sign relative to coil
B_total amplitude compared with B_coil
```

If A_ind feedback enhances instead of shields, Picard cannot fix Bind/B.

### P5: rerun physical-scale Picard

After P1/P2/P4, not large relaxation sweep now.

Add relaxation target:

```text
--self-induction-relax-target=J
--self-induction-relax-target=AInd
--self-induction-relax-target=ATotal
```

Prioritize:

```text
AInd or ATotal under-relaxation
```

Not only relax J.

Picard CSV output:

```text
iter
J_rel_raw
J_rel_relaxed
Aind_rel
Atotal_rel
P_recon
P_rel
Aind_over_Acoil
Bind_over_Bcoil
Bz_phase0_RMS
Bz_phase90_RMS
max_J_real
max_J_imag
edge_res_red_real
edge_res_red_imag
q_antisym_real
q_antisym_imag
normalization_mode
relax_target
```

---

## 6. Code Structure Review and Recommendations

### 6.1 What is done well

Structure already clearer than Round 2:

```text
electromagnetic_ophelie_team7_validation.h:
    pass report, diagnostic_only, phase90 convention

electromagnetic_ophelie_team7_probe.h:
    reference loading and probe definitions

electromagnetic_ophelie_team7_coil_path_source.h:
    filament racetrack source

electromagnetic_ophelie_edge_flux_operator_audit.h:
    conductance and partition audit

test_3d_ophelie_high_sigma_edge_flux_scaling:
    high-σ unit benchmark
```

Correct direction; keep and extend.

### 6.2 Code structure issues to fix

#### Issue 1: P1a test naming / expected misleading

`harmonic_A_imag_phi_solve` name and expected misleading. Change to:

```text
constant_A_gauge_cancellation
```

And modify pass/fail judgment.

#### Issue 2: `ComputeOphelieEdgeFluxPhiRhsFromASrcCK` name inaccurate

Now reads active A from component, not necessarily ASrc. Recommend rename:

```text
ComputeOphelieEdgeFluxPhiRhsFromAComponentCK
```

If not renaming class, at least comment:

```text
Historical name. The kernel uses component.active_a_field, not necessarily ASrc.
```

#### Issue 3: conductance ratio easily misleading

Do not only output:

```text
sum(C*r²)/(σ*Vol)
```

Change to moment audit:

```text
M_i = Σ C_ij r_ij ⊗ r_ij
trace(M_i)
eigenvalues(M_i)
anisotropy
interior / boundary split
```

#### Issue 4: P4 probe CSV empty

`team7_probe_one-way_p4_audit.csv` empty; P4 audit output path or call order issue. Must fix.

#### Issue 5: source_scale default still 0.754

Can keep smoke default, but strict validation must force:

```text
source_scale=1.0
```

And all:

```text
source_scale != 1.0
```

must:

```text
diagnostic_only=1
team7_validation_passed=0
```

---

## 7. Explicitly Forbidden Directions

Next do not:

```text
J × 0.25 empirical scaling
C_ij × empirical scaling
RHS × empirical scaling
a_sign = -1 to fix phase90
TEAM7 50 kW power scaling
TEAM7-only A–φ special equation
edge-only production
removing φ solve
declaring TEAM7 validation passed from L2 one-way
```

If operator modification later proved needed, must be:

```text
generic edge-flux operator correction
```

not TEAM7-only hack.

---

## 8. Direct Execution Summary for Cursor

```text
Current Round 3 key conclusions need correction and refinement:

1. P0 normalization sweep proved one-way post-restore insensitive to safe_rhs_l2.
   solver-local can be later default normalization.
   This does not mean TEAM7 closed, only normalization not accidental source of current 4×.

2. P1a harmonic_A_imag_phi_solve cannot be interpreted as "φ solve failed".
   constant A is pure gauge; full φ solve should cancel it; J≈0 reasonable.
   Rename case constant_A_gauge_cancellation, modify expected:
     edge-only: E=-ωA, diagnostic only
     full φ solve: J_after_phi≈0, edge_drop≈0, gauge cancellation passed

3. Do not cancel φ solve because of P1a, do not switch to edge-only production.
   TEAM7 A_coil is nonuniform, has curl; φ solve must be kept.

4. Add true non-conservative A benchmark:
   rotational_A_uniform_B:
     A=(-0.5*B0*y, 0.5*B0*x, 0)
     curl A = B0 z
   check E_theta, J_theta, P_recon, e_edge_em_mismatch, Jn/Jt.
   Better than constant A for judging high-σ edge-flux amplitude.

5. Continue TEAM7 boundary audit:
   hole_lateral_proxy cannot be n=0.
   partition with TEAM7 geometry or level-set/STL distance:
     top, bottom, outer_lateral, hole_lateral, interior
   per region output:
     J_normal, J_tangential, Jn/Jt, edge_residual, e_edge_em_mismatch, Bind/Bcoil.

6. Do not directly change Cij scale.
   conductance_ratio≈1.5e8 not clear dimensionless error.
   global Cij scaling may not change φ solution and edge reconstruction.
   if changing later, must be based on moment consistency / analytic benchmark, not tuning Bind/B to 1.

7. Picard currently diagnostic, not fix means.
   before boundary / rotational-A / Aind sign audit complete, no large relaxation sweep.
   later Picard prioritize A_ind or A_total under-relaxation, not only J.

8. P4 probe CSV empty, fix output.
   P4 audit should not produce empty probe file.

9. Continue keeping:
   imag_a_sign=+1
   phase90=minus-imag
   source_scale=1 for strict
   no TEAM7 50kW scaling
   no j_post_scale validation
   no TEAM7-only A-phi equation
```

---

## 9. Final Judgment

Round 3 not failure, but further narrowed issue scope. Most important avoid misjudgment:

```text
constant A after phi_solve J≈0 is not φ solve broken;
it is gauge cancellation.
```

Therefore cannot decide now:

```text
change Cij scale
cancel φ solve
accept TEAM7 discrete 4×
```

Correct next steps:

```text
1. fix P1a benchmark physical expected;
2. add rotational-A / uniform-B nonconservative benchmark;
3. complete TEAM7 boundary and hole wall Jn/Jt audit;
4. do A_ind / Lenz-law sign audit;
5. return to physical-scale Picard;
6. last judge whether generic edge-flux operator needs change.
```

Current main line not "tune numbers to pass TEAM7", but separately validate:

```text
gauge
boundary
self-induction sign
nonconservative A drive
source/probe reference
```

Only after these close can judge whether TEAM7 4× deviation from operator scale, boundary treatment, or source/reference alignment.
