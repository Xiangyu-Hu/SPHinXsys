# OPHELIE / TEAM7 P5 No-Flux Boundary Route Review and Next-Steps Cursor Execution Plan

**Date**: 2026-06  
**Topic**: TEAM7 P5 no-flux boundary diagnostic follow-on remediation direction  
**Scope**: SPHinXsys/SYCL OPHELIE complex edge-flux solver; TEAM7 high-σ validation; French reduced glass delivery side-line decision

---

## 0. Conclusion Summary

This round's P5 no-flux boundary work already gives a fairly clear boundary:

1. **no-flux physical target still correct**. For conductor-air boundary, correct constraint is

   \[
   n\cdot J=0,
   \]

   i.e. conduction current cannot leave conductor boundary. It is not \(B=0\) or \(A=0\). Magnetic field and vector potential still exist outside conductor, implementing coil--plate or coil--glass nonlocal magnetic coupling via Biot--Savart / integral coupling.

2. **P5.1--P5.5 failure does not negate no-flux or level-set / static-confinement route**. These tests are mainly E/J post-processing layer, edge-reconstruction layer, RHS-only layer patches; they did not truly complete \(\phi\)-LHS / kernel-support / static-confinement-level boundary operator correction.

3. **Currently not recommended to continue full P5.6 \(\phi\)-LHS static confinement directly on TEAM7 hole case**. Large effort, too many variables, easily mixes source/reference, geometry, boundary and operator issues.

4. **level-set is potential formal scheme, but last**. Current stage first do lighter, quickly validatable fixes: tangent-LS distance normalization fix, partition fix, EMF mismatch diagnostic split, L1 source/reference/probe audit, \(\phi\)-Neumann/LHS benchmark on box/slab.

5. **production default still keep boundary mode = none**. P5 no-flux modes temporarily all kept as opt-in diagnostic / regression, not strict validation, not production heat source.

---

## 1. How to Understand Current P5 Results

### 1.1 P5.0: normal connection successful

`NormalDirection` / `SignedDistance` can already be written via SPHinXsys body shape / reload flow and used for audit. Direction correct, should keep.

Continue noting:

- top/bottom/exterior normal angle error small, basically usable;
- hole_lateral normal mean angle error large but not catastrophic;
- corner / edge region normals naturally unstable, cannot mix into ordinary surface statistics;
- must distinguish `true_interior`, `boundary_shell_all`, `corner_or_edge`, avoid misclassifying boundary shell as interior.

### 1.2 P5.1/P5.2: post-hoc projection and tangent-LS did not solve TEAM7

`project-normal` can push boundary normal current near zero:

\[
E_{proj}=E-n(n\cdot E),\qquad J_{proj}=\sigma E_{proj}.
\]

But it did not re-solve \(\phi\) or fix \(\nabla\phi\) boundary operator. Therefore comparing `projected E_edge` with uncorrected

\[
E_{em}=-\nabla\phi-\omega A
\]

`e_edge_em_mismatch` worsening is not surprising. This only means:

> post-hoc projection cannot be production fix; cannot only change E/J post-processing without fixing \(\phi\) / gradient boundary consistency.

Cannot mean no-flux physical boundary condition is wrong.

### 1.3 P5.3: consistency audit runs, but partitions need fix

Current `interior/all` no-flux metrics may be polluted by boundary shell / interior mislabel. Especially if `interior` normal angle near 85°, those particles should not enter normal boundary statistics.

Must change partitions to:

```text
true_interior
boundary_shell_all
corner_or_edge
top_surface
bottom_surface
outer_lateral
hole_lateral
```

`true_interior` should satisfy signed-distance far enough from boundary, not participating in `Jn/Jt` normal boundary evaluation.

### 1.4 P5.4: no-flux ghost edge has diagnostic value, not production fix

Current ghost-edge essentially adds constraint to E-reconstruction LS matrix:

\[
E\cdot n = 0.
\]

Slightly improves surface `Jn/Jt`, but not `e_edge_em`, `Bind/Bcoil`, `P_recon`. Means it touches reconstruction-layer no-flux, not \(\phi\)-solve / \(\nabla\phi\) / operator consistency.

Recommendation:

```text
keep opt-in; default production = none.
```

Do not promote to strict validation or production.

### 1.5 P5.5: missing moment / phi-RHS ghost did not truly solve issue

#### missing moment

Current missing moment more like LS regularization than true level-set static confinement. If implementation similar:

```cpp
m_acc += moment_weight * (I - n n^T)
```

and `moment_weight` from `max_neighbor_conductance * ghost_distance^2`, magnitude likely too small to truly complete boundary missing support. Its failure does not mean level-set correction invalid.

Recommend naming/log fix:

```text
no-flux-missing-moment currently only experimental LS regularization, not true level-set missing-support correction.
```

#### phi-RHS ghost

RHS-only ghost is not complete Neumann boundary condition. EM no-flux ghost condition should be edge-drop zero:

\[
e_{ig}=(\phi_g-\phi_i)+\omega A\cdot r_{ig}=0.
\]

If only RHS added, no LHS counterpart, not complete \(\phi\)-boundary closure; result equals baseline unsurprising.

Recommend log fix:

```text
RHS-only diagnostic; incomplete without LHS counterpart.
```

---

## 2. Key Code Points to Re-Audit

### 2.1 tangent-LS distance normalization may have issue

Current tangent LS if using:

```cpp
directional_e = edge_drop / distance;
```

where `distance = |r_ij|`, needs careful re-audit.

For tangent-plane LS, if boundary particle only solves tangential components:

\[
e_{ij}\approx E_t\cdot r_{t,ij},
\]

more natural normalization:

```cpp
directional_e = edge_drop / distance_t;
```

where:

\[
r_t = r - n(n\cdot r),\qquad distance_t=|r_t|.
\]

Recommend Cursor A/B comparison:

```text
tangent-ls-distance3d:      edge_drop / |r|
tangent-ls-distance-tangent: edge_drop / |r_t|
```

Output comparison:

```text
e_edge_em_mismatch
P_recon
Jn/Jt
Bind/Bcoil
phase90 RMS
Jey error
```

Before fixing this, should not fully judge tangent-LS route failed.

### 2.2 EMF mismatch diagnostic needs normal/tangent component split

Current projected/tangent modes `e_edge_em_mismatch` may be unfair because `E_edge` projected or tangentialized but `E_em=-gradPhi-omegaA` has no corresponding projection or Neumann correction.

Recommend add comparisons:

```text
raw_Eedge_vs_raw_Eem
projected_Eedge_vs_raw_Eem
projected_Eedge_vs_projected_Eem
tangent_Eedge_vs_tangent_Eem
normal_component_mismatch
tangential_component_mismatch
```

Core questions:

- is normal mismatch eliminated by no-flux projection;
- does tangential mismatch remain;
- do TEAM7 main B/Jey errors come from normal leak or tangential amplitude.

### 2.3 boundary mode enum recommend layering

Now one `OphelieEdgeReconBoundaryMode` mixes:

```text
post-hoc diagnostic: project-normal, tangent-ls
reconstruction ghost: no-flux-ghost-edge
LS regularization: no-flux-missing-moment
phi RHS diagnostic: no-flux-phi-rhs-ghost
combo: no-flux-full
```

Recommend later split into three:

```text
edge_recon_diagnostic_mode
edge_recon_production_mode
phi_boundary_mode
```

Short term if not splitting, must clarify each mode's action layer in log.

---

## 3. Continue P5.6 Next Steps?

### 3.1 Not recommended direct full P5.6 on TEAM7

Not recommended to continue large-scale \(\phi\)-LHS static confinement directly on TEAM7 hole case. Reasons:

1. L1 source/reference/probe not yet strictly closed;
2. TEAM7 hole geometry complex, boundary, hole, corner, source/probe issues mixed;
3. P5.5 missing moment not yet true level-set correction;
4. direct LHS change high risk, high regression cost;
5. end-of-June time window limited.

### 3.2 Recommended P5.6-lite

If continuing TEAM7 effort, only narrow `P5.6-lite`:

```text
1. fix tangent-ls distance_t comparison;
2. fix true_interior / boundary_shell / corner partitions;
3. split normal/tangent EMF mismatch diagnostic;
4. parallel L1 source/reference/probe audit;
5. validate φ-LHS / Neumann / static-confinement ideas on box/slab benchmark;
6. if simple benchmark shows clear benefit, return to TEAM7.
```

level-set scheme last: only when above light fixes and benchmarks point to clear boundary operator issue, invest in formal level-set missing-support correction.

---

## 4. Execution Priority Before End of June

### P5-fix: existing implementation and diagnostic fixes, highest priority

#### P5-fix-1: tangent-LS distance_t A/B test

Implement or temporary compare:

```text
tangent-ls-distance3d
tangent-ls-distance-tangent
```

Compare TEAM7 metrics:

```text
e_edge_em_mismatch by region
Jn/Jt by region
P_recon
Bind/Bcoil
phase90 RMS
Jey error
```

#### P5-fix-2: partition fix

Add:

```text
true_interior
boundary_shell_all
corner_or_edge
top_surface
bottom_surface
outer_lateral
hole_lateral
```

Requirements:

```text
true_interior does not participate in no-flux normal angle / Jn/Jt boundary evaluation;
hole_lateral n>0 as regression;
corner_or_edge separate statistics, not mixed into top/bottom/hole.
```

#### P5-fix-3: normal/tangent EMF mismatch

Add output:

```text
Eedge_normal
Eem_normal
Eedge_tangent_norm
Eem_tangent_norm
normal_mismatch
tangential_mismatch
```

For raw/projected/tangent/ghost-edge separately.

---

### P6a: L1 source/reference/probe audit, high priority

Current TEAM7 phase0 reference should not be used directly as pure source-only reference. Must separately advance:

```text
1. clarify phase0 reference physical meaning;
2. if possible, prepare f=0 source-only reference;
3. or prepare no-plate / sigma=0 FEM reference;
4. or at least high-precision filament Biot self-reference for source-only audit;
5. check probe line coordinates, x origin, units mm/m, coil/plate/hole relative position.
```

Output:

```text
Bz_source_RMS_abs
Bz_source_RMS_rel
peak_sim/ref
peak_x_sim/ref
correlation
best_fit_scale_L2
best_fit_scale_peak
```

If L1 not closed, should not change generic operator based on L2 phase90 alone.

---

### P6b: φ-LHS / Neumann / static-confinement only on simple benchmark

Do not go directly to TEAM7 hole. First box/slab.

#### benchmark 1: constant-A gauge cancellation

Ensure new boundary LHS does not break:

\[
A=A_0,\qquad J\approx0.
\]

#### benchmark 2: linear φ / Neumann boundary

Construct exact no-flux case, check boundary residual drops.

#### benchmark 3: rotational-A slab

Closer to TEAM7 than box. Output:

```text
J magnitude
Jn/Jt
e_edge_em
P_recon
Lenz sign
```

If simple cases show φ-LHS / Neumann / static-confinement clear improvement, return to TEAM7.

---

## 5. Level-Set Route Last

User has SPHinXsys static confinement / level-set experience; level-set is potential most formal, cleanest route. But large development effort; should not be first next-step means.

Recommended order:

```text
1. tangent-ls fix;
2. partition and EMF diagnostic fix;
3. L1 source/reference/probe;
4. box/slab φ-Neumann benchmark;
5. if effective, then level-set missing-support correction.
```

Future true level-set correction should consider:

```text
missing gradient moment
missing edge LS moment
boundary no-flux closure
Neumann φ boundary correction
no real conductor ghost generated
external points do not participate in Biot JdV or Joule heat
```

---

## 6. Answers to Cursor Key Questions

### Q1: Continue P5.6 or downgrade TEAM7?

Short term continue TEAM7, but only P5.6-lite; no direct large TEAM7 φ-LHS change. If no clear improvement by end of June, July downgrade TEAM7 to benchmark side line, main line French reduced glass delivery.

### Q2: Is no-flux BC physically correct?

Yes. Correct variable is conduction current:

\[
n\cdot J=0.
\]

Not \(B=0\) or \(A=0\). Current failure is low-order discrete patch, not physical BC itself.

### Q3: If continuing coding, priority A/B/C/D?

Priority:

```text
D: L1 source/reference/probe audit
+
B/C simple benchmark validation
```

Not recommended direct A-class pairwise Laplacian mirror ghost on TEAM7.

### Q4: P5.4 ghost edge keep or rollback?

Keep opt-in, default production = none. Can be regression/diagnostic, not production or strict validation.

### Q5: interior shell partition mislabel fix first?

Yes, must fix first. Otherwise `interior/all` `Jn/Jt` and normal metrics mislead.

### Q6: hole wall e_edge_em≈0.35 L2 hard blocker?

Not L2 smoke hard blocker; but for future strict validation must be hard metric.

### Q7: Only φ-LHS worth doing?

No. E/RHS patches basically tried; if continuing boundary operator route, next must touch φ-LHS / Neumann / static confinement, but L1 source/reference/probe audit equally important, even higher priority.

### Q8: L2 add no-flux/hole hard gate?

No. L2 one-way still diagnostic_only. Soft metrics sufficient:

```text
hole_lateral e_edge_em
top/bottom e_edge_em
Jn/Jt
Bind/Bcoil
phase90 RMS
```

---

## 7. Direct Execution Summary for Cursor

```text
P5 current conclusion:
no-flux physical target n·J=0 still correct, but P5.1-P5.5 only prove E/J post-processing, E-reconstruction ghost, RHS-only ghost cannot solve TEAM7. Do not thereby negate no-flux or level-set/static confinement.

Next steps: do not push large P5.6 directly on TEAM7, do not do formal level-set immediately. level-set last.

Proceed in order:

1. Fix tangent-ls:
   current tangent-ls may use edge_drop/|r|.
   Add edge_drop/|r_t| version, A/B compare:
     tangent-ls-distance3d
     tangent-ls-distance-tangent
   Output e_edge_em, Jn/Jt, P_recon, Bind/Bcoil, phase90, Jey.

2. Fix partitions:
   Add true_interior, boundary_shell_all, corner_or_edge.
   true_interior does not participate in no-flux normal/JnJt boundary evaluation.
   hole_lateral n=0 should warning/fail diagnostic.

3. Fix EMF mismatch diagnostic:
   Add normal/tangential component mismatch:
     raw_Eedge_vs_raw_Eem
     projected_Eedge_vs_projected_Eem
     tangent_Eedge_vs_tangent_Eem
   output by region.

4. L1 source/reference/probe audit:
   current phase0 total reference cannot be pure source-only reference.
   prepare f=0/no-plate/source-only reference or high-precision filament Biot self-reference.
   check probe x coordinates, origin, units, coil/plate/hole relative position.

5. φ-LHS / Neumann / static-confinement only on box/slab benchmark:
   do not go directly to TEAM7 hole.
   must keep constant-A gauge cancellation;
   test linear φ no-flux;
   test rotational-A slab.

6. P5.4 no-flux ghost edge keep opt-in, default none.
   no-flux-missing-moment mark as experimental LS regularization, not true level-set.
   phi-RHS ghost mark as RHS-only incomplete diagnostic.

Still forbidden:
   Cij×empirical factor;
   J×constant;
   a_sign=-1;
   edge-only production;
   large Picard relaxation sweep;
   TEAM7 50kW scaling.

If no clear improvement from P5-fix + L1 source/reference + box/slab φ-LHS benchmark by end of June, July downgrade TEAM7 to benchmark side line, main line French reduced glass delivery.
```

---

## 8. Final Recommendation

Current production default keep:

```text
--ophelie-edge-recon-boundary-mode=none
```

P5 no-flux modes keep opt-in. Next focus not stacking patches, but:

```text
1. fix current tangent-LS / partition / EMF diagnostic;
2. close L1 source/reference/probe;
3. validate φ-LHS / Neumann / static-confinement on simple box/slab benchmark;
4. level-set last, invest when prior evidence clear.
```

Avoid blind large changes on TEAM7 hole complex geometry while preserving chance to break through TEAM7 before end of June.
