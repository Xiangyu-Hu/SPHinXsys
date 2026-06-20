# OPHELIE / TEAM7 Boundary No-Flux Treatment, Normal Acquisition, and Cursor Next-Steps Implementation Plan

Date: 2026-06-14  
Scope: SPHinXsysSYCL / OPHELIE complex edge-flux / TEAM7 validation line  
Purpose: Move TEAM7 boundary issues from coarse "insufficient boundary-particle neighbors" suspicion to "geometry-normal-based no-flux boundary operator correction and diagnostics".

---

## 1. Unify Current Conclusions First

After this round's discussion, boundary issue wording needs correction. Do not focus on "boundary-particle cutoff radius has insufficient neighbors". In particle methods, boundary particles lacking exterior same-type neighbors is expected; SPHinXsys fluid, heat transfer, periodic boundary, static confinement etc. do not solve this by "boundary particles having full neighbors", but by wall/dummy/contact/mirror/level-set missing integral etc. to complete boundary operators.

For current OPHELIE / TEAM7 EM problem, what really needs checking and improvement:

```text
whether current edge-flux / φ / E-J reconstruction boundary operators correctly satisfy conductor no-flux condition,
i.e. n·J = 0, not simply counting whether boundary-particle neighbor count matches interior.
```

Here no-flux analogizes adiabatic heat conduction boundary, but correct variable is conduction current, not magnetic field.

Adiabatic heat conduction boundary:

\[
q_n=-k\nabla T\cdot n=0.
\]

Conductor EM boundary:

\[
n\cdot J=0.
\]

For isotropic conductor, if

\[
J=\sigma E,
\]

then equivalent to

\[
n\cdot E=0.
\]

Where total electric field is

\[
E=-\nabla\phi-i\omega A.
\]

In current real/imag component form, imag chain can write similarly

\[
E_{imag}=-\nabla\phi_{imag}-\omega A_{real}.
\]

Therefore boundary condition equivalent to

\[
n\cdot(\nabla\phi+\omega A)=0
\]

or in edge-flux discrete form, boundary ghost edge electromotive drop satisfies

\[
e_{ig}=0.
\]

Emphasize:

```text
not B=0;
not A=0;
not "magnetic field adiabatic";
but conductive current no-flux: n·J=0.
```

Magnetic field and vector potential exist outside conductor. Magnetic coupling between coil and plate/glass is implemented via A/B in air/vacuum. OPHELIE / Biot-Savart integral form uses this; no explicit air-domain boundary needed.

---

## 2. Main Gaps in Current EM Boundary Treatment

From current TEAM7 edge-flux implementation, conductor interior φ/E/J mainly via inner relation and pairwise edge-flux. For missing boundary support, no explicit boundary operator completion like SPHinXsys fluid wall/dummy, periodic mirror, static confinement / level-set missing integral.

Current rough status:

```text
1. conductor body inner relation established;
2. edge-flux residual uses conductor interior pairs;
3. φ solve via interior pairwise conservation;
4. E_edge from interior neighbor edge-drop least-squares reconstruction;
5. no conducting ghost/dummy outside boundary participating in φ/E/J;
6. no explicit level-set missing support correction;
7. no explicit n·J=0 constraint on boundary E/J reconstruction.
```

This may cause near TEAM7 boundary:

```text
1. edge-reconstructed E inconsistent with continuum EMF E=-gradφ-ωA;
2. boundary Jn/Jt too large;
3. hole_lateral / top / bottom / outer_lateral e_edge_em_mismatch higher than interior;
4. boundary LS reconstruction no fallback, but normal/tangential components may be mixed;
5. φ equation stable at residual level, but post-processed E/J does not satisfy no-flux boundary closure.
```

Round 3 / P4.5 already saw:

```text
hole_lateral:
  e_edge_em_mismatch ≈ 0.35
  fallback_frac = 0
  j_sigmaE_mis = 0
  Jn/Jt ≈ 0.51

interior:
  e_edge_em_mismatch ≈ 0.05
```

This means:

```text
LS ill-conditioning not main cause;
J=σE storage relation itself OK;
independent issue is boundary E_edge vs E_em=-gradφ-ωA consistency, and n·J no-flux closure.
```

---

## 3. Whether Geometry Normals Are Needed?

Boundary normals needed, but not recommended to manually provide per-particle normals. Acquire by priority:

```text
1. prioritize SPHinXsys existing body shape / level-set / surface normal variables;
2. if not directly exposed, derive normal from level-set signed distance;
3. if TEAM7 geometry analytically generated, use analytic geometry normals for first prototype;
4. external normal file only as last resort.
```

Cursor should first check whether TEAM7 body / particle generation / level-set already has normal fields, e.g.:

```text
NormalDirection
SurfaceNormal
normal_direction_
level-set normal
shape normal
signed distance gradient
```

Do not hardcode normal logic in projection code. Recommend generic interface first:

```cpp
bool getBoundaryNormal(size_t i, Vecd& normal);
```

or TEAM7-specific first version:

```cpp
bool getTeam7BoundaryNormal(size_t i, Vecd& normal, Team7BoundaryRegion& region);
```

Return meaning:

```text
true  -> particle i is boundary particle, normal valid;
false -> particle i not boundary particle, no no-flux boundary correction.
```

Interface internal priority:

```text
1. SPHinXsys existing surface normal / level-set normal;
2. TEAM7 analytic normal fallback;
3. unsupported -> false.
```

---

## 4. TEAM7 First-Version Analytic Normal Fallback

TEAM7 is thick plate + hole; first diagnostic can use analytic normals without waiting for formal level-set.

### 4.1 Top / bottom surface

If z is thickness direction:

```text
top_surface:    n = (0, 0, +1)
bottom_surface: n = (0, 0, -1)
```

### 4.2 Outer lateral surface

If outer contour approximately rectangular, judge by nearest edge:

```text
x = xmin -> n = (-1, 0, 0)
x = xmax -> n = (+1, 0, 0)
y = ymin -> n = (0, -1, 0)
y = ymax -> n = (0, +1, 0)
```

If corners exist, mark corner separately; do not mix into single side-wall statistics.

### 4.3 Hole lateral wall

If hole center \(c\), particle position \(x_i\), hole wall outward normal can first take in-plane radial:

\[
n = \frac{x_i-c}{|x_i-c|}.
\]

If multiple holes, use nearest hole center.

Must ensure:

```text
hole_lateral n_particles > 0;
hole normal points outside conductor / inside hole;
top/bottom/hole/outer/corner not overlapping.
```

---

## 5. Priority Plan: EM No-Current-Flux Diagnostic by Heat-Transfer No-Flux Analogy

Current priority is not dummy wall or Cij/J empirical scaling, but:

```text
using existing surface normal / level-set normal,
first do diagnostic-only no-flux E/J boundary reconstruction,
see whether TEAM7 boundary e_edge_em_mismatch, Jn/Jt, Bz/Jey/Bind improve.
```

Priority implement three boundary modes:

```bash
--ophelie-edge-recon-boundary-mode=none
--ophelie-edge-recon-boundary-mode=project-normal
--ophelie-edge-recon-boundary-mode=tangent-ls
```

If effective, later advance:

```bash
--ophelie-edge-recon-boundary-mode=no-flux-ghost-edge
--ophelie-edge-recon-boundary-mode=levelset-missing-moment
```

But latter two not first stage.

---

## 6. Mode 1: project-normal diagnostic

Lowest-cost no-flux test. For boundary particles, use existing normal \(n\), project edge-reconstructed E onto tangent plane:

\[
E_t=E-n(n\cdot E).
\]

Then:

\[
J_t=\sigma E_t.
\]

Forces boundary particles:

\[
n\cdot J_t=0.
\]

Note: first version diagnostic only, does not overwrite production fields.

Recommend output and compare raw vs projected:

```text
E_raw, J_raw
E_projected, J_projected
Jn/Jt raw vs projected
P_recon raw vs projected
B_ind raw vs projected
Bz phase90 RMS raw vs projected
Jey error raw vs projected
Bind/Bcoil raw vs projected
e_edge_em_mismatch raw vs projected by region
```

Criteria:

```text
if after project-normal:
  Jn/Jt clearly drops;
  boundary e_edge_em_mismatch clearly drops;
  Bz/Jey/Bind clearly improve;
then no-flux boundary is important direction for current TEAM7 error.
```

Limitation:

```text
only fixes reconstruction post-processing, not φ solve itself;
if φ residual / RHS boundary closure itself has issues, projection may not suffice.
```

---

## 7. Mode 2: tangent-LS diagnostic

Compared to project-normal, tangent-LS more physical. For boundary particles, no longer reconstruct E in 3D, reconstruct in tangent basis \(t_1,t_2\):

\[
E = E_1t_1 + E_2t_2.
\]

Where \(t_1,t_2\) from surface normal \(n\), satisfying:

\[
t_1\cdot n=0,
\quad
 t_2\cdot n=0.
\]

Edge-drop equation:

\[
e_{ij}\approx E\cdot r_{ij}
\]

becomes:

\[
e_{ij}\approx E_1(t_1\cdot r_{ij}) + E_2(t_2\cdot r_{ij}).
\]

Boundary particle E/J naturally satisfies:

\[
n\cdot E=0,
\quad
n\cdot J=0.
\]

Still diagnostic-only; output raw / project-normal / tangent-ls three-way comparison.

Criteria:

```text
if tangent-LS more stable than project-normal and clearly improves e_edge_em/Jey/Bz,
then boundary-aware edge reconstruction should become generic operator direction.
```

---

## 8. Later Mode 3: no-flux ghost edge

If project-normal / tangent-LS effective, can continue no-flux ghost edge. Analogizes adiabatic dummy in heat transfer: take \(T_g=T_i\) so \(T_i-T_g=0\), thus q=0.

EM edge-flux no-flux ghost is not simply \(\phi_g=\phi_i\), but boundary edge-drop:

\[
e_{ig}=0.
\]

Imag chain current edge-drop similar:

\[
e^{imag}_{ig}=(\phi^{imag}_g-\phi^{imag}_i)+\omega A^{real}_{ig}\cdot r_{ig}.
\]

no-flux ghost should satisfy:

\[
\phi^{imag}_g=\phi^{imag}_i-\omega A^{real}_{ig}\cdot r_{ig}.
\]

Real chain if defined:

\[
e^{real}_{ig}=(\phi^{real}_g-\phi^{real}_i)-\omega A^{imag}_{ig}\cdot r_{ig},
\]

then:

\[
\phi^{real}_g=\phi^{real}_i+\omega A^{imag}_{ig}\cdot r_{ig}.
\]

Important limits of this ghost edge:

```text
1. not real conductor degrees of freedom;
2. does not participate in Biot J dV;
3. does not produce Joule heat;
4. only for reconstruction / boundary closure;
5. if entering φ solve operator later, needs separate justification.
```

---

## 9. Later Mode 4: level-set missing-support / static confinement correction

If boundary no-flux modes effective, formal scheme should prioritize level-set/static-confinement style, not real dummy conductor.

Reason: outside conductor is not conductor; cannot add conducting dummy particles or they wrongly participate in φ/J/Biot/Joule heat.

Level-set correction goal:

```text
only complete geometry-missing operator integral / moment / boundary closure,
no external real conducting degrees of freedom introduced.
```

### 9.1 Level-set v1: missing edge moment correction

For boundary edge LS matrix:

\[
M_i^{edge}=\sum_j w_{ij}\hat r_{ij}\otimes\hat r_{ij},
\]

use level-set to complete missing support region geometry moment:

\[
M_i^{edge,corr}=M_i^{particle}+M_i^{missing}.
\]

First version can only complete matrix conditioning, not construct external real b term.

### 9.2 Level-set v2: no-flux boundary closure

Further from:

\[
n\cdot J=0
\]

or

\[
e_{ig}=0
\]

construct missing edge contribution. Closer to EM static confinement.

Larger development effort, but consistent with existing static confinement experience; may become formal scheme.

---

## 10. What Quantities to Test, Not Just Neighbor Count

Next diagnostic focus not `n_neighbors`. Neighbor count can be output additionally, but not conclusion core. Really check boundary operator results.

### 10.1 Gradient consistency

Manufactured fields:

\[
\phi=x,
\quad
\phi=y,
\quad
\phi=z.
\]

Expect:

\[
\nabla x=e_x,
\quad
\nabla y=e_y,
\quad
\nabla z=e_z.
\]

Output by region:

```text
grad_x_error
grad_y_error
grad_z_error
normal_component_error
tangential_component_error
```

### 10.2 Div-grad / Laplace consistency

Test:

```text
φ = constant
φ = linear field
```

Output:

```text
Lphi_const_L2
Lphi_linear_L2
boundary_Lphi_linear
interior_Lphi_linear
```

### 10.3 Edge reconstruction consistency

Construct:

\[
e_{ij}=E_0\cdot(x_j-x_i).
\]

Test:

```text
E0 = ex
E0 = ey
E0 = ez
E0 = boundary tangent direction
E0 = boundary normal direction
```

Output by region:

```text
E_recon_error
E_normal_error
E_tangential_error
E_angle_error
```

### 10.4 Neumann no-flux consistency

For real TEAM7 solution:

```text
Jn_L2
Jt_L2
Jn_over_Jt
Jn_signed_mean
Jn_abs_mean
```

For manufactured normal/tangent drive:

```text
boundary tangent drive -> Jn should be small
boundary normal drive -> φ/no-flux should cancel normal flux
```

### 10.5 constant-A gauge cancellation

\[
A=A_0,
\quad
\nabla\times A=0.
\]

After full φ solve:

```text
edge_drop_after_phi ≈ 0
J_after_phi ≈ 0
```

### 10.6 rotational-A non-conservative benchmark

\[
A=\frac12 B_0\times r.
\]

Check:

```text
Jθ
P
Lenz sign
Jn/Jt
boundary mode sensitivity
```

### 10.7 TEAM7 raw vs no-flux corrected comparison

Each TEAM7 one-way output:

```text
raw
project-normal
tangent-ls
```

Compare:

```text
phase90 RMS
Jey error
Bind/Bcoil
P_recon
e_edge_em_mismatch by region
Jn/Jt by region
```

---

## 11. Cursor Next-Steps Execution Order

Recommend next round in order:

### P5.0: find/connect normals

```text
1. check whether SPHinXsys / TEAM7 body already has surface normal / level-set normal;
2. if yes, connect to OPHELIE TEAM7 audit;
3. if no, use TEAM7 analytic normal fallback first;
4. implement getBoundaryNormal(i, normal) interface;
5. output normal audit CSV.
```

normal audit CSV at least:

```text
particle_id
region
x,y,z
normal_x,normal_y,normal_z
boundary_flag
signed_distance_or_region_distance
```

### P5.1: project-normal diagnostic

```text
1. add CLI: --ophelie-edge-recon-boundary-mode=project-normal;
2. for boundary particles execute E <- E - n(n·E);
3. J <- σE;
4. output raw vs projected metrics;
5. do not overwrite production field.
```

### P5.2: tangent-LS diagnostic

```text
1. add CLI: --ophelie-edge-recon-boundary-mode=tangent-ls;
2. construct t1/t2 from normal;
3. 2D tangent LS reconstruct E on boundary particles;
4. output raw / project-normal / tangent-ls comparison;
5. do not overwrite production field.
```

### P5.3: boundary operator consistency tests

```text
1. grad φ consistency;
2. edge reconstruction manufactured E;
3. no-flux normal/tangent drive;
4. constant-A gauge;
5. rotational-A benchmark;
6. by-region output.
```

### P5.4: decision

If project-normal / tangent-ls clearly improve:

```text
continue no-flux ghost edge / level-set missing moment;
consider upgrading boundary-aware reconstruction to generic operator.
```

If no improvement:

```text
continue L1 source/reference/probe audit;
do not change Cij/J yet;
prepare to switch French reduced main line in July.
```

---

## 12. Still Forbidden Directions

Before P5 results, do not:

```text
1. Cij × empirical factor;
2. J × 0.25;
3. RHS × empirical factor;
4. a_sign=-1 fix phase90;
5. edge-only production;
6. TEAM7 50 kW power scaling;
7. large Picard relaxation sweep;
8. add dummy conductor to real Biot/Joule heat.
```

---

## 13. TEAM7 Strategy Before End of June

User decision: continue TEAM7 through end of June if possible; if no substantive breakthrough by July, main line switches to French literature / reduced glass bath full development.

Recommend TEAM7 before end of June focus only three things:

```text
1. no-flux boundary correction: project-normal / tangent-ls;
2. level-set normal / missing support prototype;
3. L1 source/reference/probe coordinate audit.
```

If P5 shows no-flux correction can significantly lower `e_edge_em_mismatch` or improve Bz/Jey, continue TEAM7. Otherwise from July downgrade TEAM7 to high-conductivity benchmark side line, main line back to French reduced glass EM + Joule heat + thermal/stirring + σ(T) weak coupling.

---

## 14. Direct Summary for Cursor

```text
Next steps follow heat-transfer adiabatic boundary analogy, but EM no-flux means n·J=0, not B=0/A=0.

Need boundary normal n, but user does not need to manually provide per particle.
Prioritize SPHinXsys existing surface normal / level-set normal; if not connectable yet, TEAM7 analytic geometry normal fallback first.

Implement getBoundaryNormal(i,n) interface first, output normal audit CSV.

Add diagnostic-only boundary modes:
  --ophelie-edge-recon-boundary-mode=none|project-normal|tangent-ls

project-normal:
  E <- E - n(n·E)
  J <- σE

tangent-ls:
  construct t1/t2 from n, reconstruct E only in tangent plane.

Output raw/projected/tangent-ls comparison:
  Jn/Jt by region
  e_edge_em_mismatch by region
  Bind/Bcoil
  phase90 RMS
  Jey error
  P_recon

Also run boundary operator consistency tests:
  grad φ=x,y,z
  edge reconstruction manufactured E
  no-flux normal/tangent drive
  constant-A gauge cancellation
  rotational-A uniform-B

If no-flux modes clearly improve, advance no-flux ghost edge and level-set missing moment correction.

Forbidden: Cij empirical scaling, J empirical scaling, a_sign=-1, edge-only production, TEAM7 50kW scaling, large Picard sweeps.
```

---

## 15. Final Judgment

This no-flux boundary route is the most reasonable, lowest-risk, most SPHinXsys-consistent direction for current TEAM7 effort. Directly corresponds to adiabatic heat-transfer boundary idea, but variable is conductor current:

\[
n\cdot J=0.
\]

Does not cut off external magnetic field or break Biot-Savart nonlocal magnetic coupling. First stage only needs normals and diagnostic reconstruction, not immediately real dummy conductor. If first stage effective, further use familiar level-set/static confinement to complete missing edge moment and no-flux boundary closure.
