# OPHELIE-like / Biot-Savart SYCL Branch Code Review and Next-Step Test Plan

## 1. Overall Assessment

The current code is no longer a starter skeleton but a largely formed SPHinXsys/SYCL OPHELIE-like particle-integral electromagnetic module. Structure overall matches SPHinXsys SYCL/CK style:

- Uses `LocalDynamics` / `Interaction<Inner<>>` to encapsulate dynamics operators;
- Uses `StateDynamics<MainExecutionPolicy, ...>` and `InteractionDynamicsCK<MainExecutionPolicy, ...>` for scheduling;
- Device kernels access variables via `UpdateKernel` / `InteractKernel` and `DelegatedData(ex_policy)`;
- Uses `registerStateVariable` and `addVariableToWrite` for field management;
- Relaxation uses SYCL CK pipeline, not legacy `RelaxationStepInner`;
- Test case wired into CMake via `extra_sources_3d`;
- Supports reload, VTP output, CLI parameters, power normalization, phi solver, self-induction experimental path.

Therefore: **development form is basically correct and proceeds within the SPHinXsys SYCL framework.**

However, it is not yet recommended to treat `test_3d_ophelie_team7` as formal TEAM7 validation. Several key operator-level issues need fixing/testing first:

1. Biot-Savart `B` cross-product direction may be reversed;
2. TEAM7 coil `J0` equivalent current density formula may use volume instead of current cross-section area;
3. Self-induction complex phase handling should not be trusted yet;
4. Power scaling scales only `JouleHeat`, not `E/J` synchronously, causing inconsistent output fields;
5. Phi solver has PCG/GMRES/Jacobi but lacks fixed operator-level regression.

Conclusion: **TEAM7-like geometry smoke tests can start, but add several lightweight unit/operator tests before relying on TEAM7-like results.**

---

## 2. File-Level Review

## 2.1 `electromagnetic_ophelie.h`

Aggregate entry point including:

- CLI
- Biot-Savart
- device sync
- field names
- observables
- parameters
- Laplace
- phi
- phi GMRES/PCG/Jacobi
- postprocess
- self-induction
- register fields
- source
- TEAM7 geometry

Assessment: reasonable. If code grows, split into public umbrella header and internal headers; current structure can be preserved.

---

## 2.2 `electromagnetic_ophelie_parameters.h`

Contains default geometry, physical parameters, phi solver parameters, self-induction parameters.

Strengths:

- Centralized parameters;
- Clear `omega()`;
- Supports `coil_j0_override_`;
- Supports solver kind and self-induction switch.

Needs modification / notes:

### Issue A: TEAM7 override `coilAnnulusCrossSection()` actually returns volume

`OphelieTeam7Geometry::coilAnnulusCrossSection()` currently returns:

```cpp
Pi * (Rout^2 - Rin^2) * height
```

This is annular cylinder volume, not cross-section area perpendicular to azimuthal current.

If `J_src` is volumetric current density `A/m^2`, for azimuthal current the equivalent cross-section area should be:

```cpp
(Rout - Rin) * (2 * coil_half_height)
```

Or finer weighting for the actual annular section; rectangular cross-section area is sufficient for first version.

Current TEAM7:

```cpp
params.coil_j0_override_ = N * I / (geom.coilAnnulusCrossSection() + TinyReal);
```

Dimensions are wrong. Although later `JouleHeat` is normalized to target power (total `Q` may be unaffected), absolute `A/B/E/J`, raw power, and analytic/literature comparison will be wrong.

Recommended fix:

```cpp
Real coilCurrentCrossSectionArea() const
{
    return (coil_outer_radius_ - coil_inner_radius_) * (2.0 * coil_half_height_);
}
```

Then:

```cpp
params.coil_j0_override_ = N * I / (geom.coilCurrentCrossSectionArea() + TinyReal);
```

### Issue B: Default `phi_gauge_penalty_ = 1.0`

Helps numerical stability but changes pure Neumann phi problem. Keep it, but log output and run penalty sweep in tests:

```text
lambda_phi = 0.1, 1, 10
```

Check whether `JouleHeat` is stable.

---

## 2.3 `electromagnetic_ophelie_field_names.h`

Clear field organization, already distinguishing:

- coil source field;
- coil-induced `ACoil/BCoil`;
- glass self-induced `AInd/BInd`;
- working total `ASrc/BSrc`;
- `E/J/Q`;
- `PhiImag` solver workspace.

Assessment: reasonable.

Recommendations:

- If power scaling scales only `JouleHeat`, add `JouleHeatRaw` and `JouleHeatScaled` to avoid output confusion.
- If scaled heat goes to thermal, `JouleHeat` as scaled is fine, but logs must state `E/J` are raw-current values.

---

## 2.4 `electromagnetic_ophelie_register_fields.h`

Matches SPHinXsys field registration and VTP output. First-round fields are sufficient.

Recommend adding output for:

```text
ACoilReal
BCoilReal
ASrcReal
BSrcReal
EImag
JImag
JouleHeat
Sigma
PhiImag
GradPhiImag
```

Already largely implemented.

---

## 2.5 `electromagnetic_ophelie_source.h/.hpp`

`InitializeOphelieCoilSourceCK` is correct SYCL CK style:

```cpp
e_theta = (-dy/r, dx/r, 0)
JSrcReal = J0 * e_theta
JSrcImag = 0
```

Assessment: Stage 1 complete.

Needs testing:

- At `(x > center_x, y=center_y)`, `JSrcReal` should point `+y`;
- At `(x=center_x, y>center_y)`, `JSrcReal` should point `-x`;
- When `I0` or `J0` doubles, `JSrcReal` doubles.

---

## 2.6 `electromagnetic_ophelie_biot_savart.h/.hpp`

Functionality already exceeds first-version plan, including:

- coil → glass `ACoil/BCoil`;
- glass self-induced `AInd/BInd`;
- combine `A_src = A_coil + A_ind`.

Core file.

### Critical Issue: `B` cross product direction may be reversed

Current code:

```cpp
b_sum += coeff_ * r.cross(jv) * inv_r3;
```

Where:

```cpp
r = x_observer - x_source
jv = J_source * Vol_source
```

Standard Biot-Savart for volumetric current:

```cpp
dB = mu0/(4*pi) * (J dV) cross r / |r|^3
```

i.e.:

```cpp
jv.cross(r)
```

not:

```cpp
r.cross(jv)
```

Current form reverses `B` globally. `max_BSrc > 0` acceptance will not catch it (norm unchanged), but ParaView field direction will be wrong.

Fix immediately:

```cpp
b_sum += coeff_ * jv.cross(r) * inv_r3;
```

Apply same fix to self-induced B.

### Tests to Add

Minimal direction test:

- coil source at `(R,0,0)`;
- `J` points `+y`;
- observer at origin;
- `r = observer - source = (-R,0,0)`;
- `J cross r = +z`.

So center `Bz` should be positive. Current `r.cross(J)` gives negative `z`.

### Self-Induction Phase Issue

Self-induced uses `j_imag_` but writes:

```cpp
a_ind_real_[i] = a_sum;
a_ind_imag_[i] = 0;
```

With frequency phasors, Biot-Savart from `JImag` should land in imaginary component:

```cpp
a_ind_imag_[i] = a_sum_from_j_imag;
a_ind_real_[i] = 0;
```

Otherwise self-induced A phase is wrongly placed in the real channel, changing subsequent `-i omega A` phase relation.

Recommendation: **before fix, do not treat `--self-induction` as trustworthy physics; use only as code-path smoke test.**

---

## 2.7 `electromagnetic_ophelie_postprocess.h/.hpp`

Level 0 formulas correct:

```cpp
EReal =  omega * ASrcImag
EImag = -omega * ASrcReal
J = sigma * E
Q = 0.5 * (JReal dot EReal + JImag dot EImag)
```

With `J = sigma E`, equivalent to:

```cpp
Q = 0.5 * sigma * |E|^2
```

Assessment: Stage 3 complete.

### Note: Power Normalization Scales Only `JouleHeat`

Current:

```cpp
JouleHeat *= scale
```

But `EImag/JImag` are not scaled. VTP then has `JouleHeat != 0.5 * J dot E`.

Acceptable if “only Q goes to heat equation,” but confusing for client display and diagnostics.

Pick one:

**Option A: Keep current behavior, rename/document**

```text
JouleHeat = scaled heat source
E/J = raw fields from chosen I0
```

Add output:

```text
JouleHeatRaw
JouleHeatScaled
PowerScale
```

**Option B: Scale E/J as well**

If all outputs should be self-consistent:

```cpp
field_scale = sqrt(power_scale)
EReal *= field_scale
EImag *= field_scale
JReal *= field_scale
JImag *= field_scale
JouleHeat *= power_scale
```

Equivalent to adjusting coil current amplitude to target power.

For client visualization, recommend Option B.

---

## 2.8 `electromagnetic_ophelie_laplace.h`

Implements scalar variable-coefficient pairwise Laplace:

```text
OpheliePairwiseLaplaceCK
OpheliePairwiseLaplaceDiagonalCK
```

Uses harmonic mean and pair regularization; matches phi correction needs.

Assessment: Stage 4 LHS operator already available.

Needs testing:

- constant field Laplace near 0;
- diagonal all positive;
- symmetry / energy positivity;
- constant sigma vs manufactured scalar Laplace;
- error stable on irregular relaxed particles.

---

## 2.9 `electromagnetic_ophelie_phi.h/.hpp`

Implements:

- zero scalar;
- `PhiImag` RHS from `ASrcReal`;
- scalar phi gradient;
- gauge penalty;
- Jacobi update;
- `E/J/Q` with phi;
- `solvePhiImagJacobi`;
- residual evaluation.

Assessment: Level 1 phi correction body already present.

### Key Validation Points

Current RHS:

```cpp
rhs_i -= omega * sigma_ij * g_ij.dot(a_re_i - a_re_j);
```

Corresponds to:

```text
div(sigma grad PhiImag) = - div(omega sigma AReal)
```

Sign looks consistent with current convention, but must be locked by manufactured / consistency tests.

### Recommended Phi Unit Tests

1. `ASrcReal = constant vector`  
   RHS should be near 0; `PhiImag` near constant/0.

2. `PhiImag = constant`  
   Laplace should be near 0.

3. Prescribe `PhiImag = x` or `sin(pi x)`  
   Check gradient and Laplace direction/magnitude.

4. Compare `--no-phi` and phi correction:  
   After phi, `div(JImag)` should decrease, not just `passed=1`.

---

## 2.10 `electromagnetic_ophelie_phi_gmres.h`

Host-side GMRES and PCG:

- operator apply still via SYCL InteractionDynamicsCK;
- Krylov vectors on host;
- diagonal preconditioner;
- PCG/GMRES/Jacobi optional.

Assessment: sufficient for current-scale tests; not final high-performance version but acceptable.

### Note

PCG assumes SPD operator. Current pairwise Laplace + gauge should be near SPD, but SPH volume, adjacency, boundaries, and formulation may break strict symmetry; GMRES is safer reference solver.

Recommended test order:

```text
GMRES reference
PCG fast path
Jacobi smoke only
```

---

## 2.11 `electromagnetic_ophelie_self_induction.h`

Picard-like self-induction loop already implemented:

```text
combine A
solve phi
compute grad phi
compute J
store J
loop:
    compute A_ind from J
    combine
    solve phi
    update J
    under-relax J
```

Assessment: engineering path in place, but physics phase needs fix before trust.

Current status recommendation:

```text
--self-induction is experimental smoke test only.
Do not use for client reports or TEAM7 validation.
```

Fix items:

- `JImag -> AIndImag`, not `AIndReal`;
- sync `B` cross direction fix;
- retest self-induced field scaling and relaxation.

---

## 2.12 `electromagnetic_ophelie_team7_geometry.h`

Current TEAM7-like geometry:

```text
domain 1m cube
plate: cylinder r=0.36, half-thickness=0.02, center=(0.5,0.5,0.08)
coil: annulus r_in=0.12, r_out=0.18, half-height=0.12, center=(0.5,0.5,0.5)
air: no particles, only computational domain
```

Assessment:

- OK as OPHELIE-like smoke geometry;
- not strict TEAM7 standard validation;
- `air domain` is bounding box only, not `AirBody`; README “air domain” is more accurately “air-sized surrounding domain without air particles”.

Whether `AirBody` is needed:

- for current OPHELIE-like Biot-Savart source model, no `AirBody` needed;
- for visualizing `B/A` in air, add `AirSampleBody` as sampling points only, no phi;
- full A–phi/TEAM7 PDE comparison needs a real air domain.

---

## 2.13 `electromagnetic_ophelie_relaxation.h`

Uses SYCL CK relaxation pipeline, avoiding legacy device relax hangs. Assessment: correct.

Recommendations:

- preserve;
- for TEAM7-like geometry use `--skip-relax` for smoke first;
- use reload for production figures;
- relax VTP every 100 steps is enough.

---

## 2.14 `test_3d_ophelie`

Box integration test covering:

- box glass + box shell coil;
- coil source;
- Biot-Savart;
- combine A;
- Level 0;
- optional phi;
- optional self induction;
- Q scaling;
- VTP;
- passed criteria.

Assessment: most important regression test. Priority: pin as CTest or fixed run commands.

Gaps:

- `passed` checks max/nonzero and residual only; no direction, analytic values, divJ, scaling;
- will not catch reversed B;
- will not catch `JouleHeat` vs `J/E` scaling inconsistency.

---

## 2.15 `test_3d_ophelie_team7`

Already has TEAM7-like plate + annular coil + relaxation/reload + phi solve pipeline.

Assessment: suitable as next-stage geometry smoke test, but do not skip operator tests and treat as validation.

Recommended improvements:

1. Fix `J0` area formula;
2. Fix `B` cross direction;
3. Add `--no-phi` baseline first;
4. Add `--phi-solver=GMRES/PCG` comparison;
5. Add `--sigma` sweep;
6. Optionally add `AirSampleBody` for B/A visualization only;
7. Do not treat `--self-induction` as formal output until phase is fixed.

---

## 3. Current Progress Assessment

Already completed:

```text
Stage 0: Geometry and VTP output             Completed
Stage 1: coil source initialization          Completed
Stage 2: coil -> plate/glass Biot            Completed, but B direction needs fix
Stage 3: Level 0 E/J/Q                       Completed
Stage 4: PhiImag correction                  Implemented, needs unit tests locked in
Stage 5: self-induction iteration            Path implemented, physics phase needs fix
TEAM7-like geometry + reload relax           Prototype implemented
```

Therefore not “just starting” but at **operator-level validation and TEAM7-like geometry smoke in parallel**.

---

## 4. Run TEAM7-like Case Now?

Recommendation: **yes, start, but not TEAM7-like only.**

Correct order:

```text
A. Fix two key issues first:
   1. B cross direction
   2. TEAM7 coil J0 area

B. Run minimal operator tests:
   1. coil source direction
   2. Biot-Savart B direction/sign
   3. I0 scaling: A/B ~ I, Q ~ I^2
   4. Level 0 uniform A analytic
   5. Laplace constant nullspace
   6. Phi RHS zero for constant A

C. Then TEAM7-like:
   1. --skip-relax --no-phi
   2. --skip-relax with phi GMRES/PCG
   3. relax=1 generate reload
   4. reload=1 run EM
```

TEAM7 is integration smoke, not a substitute for operator tests.

---

## 5. Recommended Additional Test Manifest

## 5.1 `test_3d_ophelie_coil_source_direction`

Goals:

- validate `e_theta` direction;
- validate `J0` scaling.

Acceptance:

```text
x>center, y=center -> J_y > 0
x=center, y>center -> J_x < 0
I0 doubled -> max(JSrcReal) doubled
```

---

## 5.2 `test_3d_ophelie_biot_savart_direction`

Goal: catch B cross-product sign.

Setup:

```text
source at (R,0,0)
J = (0,+J,0)
observer at origin
```

Expected:

```text
Bz > 0
```

Fails if code uses `r.cross(J)`.

---

## 5.3 `test_3d_ophelie_scaling`

Goal:

```text
I0 -> 2 I0
A/B/E/J -> 2x
raw Joule power -> 4x
```

Check before power normalization.

---

## 5.4 `test_3d_ophelie_level0_uniform_A`

Prescribe `ASrcReal = constant Vecd(A0,0,0)`, skip Biot:

```text
EImag = -omega * ASrcReal
JImag = sigma * EImag
Q = 0.5 * sigma * |EImag|^2
```

---

## 5.5 `test_3d_ophelie_phi_laplace_constant`

Goal:

```text
PhiImag = constant
Laplace(PhiImag) ~ 0
```

Also check diag positive.

---

## 5.6 `test_3d_ophelie_phi_rhs_constant_A`

Goal:

```text
ASrcReal = constant
RHS = -div(omega sigma A) ~ 0
```

If RHS not near zero, indicates unstable `g_ij` or sign/boundary handling.

---

## 5.7 `test_3d_ophelie_phi_reduces_divJ`

Goal:

Compare:

```text
Level 0: JImag = -sigma omega AReal
Phi:    JImag = sigma(-gradPhiImag - omega AReal)
```

Acceptance:

```text
||divJ_phi|| < ||divJ_level0||
```

Core test that phi correction actually helps.

---

## 6. TEAM7-like Test Schedule

## 6.1 Round 1: No relax, quick smoke

```bash
cd ~/SPHinXsysSYCL/build
ninja test_3d_ophelie_team7

./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_team7/bin/test_3d_ophelie_team7 \
  --state_recording=1 --skip-relax --no-phi --sigma=1e4
```

Check:

```text
passed=1
max_ASrc > 0
max_BSrc > 0
max_EImag > 0
max_JImag > 0
P_scaled ≈ 50000
VTP: JSrcReal tangential about coil
VTP: BSrc direction correct at coil center
```

## 6.2 Round 2: Phi correction

```bash
./tests/.../test_3d_ophelie_team7 \
  --state_recording=1 --skip-relax --phi-solver=GMRES --sigma=1e4
```

Then:

```bash
./tests/.../test_3d_ophelie_team7 \
  --state_recording=1 --skip-relax --phi-solver=PCG --sigma=1e4
```

Check:

```text
phi_rel_res < 10 * tolerance
max_PhiImag > 0
max_PhiRhs > 0
JouleHeat >= 0
```

After divJ diagnostic, require clear divJ drop after phi.

## 6.3 Round 3: Relax/reload

Generate reload:

```bash
./tests/.../test_3d_ophelie_team7 \
  --relax=1 --state_recording=0 --relax-steps=400 --relax-vtp-every=100
```

Run with reload:

```bash
./tests/.../test_3d_ophelie_team7 \
  --reload=1 --state_recording=1 --sigma=1e4
```

Check:

```text
Body names: PlateBody, CoilSourceBody
particle count stable
VTP OK
passed=1
```

## 6.4 Air Domain?

No real `AirBody` needed currently.

If client/advisor wants “field in air region,” add:

```text
AirSampleBody
```

Sampling `ASrc/BSrc` only; no phi; no relax.

---

## 7. Answers to Current Cursor Questions

## 7.1 Can Reference Code Directory Be Deleted?

Yes, after:

1. Current formal code is git add / commit;
2. Development plan docs moved to `docs/ophelie/` or other formal location.

Recommendation:

```bash
mkdir -p docs/ophelie
cp FRENCH_OPHELIE_PARTICLE_DEVELOPMENT_PLAN.md docs/ophelie/
git add docs/ophelie/FRENCH_OPHELIE_PARTICLE_DEVELOPMENT_PLAN.md
```

Then:

```bash
rm -rf biot_reference_selection_for_cursor
rm -rf french_ophelie_particle_plan_and_starter
```

## 7.2 Can build_ophelie / build_ophelie_cpu Be Deleted?

Yes. Keep current in-progress `build/` only.

```bash
rm -rf build_ophelie build_ophelie_cpu
```

## 7.3 Formal Code Should Enter Git First

Recommend immediately:

```bash
git add tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie
git add tests/extra_source_and_tests/3d_examples/test_3d_ophelie
git add tests/extra_source_and_tests/3d_examples/test_3d_ophelie_team7
git add docs/ophelie/FRENCH_OPHELIE_PARTICLE_DEVELOPMENT_PLAN.md
git status
```

Then initial commit.

---

## 8. Clear Next Tasks for Cursor

### Required Fixes

```text
1. Fix Biot-Savart B cross product:
   r.cross(jv) -> jv.cross(r)

2. Fix TEAM7 coil J0:
   use (Rout - Rin) * (2 * half_height) as current cross-section area,
   not annular cylinder volume.

3. Clarify power scaling:
   either add JouleHeatRaw/JouleHeatScaled,
   or scale E/J by sqrt(scale) for self-consistent output.

4. Mark self-induction experimental for now.
   If continuing, fix JImag -> AIndImag phase channel.
```

### Required Tests

```text
1. coil source direction
2. Biot-Savart B direction/sign
3. I0 scaling
4. Level0 uniform A analytic
5. Laplace constant nullspace
6. Phi RHS constant A
```

### TEAM7-like Can Proceed in Parallel

```text
1. --skip-relax --no-phi first
2. then --skip-relax + phi
3. then relax/reload
4. then consider AirSampleBody
```

---

## 9. Summary

Current code already implements the OPHELIE-like particle EM heating main chain in SPHinXsys/SYCL form, well ahead of the starter plan. Priority is not more geometry but fixing Biot-Savart direction, J0 dimensions, power scaling consistency, and adding small operator tests. After those, TEAM7-like relax/reload is more reliable.
