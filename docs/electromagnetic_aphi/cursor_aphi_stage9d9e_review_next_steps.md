# A-phi CK Matrix-Free Solver: Stage 9D/9E Review and Next Steps

## Overall judgement

Stage 9D and Stage 9E are a meaningful milestone. The current code remains aligned with the SPHinXsys CK/SYCL design: particle-level operations are implemented through `LocalDynamics`, `Interaction<Inner<>>`, `InteractionDynamicsCK`, `StateDynamics`, `ReduceDynamicsCK`, and registered discrete variables. The only substantial host-side dense computation is the small GMRES Hessenberg least-squares solve, which is acceptable because the large particle-vector operations remain in CK kernels.

The project has now moved from “can the A-phi matrix-free operator be assembled and solved?” to “can the solver and validation chain support quantitative TEAM7/cold-crucible claims?”

Current production default:

```text
AphiMatrixFreeSolveCK
+ restarted right-preconditioned GMRES
+ restart dimension m = 50
+ CoupledPointBlock8x8 local preconditioner
+ phi gauge penalty
+ Standard SPH Vol_j pairwise Laplace weight
```

The default Laplace weight should not be changed at this stage.

---

## Stage 9D/9E key conclusions

### Stage 9D-0: GMRES/PC cleanup

Completed items:

- GMRES final residual recomputation after max-outer exit was fixed.
- GMRES restart dimension was extended up to 80.
- TEAM7-like restart sweep was performed.
- no-PC / decoupled PC / coupled PC comparison was added.
- component-wise and region-wise residual diagnostics were added.

Main conclusion:

```text
m=30 was too restrictive for the high-contrast TEAM7-like case.
Increasing m to 50 or 80 reduces the residual from the old ~5e-4 plateau to ~2e-5.
```

This shows that the old residual floor was partly a reporting/restart limitation, not a hard operator failure.

### Stage 9D-1: CoupledPointBlock8x8 PC

The coupled 8x8 point-block Jacobi PC is a reasonable CK-friendly next PC. It includes local A-phi coupling and is still embarrassingly parallel at particle level.

However, it does not yet break the full-contrast TEAM7-like case below 1e-5:

```text
best current result: true relative residual ≈ 1.96e-5
```

This is close enough for engineering use, especially after Joule sensitivity testing, but not enough to claim strict 1e-5 closure.

### Stage 9D-lite and 9D-2: Joule source plumbing and tolerance sensitivity

The one-way EM-to-Joule data path is now implemented:

```text
grad phi -> frequency-domain E proxy -> current density -> Joule source -> regional Joule power
```

The important result is:

```text
Conductor Joule power changes by < 0.01% when EM tolerance is tightened from 5e-4 to the best available ~2e-5.
```

This supports accepting the current ~2e-5 EM residual level for engineering Joule-power diagnostics.

### Stage 9E-1: Flux-matched interface MMS

The new interface-aware MMS is a significant improvement over the earlier separable variable-sigma MMS:

```text
continuous field error ≈ 1.8e-6
interface-band residual ≈ 1.9e-6
true relative residual ≈ 5e-6
```

This closes a major validation gap: the variable-sigma interface is no longer validated only by a discrete self-consistency test.

### Stage 9E-2: Physical TEAM7 dimensions smoke

The physical-size box smoke test runs successfully:

```text
box: 1.2 x 1.0 x 0.3 m
true relative residual ≈ 4.7e-4 at engineering tolerance
```

This is a geometry/scaling smoke test, not yet a quantitative TEAM7 validation.

### Stage 9E-3: dp sweep

Current dp sweep:

```text
dp = 0.2   too coarse, fails
dp = 0.1   recommended baseline
dp = 0.075 passes but needs more outer restarts
```

The conclusion “dp=0.1 is the current engineering baseline” is acceptable for the simplified TEAM7-like benchmark. It is not yet a formal convergence result for a physical observable.

---

## Code style and framework compliance

The new code is still consistent with SPHinXsys CK/SYCL style.

Positive points:

1. Large particle-vector operations remain CK kernels.
2. New fields are registered as discrete variables.
3. The Joule path is implemented through `InteractionDynamicsCK` and `StateDynamics`.
4. The coupled 8x8 PC is local per particle and therefore CK-friendly.
5. The high-level solver wrapper `AphiMatrixFreeSolveCK` is the right abstraction for production.
6. Diagnostic solvers are clearly separated from the production GMRES path.

Acceptable host-side parts:

- Eigen least-squares solve for GMRES Hessenberg matrix.
- Host post-processing for benchmark metrics.

These do not violate the CK/SYCL design because the heavy vector operations remain on the CK/SYCL path.

---

## Important issues and recommended fixes

### 1. Add diagnostics for 8x8 PC fallback and conditioning

`aphiSolve8x8()` can fail and silently fall back to the decoupled PC. That is acceptable for robustness, but it must be diagnosable.

Add optional diagnostic outputs:

```text
number of 8x8 solve failures
minimum pivot magnitude
maximum local condition proxy
region-wise fallback count: air / coil / conductor
```

At minimum, add a counter or a diagnostic variable such as:

```text
Jacobi8x8FallbackFlag
Jacobi8x8MinPivot
```

Then add a test or log in the TEAM7-like solver study.

Reason: if coupled PC is not improving much beyond decoupled PC, silent 8x8 fallback or ill-conditioned local blocks could be part of the reason.

### 2. Keep Joule heating labeled as “operator-consistent proxy” until validated

The implemented Joule path is appropriate for 9D-lite, but before using it for quantitative heating claims, verify the field formula and sign convention against a simple analytic case.

Current formula:

```text
E_real = omega * A_imag - grad(phi_real)
E_imag = -omega * A_real - grad(phi_imag)
q = 0.5 * sigma * (|E_real|^2 + |E_imag|^2)
```

This is plausible for the selected time-harmonic convention, but the gradient stencil sign must be verified against the A-phi operator convention.

Add a small Joule verification case:

```text
uniform conductor
A = 0
phi = -E0 x
expected Joule source = 0.5 sigma |E0|^2
```

This should be done before full thermal coupling.

### 3. Do not overinterpret Joule tolerance sensitivity

The current tolerance sensitivity result is strong for the simplified TEAM7-like setup, but it is based on total/regional Joule power. Before publication-level claims, also compare:

```text
L2 difference of Joule field
max-norm difference of Joule field
interface/conductor near-surface Joule difference
integrated conductor power
```

Total power can be insensitive while local source distribution still changes.

### 4. Physical TEAM7 dimension smoke is not quantitative TEAM7 validation

The 1.2 x 1.0 x 0.3 m smoke test confirms that scaling and thin-domain setup run. It is not yet the real TEAM7 benchmark.

Do not call it “TEAM7 validation” yet. Call it:

```text
TEAM7-dimension smoke test
TEAM7-like physical-size solver benchmark
```

Quantitative validation still needs reference geometry, reference material parameters, source current definition, and comparison metrics.

### 5. dp sweep should move from residual-based to observable-based

The current dp sweep uses solver residual and field maxima. That is useful, but the next sweep should report observable convergence:

```text
conductor Joule power
coil Joule power
max/mean Joule source in conductor
field norm in conductor
field norm near conductor surface/interface
```

This is more relevant than residual alone for cold-crucible applications.

### 6. GMRES orthogonality may become an issue for m=80

Current GMRES uses modified Gram-Schmidt-like orthogonalization. For m=80 and high-contrast systems, loss of orthogonality may become relevant.

Recommended diagnostic:

```text
max |<v_i, v_j>| for i != j after Arnoldi
```

If this grows, add one extra reorthogonalization pass or switch to a Givens-rotation GMRES residual update.

This is not urgent, but it becomes important if m > 80 or if strict 1e-5 closure is required.

---

## Should we force 1e-5 now?

No, not as the main priority.

Current best full-contrast TEAM7-like result is around:

```text
true relative residual ≈ 2e-5
```

The Joule-power sensitivity result shows that tightening EM tolerance from 5e-4 to ~2e-5 changes integrated conductor Joule power by less than 0.01%.

Therefore:

```text
Engineering EM tolerance: accept ~2e-5 for Joule-power diagnostics.
Strict 1e-5: optional research target, not a blocker for the next development stage.
```

Do not spend the next phase only trying to reduce 2e-5 to 1e-5 unless there is a specific benchmark or reviewer-driven need.

---

## Current project status

The A-phi CK solver is now at the following maturity level:

```text
Core matrix-free operator: mature prototype
GMRES solver path: production default for tests
Gauge handling: functional via phi penalty
Variable-sigma interface validation: substantially improved
TEAM7-like engineering benchmark: functional
Joule source plumbing: functional diagnostic path
Full thermal feedback: not started
Quantitative TEAM7 validation: not started
Cold-crucible geometry: not started
```

A reasonable maturity estimate:

```text
A-phi solver prototype: 80-85%
Publication-method prototype: 65-70%
Cold-crucible engineering tool: 35-45%
```

---

## Recommended next development sequence

### Stage 10A: finalize EM/Joule verification before full thermal coupling

Do this first:

1. Add 8x8 PC fallback/conditioning diagnostics.
2. Add Joule analytic verification with uniform E field.
3. Add observable-based tolerance sensitivity:
   - conductor Joule power,
   - local Joule max,
   - Joule L2 field difference,
   - conductor/interface band metrics.
4. Add physical-dimension dp sweep for the TEAM7-dimension smoke case.

Goal: make the EM-to-Joule pipeline defensible.

### Stage 10B: quantitative TEAM7 preparation

Then prepare true TEAM7 validation:

1. Define exact TEAM7 geometry and dimensions.
2. Define coil current/source convention.
3. Define material parameters and frequency.
4. Choose comparison metrics:
   - induced current/Joule distribution in plate,
   - total plate Joule power,
   - field magnitude along sampling lines,
   - region-integrated EM energy or source response.
5. Compare against reference data or an FEM/MFEM baseline.

This is the natural next paper-quality validation step.

### Stage 10C: cold-crucible geometry scaffold

After TEAM7 setup is stable, add cold-crucible geometry scaffold:

1. Geometry and region tagging only.
2. Coil/source tagging.
3. Conductor/melt/crucible/air material tagging.
4. EM-only solve.
5. Joule source output.
6. No temperature feedback at first.

### Stage 10D: one-way thermal coupling

Only after EM-only cold-crucible solve is stable:

1. Map Joule source into thermal equation.
2. Run one-way thermal response with fixed sigma/nu.
3. Check energy accounting.
4. Verify heat source localization.

### Stage 10E: two-way thermal feedback

This should be last:

1. sigma(T), nu(T) update.
2. EM/thermal staggered coupling.
3. relaxation or under-relaxation.
4. convergence of thermal feedback loop.

---

## Recommended message back to Cursor

Use this as the next instruction block:

```text
Stage 9D and 9E are accepted as a milestone. The code still follows the SPHinXsys CK/SYCL style: large particle operations are implemented through LocalDynamics, InteractionDynamicsCK/StateDynamics, ReduceDynamicsCK and registered discrete variables; host Eigen is only used for the small GMRES Hessenberg least-squares solve and is acceptable.

Do not keep pushing strict 1e-5 as the main blocker. The current ~2e-5 full-contrast TEAM7-like residual, together with the Joule-power tolerance sensitivity (<0.01% variation from 5e-4 to ~2e-5), is sufficient for engineering EM/Joule diagnostics.

Next stage should not be full two-way thermal coupling yet. First do Stage 10A:
1. Add 8x8 coupled PC fallback/conditioning diagnostics.
2. Add an analytic Joule verification case: uniform conductor, A=0, phi=-E0 x, expected q=0.5 sigma |E0|^2.
3. Extend Joule tolerance sensitivity from total power to local/source-field metrics: Joule L2 field difference, max Joule difference, conductor/interface band Joule metrics.
4. Run physical TEAM7-dimension dp sweep and report observable convergence, not only residual convergence.

After Stage 10A, proceed to quantitative TEAM7 preparation. Cold-crucible geometry scaffold should come after the TEAM7 validation path is defined, and full thermal feedback should remain deferred until EM-only Joule source output is stable.
```
