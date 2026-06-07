# A-phi CK Matrix-Free Solver: Stage 9C-3 Review and Stage 9D Plan

## Executive conclusion

The current Stage 9A--9C development is on the right track and remains consistent with the SPHinXsys CK/SYCL programming model. The production path is now:

- operator: standard SPH `Vol_j` matrix-free A--phi operator;
- solver: restarted right-preconditioned GMRES;
- default preconditioner: local block-Jacobi;
- gauge handling: phi penalty regularization;
- BiCGStab/PCG: diagnostic only.

The Stage 9C-3 TEAM7-like test at full conductivity contrast (`sigma_air=1e-4`, `sigma_conductor=1`) reaching about `5e-4` but not `1e-5` should be interpreted primarily as a solver/preconditioner limitation, not as operator plumbing failure. The homogeneous impressed-current case, the variable-sigma discrete MMS case, and the scalar/vector/full manufactured GMRES convergence tests already support the correctness of the CK operator/solver pipeline.

## Code-structure review

The new code is still following the SPHinXsys CK/SYCL style:

- material/source setup uses `StateDynamics<MainExecutionPolicy, ...>`;
- matrix-free operator application uses `AphiApplyDynamicsBundle<ExecutionPolicy>` and `InteractionDynamicsCK`-style kernels;
- reductions use `ReduceDynamicsCK`;
- per-particle state is registered through SPHinXsys discrete variables;
- host-side work is limited to small GMRES Hessenberg least-squares solves and benchmark post-processing.

This is acceptable. A small dense Eigen solve for the GMRES Hessenberg system is not a violation of the CK/SYCL design, because the large particle-vector operations remain in CK kernels.

## Important implementation issue to fix immediately

In `AphiGMRESSolverCK::solve()`, when the solver exits because `max_outer_iterations` is reached, the reported `final_relative_residual` may be stale. The residual is recomputed at the beginning of each outer cycle, then an Arnoldi cycle updates the solution. If the loop ends immediately after that update, the result still contains the residual from before the last update.

Fix this by recomputing residual and true residual after the loop before returning `MaxOuterIterationsReached`:

```cpp
apply_solution_to_lhs_.exec();
compute_recursive_residual_.exec();
const Real final_norm = std::sqrt(norm_squared_r_.exec());
result.final_residual_norm = final_norm;
result.final_relative_residual = final_norm / initial_norm;
recomputeFinalTrueResidual(result, initial_norm);

if (result.final_relative_residual < solver_options_.relative_tolerance ||
    final_norm < solver_options_.absolute_tolerance)
{
    result.converged = true;
    result.breakdown = false;
    result.breakdown_code = AphiGMRESBreakdownCode::None;
    return result;
}

result.converged = false;
result.breakdown = false; // max iteration is not numerical breakdown
result.breakdown_code = AphiGMRESBreakdownCode::MaxOuterIterationsReached;
return result;
```

This matters especially for `team7_like_coil_plate_tight_tol_diagnostic`, because the current reported residual floor may be pessimistic by one final GMRES update.

## Answer to Question 1: why does TEAM7-like full contrast stall near 5e-4?

The most likely cause is the combination of high conductivity contrast, restarted GMRES, and a weak local block-Jacobi preconditioner.

It is less likely to be a basic discretization/plumbing error because:

- homogeneous impressed-current benchmark converges;
- variable-sigma discrete MMS converges;
- rough adjointness/Rayleigh/tiny-matrix diagnostics were already acceptable;
- true residual and recursive residual agree in the GMRES tests.

It is also not mainly a continuation issue. Sigma or omega continuation can reduce the initial residual by warm-starting, but it does not change the final operator spectrum. If the final high-contrast operator is poorly preconditioned, continuation will not remove the late-stage residual plateau.

The current block-Jacobi preconditioner includes only local Laplace/reaction/phi-penalty diagonal information. It ignores:

1. neighbor-level high-contrast interface coupling;
2. off-diagonal A--phi grad/div coupling;
3. low-frequency/global error modes across air/coil/conductor regions;
4. material-region imbalance in the residual norm;
5. non-normality introduced by reaction and A--phi coupling.

This is exactly the kind of situation where restarted GMRES can decrease monotonically but very slowly, or appear to plateau near a residual level determined by the missing spectral information in the preconditioner.

## Answer to Question 2: what should be tried first to reach 1e-5?

Do not change the default Laplace weight first. The priority order should be:

### Step 1: fix the GMRES final-residual reporting issue

Before drawing any final conclusion about the `5e-4` floor, fix the max-outer final residual recomputation described above and rerun the tight diagnostic.

### Step 2: allow larger restart dimensions

The current `AphiGMRESWorkspaceCK` hard-caps `AphiGMRESMaxRestartDimension = 30`, and `AphiMatrixFreeSolveCK` registers a max-30 workspace regardless of the requested solver option. This prevents testing whether the plateau is caused by restart truncation.

Refactor the workspace registration so that `AphiMatrixFreeSolveCK` registers exactly the requested restart dimension:

```cpp
gmres_workspace_(buildAphiGMRESWorkspaceNames(solver_options.gmres.restart_dimension)),
register_gmres_workspace_(sph_body, solver_options.gmres.restart_dimension)
```

Then test:

- `m = 30`, `50`, `80`;
- `max_outer = 100` or more;
- same TEAM7-like full-contrast setup;
- report final true residual after the final update.

If increasing `m` significantly lowers the residual, the main issue is restart truncation rather than the preconditioner alone.

### Step 3: add a no-PC and PC comparison specifically for TEAM7-like full contrast

The existing no-PC comparison was mostly for simpler manufactured/scalar cases. Add a TEAM7-like diagnostic:

- GMRES(m=30), no PC;
- GMRES(m=30), block-Jacobi PC;
- GMRES(m=50/80), block-Jacobi PC;
- same initial RHS and tolerance.

This tells whether block-Jacobi helps, hurts, or only helps early iterations.

### Step 4: implement a coupled point-block Jacobi preconditioner

The next real preconditioner should still be local and CK-friendly. Recommended form: a per-particle coupled block inverse for the real-valued A--phi block.

For each particle, the local unknown vector is approximately:

```text
[A_real_x, A_real_y, A_real_z,
 A_imag_x, A_imag_y, A_imag_z,
 phi_real, phi_imag]
```

The current PC only handles the A real/imag reaction block and the phi scalar diagonal. The improved PC should also include the local diagonal contributions from grad(phi) and div(sigma A) coupling:

- A equation depends locally on `phi_i` through `sigma_i sum_j g_ij phi_i`;
- phi equation depends locally on `A_i` through `omega sum_j sigma_ij g_ij dot(A_i)`.

This gives an 8x8 local block per particle. The block can be assembled inside a CK local update and inverted locally. As a first implementation, only include the diagonal/local part of the coupling, not neighbor off-diagonal coupling.

This is likely the best next PC because it directly targets the missing A--phi coupling in the current block-Jacobi PC.

### Step 5: consider residual/block scaling

If the coupled local PC is not enough, add block scaling/equilibration:

- separate scaling for A and phi components;
- possibly material-region residual scaling for air/coil/conductor;
- report component-wise residual norms, not only the total block norm.

This will show whether the `5e-4` residual is dominated by phi, A in air, A in conductor, or interface particles.

## Answer to Question 3: should Stage 9D thermal coupling wait for 1e-5?

Split Stage 9D into two levels.

### 9D-lite can start at 5e-4

It is reasonable to start a limited Joule-heating plumbing stage at the current `5e-4` level:

- compute induced electric-field proxy from A--phi;
- compute Joule source in conductor/coil regions;
- output total Joule power and regional source distributions;
- verify source non-negativity and energy accounting;
- do not claim quantitative validated heating yet.

This is useful because it tests data flow from EM to thermal variables and identifies what field quantities must be saved.

### Full thermal coupling should wait

Do not start full temperature-dependent sigma/nu feedback, cold-crucible heating, or paper-level quantitative heat-transfer claims until one of these is true:

1. GMRES reaches around `1e-5` on the full-contrast TEAM7-like case with the improved PC; or
2. a tolerance sweep proves the Joule heating integral is insensitive to EM residual, for example:
   - compare EM tolerances `5e-3`, `5e-4`, `1e-4`, and if possible `1e-5`;
   - require total conductor Joule power to vary by less than about 1--2% between the accepted tolerances.

Therefore: start 9D-lite now, but keep full two-way thermal coupling blocked until stronger PC or observable-level tolerance evidence exists.

## Answer to Question 4: does 9C-2 need an interface-aware manufactured solution?

For solver and operator plumbing, the current discrete MMS is enough. It validates that the same discrete operator can recover a field from `b = K(u_exact)` in the presence of a sigma jump.

However, this is not a physical interface-condition validation. A future interface-aware manufactured solution is still valuable before publication-level claims. It should use piecewise analytic fields satisfying at least:

- continuity of phi or the chosen potential variable;
- continuity of normal conductive/electromagnetic flux, e.g. `sigma partial_phi/partial_n` for the scalar part;
- compatible source terms on each side of the interface.

Suggested next validation after Stage 9D-lite:

- a 1D or quasi-1D two-material interface patch test;
- piecewise linear or sinusoidal phi with flux-continuity enforced;
- report residual localization near the interface and convergence with dp.

So: not required before the next engineering development step, but should be planned before a formal paper/benchmark claim.

## Recommended Stage 9D / 9E order

### Stage 9D-0: GMRES/PC cleanup ✅ (2026-05-21)

1. ✅ Fix max-outer final residual recomputation.
2. ✅ Make GMRES restart dimension configurable beyond 30 (max 80).
3. ✅ Add TEAM7-like restart sweep: `m=30/50/80`.
4. ✅ Add TEAM7-like no-PC vs decoupled vs coupled diagnostic.
5. ✅ Add component-wise + region-wise residual diagnostics.

### Stage 9D-1: improved PC ✅ (2026-05-21)

1. ✅ Implement coupled point-block Jacobi 8×8 with local A--phi coupling.
2. ✅ Compare on TEAM7-like full contrast (+ solver study on decoupled baseline).
3. ⚠️ Target `1e-5` true relative residual: **best ~1.96e-5** (coupled m=80, outer=150); not fully met.

### Stage 9D-lite: one-way Joule heating plumbing ✅ (2026-05-21)

1. ✅ Compute E proxy from frequency-domain A--phi.
2. ✅ Compute conductor/coil Joule source + regional power.
3. ✅ Verify source non-negativity.
4. ✅ **Next:** tolerance sensitivity — **done (9D-2)**; conductor Joule stable <0.01% for tol 5e-4→~2e-5.

### Stage 9D-2: Joule EM-tolerance sensitivity ✅ (2026-05-21)

1. ✅ Sweep tol = 5e-3 / 5e-4 / 1e-4 / 1e-5 on TEAM7-like case.
2. ✅ Record total/conductor Joule power + EM rel residual.
3. ✅ Pass if conductor power variation vs 5e-4 ref ≤ 2% — **actual ~0.005%**.

### Stage 9E: stronger validation ✅ (2026-05-21)

1. ✅ interface-aware manufactured solution (flux-matched piecewise φ);
2. ✅ real TEAM7 dimensions smoke (1.2×1.0×0.3 m, fraction layout);
3. ✅ convergence with dp (0.2 / 0.1 / 0.075 diagnostic);
4. only then move toward cold-crucible geometry and thermal feedback.

**Record:** `tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_STAGE9E_VALIDATION_RECORD.md`

**Key results:**
- 9E-1: discrete defect ~5e-6, continuous error ~2e-6 (vs 9C-2 separable O(1) continuous error)
- 9E-2: physical box @ tol=5e-4, outer=2, true_rel ~4.7e-4
- 9E-3: dp=0.1 recommended; dp=0.2 too coarse to converge

**1e-5 status:** not closed at full σ contrast (best ~2e-5); **engineering acceptance at ~2e-5** documented in 9D-2.

## Do not do yet

- Do not replace the default `Vol_j` Laplace weight.
- Do not add Contact.
- Do not enter full cold-crucible geometry yet.
- Do not claim the current 9C-3 test is quantitatively validated TEAM7.
- Do not treat `5e-4` as a universal engineering tolerance until Joule-power sensitivity is checked.
