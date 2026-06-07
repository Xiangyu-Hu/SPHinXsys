# A-phi CK Matrix-Free Solver: Stage 9 Plan after Stage 8C

## 1. Current conclusion

Stage 8C is a major milestone. The scalar phi Laplace + penalty, vector Helmholtz, and full A-phi + penalty manufactured systems all converge with restarted right-preconditioned GMRES(m), while BiCGStab and PCG still diverge on the same operator family. Therefore the current evidence supports the following conclusion:

> The matrix-free operator plumbing is basically correct; the remaining failure is primarily a Krylov-method/preconditioner compatibility issue, not an operator assembly or RHS/residual wiring issue.

The default Laplace weight should not be changed at this stage. The current diagnostics show rough-field adjointness around 1e-7, positive random Rayleigh quotients, compatible constant-mode RHS for pure Laplace, and a small host-assembled scalar matrix whose weighted symmetric part is SPD after phi penalty.

## 2. Immediate policy decision

Use restarted right-preconditioned GMRES as the default solver for the A-phi CK matrix-free path.

Keep BiCGStab and PCG as diagnostic solvers, not as the production/default solver. PCG is useful only for scalar SPD diagnostics. BiCGStab can remain as an experimental non-symmetric Krylov method but should not block the A-phi development path.

## 3. Code cleanup required before Stage 9 applications

### 3.1 Rename tests and solver roles

Recommended names / comments:

- `test_3d_aphi_ck_pcg_scalar_phi_laplace_penalty_diagnostic`: diagnostic only; not a required convergence test.
- `test_3d_aphi_ck_bicgstab_scalar_phi_laplace_penalty_convergence`: mark as expected-failing or diagnostic until BiCGStab is revisited.
- `test_3d_aphi_ck_gmres_scalar_phi_laplace_penalty_convergence`: production Krylov baseline.
- `test_3d_aphi_ck_gmres_vector_helmholtz_convergence`: production Krylov baseline.
- `test_3d_aphi_ck_gmres_full_aphi_penalty_convergence`: production Krylov baseline for the coupled A-phi block.

Do not present BiCGStab/PCG failures as blocking failures once GMRES convergence is established.

### 3.2 Add an AphiKrylovSolverKind enum

Introduce a small dispatch layer so tests and later examples can choose the solver without duplicating setup code:

```cpp
enum class AphiKrylovSolverKind
{
    GMRES,
    BiCGStabDiagnostic,
    PCGScalarDiagnostic
};
```

For the production A-phi path, default to:

```cpp
AphiKrylovSolverKind::GMRES
```

### 3.3 Add a compact solver wrapper

Create a high-level wrapper, for example:

```cpp
template <class ExecutionPolicy>
class AphiMatrixFreeSolveCK
{
  public:
    AphiMatrixFreeSolveCK(SPHBody &body,
                          Inner<> &inner_relation,
                          const AphiVariableNames &names,
                          const AphiLhsAssemblyOptions &operator_options,
                          const AphiMatrixFreeSolverOptions &solver_options);

    AphiMatrixFreeSolverResult solve();
};
```

The wrapper should own or receive:

- `AphiGMRESWorkspaceNames`
- GMRES workspace registration requirement
- solver options
- final result reporting

This wrapper should call `AphiGMRESSolverCK` by default.

### 3.4 Keep GMRES right-preconditioned semantics explicit

Documentation comment:

```cpp
// Restarted right-preconditioned GMRES for K(X)=B.
// z_j = M^{-1} v_j, w = K(z_j), and X <- X + sum_j y_j z_j.
// The recursive residual is recomputed at the beginning of every restart.
```

### 3.5 Add true-residual verification to GMRES tests

The current tests already print true residual in some cases. Make it part of the pass condition for the three GMRES convergence tests:

```cpp
const bool passed = result.converged &&
                    !result.breakdown &&
                    result.final_relative_residual < tolerance &&
                    result.final_true_relative_residual < 10.0 * tolerance &&
                    result.final_recursive_true_gap < 1.0e-4;
```

If `final_true_relative_residual` is not updated on the exact final iteration, force a final true-residual recomputation before returning `converged=true`.

## 4. GMRES implementation improvements

### 4.1 Always recompute final true residual before returning convergence

In `AphiGMRESSolverCK::solve()`, before any convergence return, call a final true residual recomputation if `recompute_true_residual=true`. This avoids reporting a stale `final_true_relative_residual` from the previous interval.

Pseudo-code:

```cpp
if (relative_residual < solver_options_.relative_tolerance ||
    beta < solver_options_.absolute_tolerance)
{
    if (solver_options_.recompute_true_residual)
    {
        apply_solution_to_lhs_.exec();
        compute_true_residual_.exec();
        const Real true_norm = std::sqrt(norm_squared_true_r_.exec());
        result.final_true_residual_norm = true_norm;
        result.final_true_relative_residual = true_norm / initial_norm;

        form_residual_gap_.exec();
        const Real gap_norm = std::sqrt(norm_squared_gap_.exec());
        result.final_recursive_true_gap = gap_norm / (true_norm + TinyReal);
    }
    result.converged = true;
    return result;
}
```

### 4.2 Add optional modified Gram-Schmidt reorthogonalization

Restarted GMRES currently uses one-pass modified Gram-Schmidt. Add an optional second pass:

```cpp
bool use_reorthogonalization = true;
```

When enabled, run a second MGS pass and accumulate the correction into `H(i,j)`. This is useful for larger domains and more difficult A-phi coupling.

### 4.3 Add estimated residual from the least-squares solve

After QR solve,

```cpp
residual_estimate = || beta e1 - H y || / initial_norm;
```

Store it in the result. This helps compare Arnoldi estimated residual and true residual.

### 4.4 Add restart-size sensitivity tests

Add a diagnostic test or command-line option to compare:

```text
m = 10, 20, 30
```

Expected behavior:

- scalar phi should converge even with smaller m;
- vector Helmholtz and full A-phi may need m=20 or 30;
- m=30 can remain the default.

### 4.5 Add no-preconditioner GMRES diagnostic

This is not required for production, but useful for understanding the block-Jacobi effect:

```text
GMRES(m) + no PC
GMRES(m) + block-Jacobi PC
```

Record outer iterations, Arnoldi steps, and final true residual. If no-PC GMRES also converges but much slower, this confirms block-Jacobi is helpful for GMRES even though it does not rescue PCG/BiCGStab.

## 5. Do not prioritize SymmetricPairVolume now

Do not implement or switch to a symmetric pair-volume weight as the next main step. It is lower priority because:

1. rough-field adjointness is already around 1e-7;
2. random Rayleigh quotients are positive;
3. the tiny scalar host matrix with penalty is SPD;
4. GMRES converges on scalar, vector, and full A-phi systems;
5. changing the default SPH Vol_j operator would alter the physical/discrete operator and force revalidation.

If implemented later, it must be optional:

```cpp
enum class AphiLaplaceWeightMode
{
    StandardSPHVolJ,
    SymmetricPairVolumeDiagnostic
};
```

Default must remain `StandardSPHVolJ`.

## 6. Stage 9 development route

### Stage 9A: consolidate solver infrastructure

- Keep GMRES as the default matrix-free Krylov solver.
- Add the high-level solver wrapper.
- Add final true-residual pass conditions.
- Add restart-size and no-PC diagnostic tests.
- Keep BiCGStab/PCG as diagnostic, not blocking.

### Stage 9B: manufactured-solution robustness sweep

Run GMRES convergence for:

- `dp = 0.2, 0.1, 0.075` if affordable;
- `omega = 0.1, 1.25, 10`;
- `sigma = 0, 0.1, 2.0` where meaningful;
- `phi_gauge_penalty = 1, 10, 100`.

Record:

- initial residual norm;
- final recursive residual;
- final true residual;
- Arnoldi steps;
- outer iterations;
- whether convergence is monotonic by restart.

### Stage 9C: prepare TEAM7-style non-contact electromagnetic benchmark, but not full cold-crucible yet

Only after Stage 9A/9B are stable, proceed to a simple benchmark with realistic A-phi terms. Do not jump directly to the full cold-crucible geometry.

Recommended order:

1. homogeneous box with impressed current source;
2. two-material conductor/air-like manufactured interface test;
3. simplified TEAM7-like plate/coil source without Contact;
4. only then discuss cold-crucible-specific geometry and heating coupling.

### Stage 9D: thermal coupling remains out of scope until electromagnetic solve is stable

Do not add Joule heating or temperature feedback until the electromagnetic solve has passed at least one nontrivial benchmark.

## 7. Interpretation of PCG/BiCGStab failures

Recommended wording in the progress document:

> Stage 8C shows that the operator and residual plumbing are sufficiently consistent for Krylov solution: restarted right-preconditioned GMRES converges for scalar phi Laplace + penalty, vector Helmholtz, and full A-phi + penalty manufactured systems. The remaining divergence of PCG and BiCGStab is therefore interpreted as a solver/preconditioner compatibility issue rather than an operator assembly failure. PCG is not adopted as the production solver because its convergence assumptions are too restrictive for the current matrix-free A-phi discretization and block-Jacobi preconditioning. BiCGStab remains diagnostic only. GMRES becomes the default solver for subsequent A-phi development.

## 8. Minimum acceptance before moving to TEAM7-like tests

The following should pass:

```text
test_3d_aphi_ck_gmres_scalar_phi_laplace_penalty_convergence
test_3d_aphi_ck_gmres_vector_helmholtz_convergence
test_3d_aphi_ck_gmres_full_aphi_penalty_convergence
test_3d_aphi_ck_scalar_phi_random_adjointness_diagnostic
test_3d_aphi_ck_scalar_phi_random_curvature_diagnostic
test_3d_aphi_ck_scalar_phi_tiny_host_matrix_diagnostic
test_3d_aphi_ck_block_jacobi_diagonal_diagnostic
```

The GMRES tests should check final true residual, not only recursive residual.
