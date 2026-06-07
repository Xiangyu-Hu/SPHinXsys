# alternate_krylov — Stage 8 Off-Main-Line Krylov Research

This directory holds **BiCGStab** and **PCG** solver implementations. They were used in Stage 8 for comparison/diagnostics and verified **not to be the A-phi production main line** (production default is GMRES).

Unrelated to `legacy/`: `legacy/` refers to a completely different historical code path; here it is still the same A-phi matrix-free operator, only with a different Krylov backend.

## Files

| File | Description |
|---|---|
| `aphi_bicgstab_solver_ck.*` | BiCGStab + block-Jacobi PC (after Stage 8A cleanup) |
| `aphi_pcg_solver_ck.*` | PCG sanity check for scalar φ Laplace + penalty |

## Entry Points Still Referencing This Directory

- `AphiKrylovSolverKind::BiCGStabDiagnostic` / `PCGScalarDiagnostic` branches in `aphi_matrix_free_solve_ck.*`
- Several `test_3d_aphi_ck_bicgstab_*` / `test_3d_aphi_ck_pcg_*` diagnostic cases

New production code **should not** include this directory; use only when diagnostic tests or matrix-free dispatch need it.
