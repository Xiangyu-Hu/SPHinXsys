# Matrix-Free Complex Helmholtz Test

## Purpose

This is the first relation-based runnable verification case for the restarted Laplace A-phi matrix-free line.
It isolates a single complex scalar field before we enter the full coupled `A-phi` solve.

## Current Version

The current version uses:

- a simple box `SolidBody`
- lattice particles
- a real `InnerRelation`
- a `MatrixFreePairwiseGraph` extracted from `MatrixFreeAPhiDiscreteView`

The default geometry is intentionally thin in the transverse directions so the baseline behaves like a quasi-1D manufactured verification problem while still using a real SPH relation.

## Target Equation

We solve

`-Laplace(u) + alpha * u = f`

with a manufactured complex solution.

## Why This Case Matters

- it is the first matrix-free scalar verification that uses a real SPH body and inner relation
- it validates the shared pairwise-graph helper path
- it bridges the gap between 1D synthetic prototypes and later full relation-based `A-phi` assembly

## Current Outputs

- `em_aphi_matrix_free_complex_helmholtz_summary_<tag>.csv`
- `em_aphi_matrix_free_complex_helmholtz_profile_<tag>.csv`
- stdout summary line with residual and error statistics

## Main Environment Variables

- `EM_APHI_MATRIX_FREE_HELMHOLTZ_DP`
- `EM_APHI_MATRIX_FREE_HELMHOLTZ_LENGTH`
- `EM_APHI_MATRIX_FREE_HELMHOLTZ_HEIGHT`
- `EM_APHI_MATRIX_FREE_HELMHOLTZ_WIDTH`
- `EM_APHI_MATRIX_FREE_HELMHOLTZ_ALPHA_REAL`
- `EM_APHI_MATRIX_FREE_HELMHOLTZ_ALPHA_IMAG`
- `EM_APHI_MATRIX_FREE_HELMHOLTZ_RHS_MODE`
- `EM_APHI_MATRIX_FREE_HELMHOLTZ_MAX_ITERATIONS`
- `EM_APHI_MATRIX_FREE_HELMHOLTZ_RELAXATION`
- `EM_APHI_MATRIX_FREE_HELMHOLTZ_ABS_TOL`
- `EM_APHI_MATRIX_FREE_HELMHOLTZ_OUTPUT_TAG`
