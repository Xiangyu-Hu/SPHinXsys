# Matrix-Free Gauge Projection Test

## Purpose

This is the first relation-based runnable gauge-projection verification case for the restarted Laplace A-phi matrix-free line.

## Current Version

The current version uses:

- a thin box `SolidBody`
- lattice particles
- a real `InnerRelation`
- a `MatrixFreePairwiseGraph` extracted from `MatrixFreeAPhiDiscreteView`

The default geometry is intentionally thin in the transverse directions so the baseline stays quasi-1D while still using a real SPH relation.

## What It Verifies

- scalar `chi` solve for the gauge step
- gauge transform update of `A_x` and `phi`
- reduction of the operator-consistent projection residual
- reduction of the divergence-like diagnostic built from the same pairwise graph
- approximate invariance of `E_x = -i * omega * A_x - grad(phi)_x`

## Current Outputs

- `em_aphi_matrix_free_gauge_projection_summary_<tag>.csv`
- `em_aphi_matrix_free_gauge_projection_profile_<tag>.csv`
- stdout summary line with residual and diagnostic statistics

## Main Environment Variables

- `EM_APHI_MATRIX_FREE_GAUGE_DP`
- `EM_APHI_MATRIX_FREE_GAUGE_LENGTH`
- `EM_APHI_MATRIX_FREE_GAUGE_HEIGHT`
- `EM_APHI_MATRIX_FREE_GAUGE_WIDTH`
- `EM_APHI_MATRIX_FREE_GAUGE_OMEGA`
- `EM_APHI_MATRIX_FREE_GAUGE_CHI_REAL`
- `EM_APHI_MATRIX_FREE_GAUGE_CHI_IMAG`
- `EM_APHI_MATRIX_FREE_GAUGE_AX0_REAL`
- `EM_APHI_MATRIX_FREE_GAUGE_AX0_IMAG`
- `EM_APHI_MATRIX_FREE_GAUGE_MAX_ITERATIONS`
- `EM_APHI_MATRIX_FREE_GAUGE_RELAXATION`
- `EM_APHI_MATRIX_FREE_GAUGE_ABS_TOL`
- `EM_APHI_MATRIX_FREE_GAUGE_OUTPUT_TAG`
