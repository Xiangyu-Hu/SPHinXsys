# A-phi SPHinXsys Core

This folder contains only the later migrated A-phi core that follows the
SPHinXsys SYCL particle-discretization / computing-kernel organization.

Contents:
- `electromagnetic_aphi_matrix_free_*`: matrix-free Laplace-type A-phi solver,
  residual, gauge, graph, workspace, and SYCL kernel wiring.
- `electromagnetic_aphi_sph_*`: SPH-native particle operators used by the
  matrix-free solver, including assembly, gradient, divergence, Laplace,
  weighted gradient, diagnostics, and native context helpers.
- `records/`: migration logs and rewrite records for this SPHinXsys-core path.

This is the folder to hand to Cursor if you want it to explain the migrated
SPHinXsys/SYCL A-phi implementation.

For most current work, pair this folder with:
- `../aphi_case_support/`

Archived older implementations live under:
- `../legacy_aphi_archive/`
