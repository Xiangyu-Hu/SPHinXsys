# A-phi Case Support

This folder contains case-specific and source-driving support code that is
used by A-phi tests, but is not part of the migrated SPHinXsys matrix-free
solver core itself.

Contents:
- `electromagnetic_team7_aphi_*`: Team7-like dynamics helpers.
- `electromagnetic_multiturn_coil_drive*`: multiturn coil/source support.
- `electromagnetic_team7_like_case_utils*`: Team7-like geometry and material
  helpers shared by baseline and matrix-free tests.
- `em_aphi_residual_diagnostics_snippet.hpp`: historical snippet/reference
  helper related to case diagnostics.

For archived baseline-only paths, see:
- `../legacy_aphi_archive/`
