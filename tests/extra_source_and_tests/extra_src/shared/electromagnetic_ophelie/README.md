# electromagnetic_ophelie source layout

> Edge-flux Stage 1 production path + shared infrastructure stay in **this directory root**; legacy routes, diagnostics, Team7, and French extensions live in subdirectories by category.
>
> CMake (`extra_src/CMakeLists.txt`) recursively scans all subdirectories under `shared/` and adds them to the include path, so `#include "electromagnetic_ophelie_*.h"` **still uses flat filenames**—no subdirectory prefix in includes.

## Root directory (edge-flux + required shared code)

| File | Purpose |
|------|------|
| `electromagnetic_ophelie.h` | Umbrella header |
| `electromagnetic_ophelie_edge_flux.h` / `.hpp` | **Edge-flux production**: C_ij, e_ij, Residual/JouleHeat/Recon CK |
| `electromagnetic_ophelie_edge_flux_diagnostics.h` | Edge acceptance diagnostics (legacy + production) |
| `electromagnetic_ophelie_phi.h` / `.hpp` | φ solve (div-grad fallback + edge-flux RHS) |
| `electromagnetic_ophelie_phi_gmres.h` | GMRES solver |
| `electromagnetic_ophelie_phi_krylov_ck.h` | Krylov device vector ops |
| `electromagnetic_ophelie_phi_device_vector_ops.h` | Device vector ops test helpers |
| `electromagnetic_ophelie_phi_solvability.h` | RHS solvability diagnostics (French pipeline) |
| `electromagnetic_ophelie_phi_boundary.h` | φ boundary conditions (production) |
| `electromagnetic_ophelie_laplace.h` | Laplace / div-grad operator |
| `electromagnetic_ophelie_french_literature.h` | French reduced literature acceptance pipeline |
| `electromagnetic_ophelie_french_reduced_geometry.h` | French reduced geometry |
| `electromagnetic_ophelie_multiloop_source.h` | Multiloop coil source |
| `electromagnetic_ophelie_current_moment_source.h` | Current-moment source (multiloop dependency) |
| `electromagnetic_ophelie_biot_savart.h` / `.hpp` | Biot–Savart A_src |
| `electromagnetic_ophelie_source.h` / `.hpp` | Generic source terms |
| `electromagnetic_ophelie_postprocess.h` / `.hpp` | E/J/Joule post-processing |
| `electromagnetic_ophelie_cli.h` | CLI (`--ophelie-current-form`, etc.) |
| `electromagnetic_ophelie_parameters.h` | Parameters and acceptance thresholds |
| `electromagnetic_ophelie_field_names.h` | Field name constants |
| `electromagnetic_ophelie_register_fields.h` | Field registration |
| `electromagnetic_ophelie_device_sync.h` | Device sync helpers |
| `electromagnetic_ophelie_observables.h` | Observables |
| `electromagnetic_ophelie_diagnostics.h` | Generic divJ / power diagnostics |

## `legacy/div_grad/` — div-grad baseline fallback

| File | Purpose |
|------|------|
| `electromagnetic_ophelie_phi_gradient.h` | Particle-gradient div-grad path (`--phi-projection-operator=div-grad`) |

## `diagnostics/` — MMS / operator audit / non-production diagnostics

| File | Purpose |
|------|------|
| `electromagnetic_ophelie_phi_mms_helpers.h` | φ MMS manufactured-solution helpers |
| `electromagnetic_ophelie_phi_operator_diagnostics.h` | compatible-div-grad and other operator linear consistency |
| `electromagnetic_ophelie_phi_boundary_diagnostics.h` | Neumann boundary MMS |
| `electromagnetic_ophelie_phi_rhs_diagnostics.h` | RHS flux sign audit |
| `electromagnetic_ophelie_vector_divergence_diagnostics.h` | Vector D/G MMS (sign diagnostics) |
| `electromagnetic_ophelie_aind_diagnostic.h` | A_ind diagnostics (Stage 2 precursor) |

## `team7/` — Team7 / racetrack case-specific

| File | Purpose |
|------|------|
| `electromagnetic_ophelie_team7_geometry.h` | Team7 geometry |
| `electromagnetic_ophelie_team7_native_geometry.h` | Team7 native geometry |
| `electromagnetic_ophelie_team7_probe.h` | Team7 probes |
| `electromagnetic_ophelie_racetrack_source.h` | Racetrack source (also referenced from CLI) |

## `french_extras/` — French glass relax and related case extensions

| File | Purpose |
|------|------|
| `electromagnetic_ophelie_french_glass_mesh_relax.h` | Glass mesh relax |
| `electromagnetic_ophelie_relaxation.h` | Particle relaxation |
| `electromagnetic_ophelie_progress.h` | Progress output |

## `stage2/` — Self-induction / coupling (not yet production-gated)

| File | Purpose |
|------|------|
| `electromagnetic_ophelie_self_induction.h` | Self-induction Picard iteration |

## CLI route map

| Route | CLI |
|------|-----|
| **Edge-flux production** | `--ophelie-current-form=edge-flux` |
| div-grad baseline | `--ophelie-current-form=particle-gradient --phi-projection-operator=div-grad` |
| edge-flux phi-only diagnostic | `--ophelie-current-form=particle-gradient --phi-projection-operator=edge-flux` |

## Related documentation

- Archive index: `docs/ophelie/archive_phi_divgrad_residual_floor/README.md`
- Edge-flux plan: `OPHELIE_EDGE_FLUX_SOLVER_DEVELOPMENT_PLAN.md`
- Run/discussion bundle: `discussion_bundle/`
