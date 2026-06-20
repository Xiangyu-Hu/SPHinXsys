# OPHELIE div-grad / phi Discretization Diagnostic Archive Index

> This directory is an **index only**; it does not move original documents under the repository root, to avoid broken links.

## Archive background

Before June 2026 the primary path was `div-grad baseline` (`--phi-projection-operator=div-grad`) plus extensive phi discretization diagnostics (Neumann, compatible D/G, vector divergence MMS, etc.).

On French reduced + multiloop Biot `A_src`, `phi_eq_res_vol` stayed at **~0.49** for a long time; engineering bugs such as RHS overwrite, insufficient GMRES iterations, and simple grad corrections were ruled out as the main cause. The more likely reason is a **discrete range floor** of the **DivSigmaGrad / DivSigmaA elimination form** for non-gradient Biot fields.

**Route pivot 2026-06-02**: start **SPH pairwise edge-flux current formulation** (see repository-root `OPHELIE_EDGE_FLUX_SOLVER_DEVELOPMENT_PLAN.md`).

## Current role of each route

| Route | CLI | Status |
|------|-----|------|
| div-grad baseline | `--phi-projection-operator=div-grad` (default) | **fallback / reference / regression** |
| compatible-div-grad | `--phi-projection-operator=compatible-div-grad` | diagnostic-only |
| edge-flux phi-only | `--phi-projection-operator=edge-flux` only (no current-form) | diagnostic-only; prints warning |
| **edge-flux production** | `--ophelie-current-form=edge-flux` | **Stage 1 new primary path** |

In dev log §4.12, “freeze edge-flux” refers to the **old phi-only prototype**, not the new full edge-flux solver.

## Archived document manifest (still at repository root)

| File | Description |
|------|------|
| `OPHELIE_PHI_DEVELOPMENT_LOG.md` | Phi discretization sprint history and decisions |
| `OPHELIE_PHI_DISCRETIZATION_NEXT_SPRINT_PLAN.md` | div-grad tuning / Neumann / grad correction plan |
| `OPHELIE_DIVERGENCE_OPERATOR_VALIDATION_AND_EDGE_FLUX_PLAN.md` | Vector MMS + edge diagnostics plan (**superseded**, historical reference) |
| `OPHELIE_GUIDE_MISUNDERSTANDING_AUDIT_AND_CURSOR_IMPLEMENTATION_PLAN.md` | GPT misunderstanding audit |
| `OPHELIE_PHI_BOUNDARY_AND_NEXT_DEVELOPMENT_PLAN.md` | Boundary Neumann plan |

## Current primary-path documents

- **`OPHELIE_EDGE_FLUX_SOLVER_DEVELOPMENT_PLAN.md`** — edge-flux Stage 1–3 implementation guide
- **`tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/README.md`** — source tree layout (main / legacy / diagnostics / team7 / french_extras / stage2)

## Diagnostic tests still worth keeping

- `test_3d_ophelie_vector_divergence_mms` — D/G operator MMS (sign diagnostic, not edge production acceptance)
- `test_3d_ophelie_phi_neumann_*` — boundary MMS
- `test_3d_ophelie_phi_compatible_operator_mms` — compatible operator MMS

## Key historical conclusions (for papers / debug citations)

1. div-grad baseline: `phi_eq_res≈0.49`, `divJ_L2_red≈2.06`, 50 kW literature **PASS** (particle post-processing self-consistent).
2. edge-flux phi-only: `phi_eq_res≈1e-4`, `edge_res_red≈341`, but particle `divJ_L2_red≈0.93` — because **E/J still use particle grad**, not edge equation failure.
3. vector MMS: `linear_xyz` `sign_alpha≈-0.91` is an **MMS sign diagnostic**, not directly for edge-flux production sign.
