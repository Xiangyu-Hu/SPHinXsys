# OPHELIE / RH200 Documentation Index

> **English-only documentation** (except client LaTeX report at repo root).  
> Full archived plans and discussion records under `archive/` — nothing was deleted.

## Start here

| Document | Scope |
|----------|--------|
| [01_MASTER_DEVELOPMENT_TIMELINE.md](01_MASTER_DEVELOPMENT_TIMELINE.md) | Cross-cutting timeline, milestones, current production status |
| [02_PHI_OPERATOR_AND_BOUNDARY.md](02_PHI_OPERATOR_AND_BOUNDARY.md) | Scalar potential φ, operators, boundaries, MMS, GMRES |
| [03_EDGE_FLUX_PRODUCTION.md](03_EDGE_FLUX_PRODUCTION.md) | Edge-flux current equation, power calibration, complex mode |
| [04_FRENCH_REDUCED_AND_THERMAL.md](04_FRENCH_REDUCED_AND_THERMAL.md) | French reduced geometry, literature mode, thermal one-way |
| [05_TEAM7_BENCHMARK.md](05_TEAM7_BENCHMARK.md) | TEAM7 side-track: L2 validation, P5 no-flux, high-σ scaling |
| [06_RH200_GLASS_EM_STIRRING.md](06_RH200_GLASS_EM_STIRRING.md) | RH200 midterm case: EM → Joule → Euler grid → SPH stirring |
| [07_GPT_DISCUSSION_ARCHIVE.md](07_GPT_DISCUSSION_ARCHIVE.md) | Merged GPT discussion summaries (edge-flux, TEAM7 rounds) |
| [08_EDGE_FLUX_UNDIRECTED_POWER_CLOSURE.md](08_EDGE_FLUX_UNDIRECTED_POWER_CLOSURE.md) | Stage 1 B1 reopen: volume-consistent / undirected edge power (advisor pack) |
| [09_FRENCH_LITERATURE_REPRODUCTION_STATUS_2026-09.md](09_FRENCH_LITERATURE_REPRODUCTION_STATUS_2026-09.md) | 2026-09 status vs Jacoutot CES 2008 + roadmap options for advisor/ChatGPT |
| [10_THERMO_EM_COUPLING_ROUND1.md](10_THERMO_EM_COUPLING_ROUND1.md) | Round-1 periodic thermo–EM coupling (natural first); audit + CLI |
| [11_THERMO_EM_COUPLING_AFTER_CHATGPT_TASK.md](11_THERMO_EM_COUPLING_AFTER_CHATGPT_TASK.md) | 2026-09-14 close-out: hydro fix, natural 50 s, stirring 5 rev, spatial stats |
| [12_THERMO_EM_COUPLING_CHATGPT_DISCUSSION_2026-09.md](12_THERMO_EM_COUPLING_CHATGPT_DISCUSSION_2026-09.md) | Paste-to-ChatGPT pack: post-brief status + open questions |
| [`archive/DOCS_001_050_MERGED.md`](archive/DOCS_001_050_MERGED.md) | Docs 01–50 combined archive (Doc13–50 no longer standalone) |
| [`archive/DOCS_051_100_MERGED.md`](archive/DOCS_051_100_MERGED.md) | Docs 51–100 combined archive |
| [`archive/DOCS_101_150_MERGED.md`](archive/DOCS_101_150_MERGED.md) | Docs 101–150 combined archive |
| [`archive/DOCS_151_200_MERGED.md`](archive/DOCS_151_200_MERGED.md) | Docs 151–200 combined archive |
| [`archive/README_MERGED_DOCS.md`](archive/README_MERGED_DOCS.md) | Merge policy (201+ stay standalone) |
| [packs/gpt_doc76_to_83_discussion_2026-10-03.tar.gz](packs/gpt_doc76_to_83_discussion_2026-10-03.tar.gz) | Pack: Docs76–83 + A/B/C/French code + results |
| [OPHELIE_CURRENT_CODE_AUDIT.md](OPHELIE_CURRENT_CODE_AUDIT.md) | 2026-08 code audit vs Jacoutot papers |
| [OPHELIE_STAGE3_4_CLOSEOUT_AND_FRENCH_NATURAL_NEXT_PLAN_2026-08-11.md](OPHELIE_STAGE3_4_CLOSEOUT_AND_FRENCH_NATURAL_NEXT_PLAN_2026-08-11.md) | Binding plan: φ dual gate, σ_nat, Stage 3.5 → Natural thermal / Stage 4.x |

## French Natural track (Stage 3.5–4.1)

| Test | Role |
|------|------|
| `test_3d_ophelie_french_natural_glass_relax` | Reload mesh for Natural geometry |
| `test_3d_ophelie_french_natural_em` | Literature EM + power gate |
| `test_3d_ophelie_french_complex_joule_to_heat_one_way` | σ(T), thermal BC, Stage 3.x |
| `test_3d_ophelie_french_natural_convection_frozen_q` | Stage 4.1 frozen-Q + WCSPH + Boussinesq |
| `test_3d_ophelie_boussinesq_stratified_hydrostatic` | Doc77-A: closed-box stratified hydrostatic (production Boussinesq) |
| `test_3d_ophelie_thermal_1d_conduction` | Doc77-B: 1D conduction (ΔT Laplace path) |
| `test_3d_ophelie_devahl_davis_cavity` | Doc77-C: de Vahl Davis cavity (freeze-T smoke; Nu refine open) |
| `test_3d_ophelie_viscosity_reference` | Doc180: SYCL/CK ViscousForceCK vs production InteractionSplit vs analytic |
| `test_3d_ophelie_french_natural_reduced_pressure` | French reduced-pressure production path (Doc126+); Doc180 full EM–thermal–flow; Doc183 field export |

| Doc | Role |
|-----|------|
| [`reference/doc180_runs/`](reference/doc180_runs/) | Results: P1/P2 refs, coupled cold dt005/dt0025, readiness |
| [`reference/doc183_runs/`](reference/doc183_runs/) | Field exports, diagnostics figures, A/B/C compares |
| [`reference/doc186_runs/`](reference/doc186_runs/) | Classifier tests, U0/U1/U2, bulk/flux reviews |
| [`reference/doc200_runs/`](reference/doc200_runs/) | Doc200 results + pack under `PACKAGES/` |
| [201_OPHELIE_UNIFIED_COMPLEX_EM_AND_LIVE_BUOYANCY_CURSOR_PLAN_2026-10-30.md](201_OPHELIE_UNIFIED_COMPLEX_EM_AND_LIVE_BUOYANCY_CURSOR_PLAN_2026-10-30.md) | Doc201: unified complex EM + live buoyancy / acoustic top |
| [202_DOC201_ROUND_SUMMARY_AND_CLOSE_2026-10-31.md](202_DOC201_ROUND_SUMMARY_AND_CLOSE_2026-10-31.md) | Doc201 closeout; §13 answers |
| [`reference/doc201_runs/`](reference/doc201_runs/) | Single-layer evidence; [`DELIVERABLE_CHECKLIST.md`](reference/doc201_runs/DELIVERABLE_CHECKLIST.md); brief + NEXT |
| [202_OPHELIE_B2_LID_SIDE_SPLIT_ACOUSTIC_CURSOR_PLAN_2026-10-31.md](202_OPHELIE_B2_LID_SIDE_SPLIT_ACOUSTIC_CURSOR_PLAN_2026-10-31.md) | Doc202: B2 lid/side acoustic split (non-dens) |
| [`reference/doc202_runs/`](reference/doc202_runs/) | H1/H1b split reviews; [`H1_SPLIT_REVIEW.md`](reference/doc202_runs/H1_SPLIT_REVIEW.md) |
| [203_DOC201_DOC202_ROUND_PACK_2026-10-31.md](203_DOC201_DOC202_ROUND_PACK_2026-10-31.md) | Round pack pointer; `packs/Doc201_Doc202_UNIFIED_EM_LIVE_B2_2026-10-31.tar.gz` |
| [212_OPHELIE_TARGETED_HYDRO_EM_REPAIR_CURSOR_PLAN_2026-11-06.md](212_OPHELIE_TARGETED_HYDRO_EM_REPAIR_CURSOR_PLAN_2026-11-06.md) | Doc212: cavity CFL + shared complex EMF; phase/edge PASS; field PARTIAL |
| [213_OPHELIE_LEGACYFLUX_EMF_SKIN_CLOSURE_CURSOR_PLAN_2026-11-07.md](213_OPHELIE_LEGACYFLUX_EMF_SKIN_CLOSURE_CURSOR_PLAN_2026-11-07.md) | Doc213: EMF-internal skin/field closure (current) |
| [`reference/doc213_runs/`](reference/doc213_runs/) | G0 baselines / EMF surgery evidence（本地 runs） |

**Geometry freeze:** glass Ø **500 mm** (thesis; user-confirmed 2026-09-30). Crucible Ø650 mm ≠ melt diameter.  
σ(T) tables: [`reference/jacoutot_sigma_t/`](reference/jacoutot_sigma_t/)

## Client-facing report

- LaTeX/PDF: [`RH200_OPHELIE_CLIENT_PROGRESS_REPORT.tex`](../../RH200_OPHELIE_CLIENT_PROGRESS_REPORT.tex) (repo root)

## Archived originals

| Archive folder | Contents |
|----------------|----------|
| [`archive/plans/`](archive/plans/) | Former repo-root `OPHELIE_*.md`, `RH200_*.md`, `FRENCH_*.md` (30 files) |
| [`archive/discussion/`](archive/discussion/) | Stage work records, GPT reviews, discussion markdown |
| [`archive/discussion_packs/`](archive/discussion_packs/) | TEAM7 P5 / P56 lite pack markdown copies |
| [`archive/legacy_packages/`](archive/legacy_packages/) | Former GPT discussion packages (now English) |

> Local ephemeral GPT upload packs (`docs/ophelie/packs/`) are **not** part of the repository.

## TEAM7 reference data

| Path | Content |
|------|---------|
| [`tests/.../reference_data/team7/`](../../tests/extra_source_and_tests/3d_examples/reference_data/team7/) | Probe definitions, COMSOL/muFEM reference CSVs |
| [`docs/TEAM7-reference/README.md`](../TEAM7-reference/README.md) | Pointer only (no source code) |

## Data / discussion text

- `discussion_bundle/` — GPT prompts and French/TEAM7 summary `.txt` only (no code, no run CSVs)

## Source code index

- EM core: `tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/README.md`
- RH200 case: `tests/extra_source_and_tests/3d_examples/test_3d_ophelie_rh200_glass_em_stirring/README.md`
- TEAM7 case: `tests/extra_source_and_tests/3d_examples/test_3d_ophelie_team7_complex_edge_flux/TEAM7_VALIDATION_LOG.md`
