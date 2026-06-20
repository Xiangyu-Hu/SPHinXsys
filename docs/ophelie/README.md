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

## Client-facing report

- LaTeX/PDF: [`RH200_OPHELIE_CLIENT_PROGRESS_REPORT.tex`](../../RH200_OPHELIE_CLIENT_PROGRESS_REPORT.tex) (repo root)

## Archived originals

| Archive folder | Contents |
|----------------|----------|
| [`archive/plans/`](archive/plans/) | Former repo-root `OPHELIE_*.md`, `RH200_*.md`, `FRENCH_*.md` (30 files) |
| [`archive/discussion/`](archive/discussion/) | Former `discussion_bundle/*.md` markdown |
| [`archive/discussion_packs/`](archive/discussion_packs/) | TEAM7 P5 / P56 lite pack markdown copies |
| [`archive/legacy_packages/`](archive/legacy_packages/) | Former GPT discussion packages (now English) |

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
