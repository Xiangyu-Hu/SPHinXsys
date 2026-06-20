# TEAM7 L2 — ChatGPT upload file manifest (Round 3)

> **Primary discussion document (required; upload first)**: 
> [`OPHELIE_GPT_DISCUSSION_TEAM7_L2_EDGE_FLUX_ROUND3.md`](OPHELIE_GPT_DISCUSSION_TEAM7_L2_EDGE_FLUX_ROUND3.md)

> **One-click zip (recommended)**: 
> `discussion_bundle/team7_l2_gpt_upload_round3.zip`

> **Round 2 background (optional)**: 
> [`OPHELIE_GPT_DISCUSSION_TEAM7_L2_EDGE_FLUX.md`](OPHELIE_GPT_DISCUSSION_TEAM7_L2_EDGE_FLUX.md)

---

## A. Required — discussion and milestones (6)

| # | Absolute path | Description |
|---|----------|------|
| 1 | `discussion_bundle/OPHELIE_GPT_DISCUSSION_TEAM7_L2_EDGE_FLUX_ROUND3.md` | **Round 3 main document**: P0–P4, 12 GPT issues |
| 2 | `discussion_bundle/OPHELIE_GPT_DISCUSSION_TEAM7_L2_EDGE_FLUX.md` | Round 2 main document (background) |
| 3 | `tests/.../TEAM7_VALIDATION_LOG.md` | dedicated validation milestones (P0–P4) |
| 4 | `OPHELIE_TEAM7_HIGH_SIGMA_EDGE_FLUX_BOTTLENECK_AND_NEXT_PLAN.md` | GPT push master plan §12 P4 |
| 5 | `tests/.../HIGH_SIGMA_EDGE_FLUX_P1A_LOG.md` | P1a detailed log |
| 6 | `discussion_bundle/TEAM7_L2_GPT_FIRST_PROMPT_ROUND3.txt` | first prompt |

---

## B. Required — Round 3 run CSV / logs

### B.1 P0 normalization

| # | File | Description |
|---|------|------|
| 7 | `team7_l2_outputs/team7_normalization_sweep_one-way_f50_fil_*.csv` | safe_rhs_l2 sweep |
| 8 | `team7_l2_outputs/logs/team7_run_normalization_sweep.log` | full stdout |

### B.2 P1a high-σ benchmark

| # | File | Description |
|---|------|------|
| 9 | `team7_l2_outputs/high_sigma_edge_flux_scaling_P1a_v3.csv` | three cases × σ × f |
| 10 | `team7_l2_outputs/logs/p1a_high_sigma_sweep_v6.log` | final v6 log |

### B.3 P2β / P2 Picard (solver-local)

| # | File | Description |
|---|------|------|
| 11 | `team7_l2_outputs/logs/p2_picard_solver_local_r005.log` | physical-scale Picard 6 iters |

### B.4 P4 operator audit

| # | File | Description |
|---|------|------|
| 12 | `team7_l2_outputs/team7_edge_conductance_audit_one-way_p4_audit.csv` | P4.1 C_ij stats |
| 13 | `team7_l2_outputs/team7_edge_partition_audit_one-way_p4_audit.csv` | P4.4 partition J/B/e_edge |
| 14 | `team7_l2_outputs/team7_probe_one-way_p4_audit.csv` | same-run probe (if copied) |
| 15 | `team7_l2_outputs/logs/p4_operator_audit_oneway.log` | full stdout |

### B.5 Round 2 baseline (context; recommended upload)

| # | File | Description |
|---|------|------|
| 16 | `team7_l2_outputs/TEAM7_Bz_A1_B1_reference_mT.csv` | Bz reference |
| 17 | `team7_l2_outputs/TEAM7_Jey_A1_surface_reference_Am2.csv` | Jy reference |
| 18 | `team7_l2_outputs/team7_probe_one-way_f50_fil_*.csv` | 50 Hz probe decomposition |
| 19 | `team7_l2_outputs/team7_omega_scaling_one-way_f50_fil_*.csv` | P3c ω sweep |
| 20 | `team7_l2_outputs/team7_j_vs_b_probe_split_one-way_f50_fil_*.csv` | J/B split |
| 21 | `team7_l2_outputs/team7_picard_r0.05.csv` etc. | historical field-scale Picard |
| 22 | `team7_l2_outputs/logs/team7_run_filament_one_way_src1.log` | Round 2 baseline log |
| 23 | `team7_l2_outputs/logs/team7_run_omega_sweep.log` | ω sweep log |

---

## C. Required — Round 3 source (18)

| # | Absolute path | Description |
|---|----------|------|
| 24 | `.../diagnostics/electromagnetic_ophelie_edge_flux_operator_audit.h` | **P4 new diagnostic** |
| 25 | `.../test_3d_ophelie_team7_complex_edge_flux.cpp` | TEAM7 test + `--team7-operator-audit` |
| 26 | `.../test_3d_ophelie_high_sigma_edge_flux_scaling.cpp` | **P1a new test** |
| 27 | `.../test_3d_ophelie_high_sigma_edge_flux_scaling/CMakeLists.txt` | P1a build |
| 28 | `.../team7/electromagnetic_ophelie_team7_validation.h` | Picard, normalization sweep |
| 29 | `.../diagnostics/electromagnetic_ophelie_aind_diagnostic.h` | one-way + restore audit |
| 30 | `.../electromagnetic_ophelie_edge_flux.h` | edge-flux + solver-local recon |
| 31 | `.../electromagnetic_ophelie_parameters.h` | `edge_flux_solver_local_rhs_scale_` |
| 32 | `.../electromagnetic_ophelie_phi_component.h` | φ PCG solver-local scaling |
| 33 | `.../diagnostics/electromagnetic_ophelie_phi_rhs_diagnostics.h` | scale helpers |
| 34 | `.../electromagnetic_ophelie_phi_gmres.h` | GMRES path |
| 35 | `.../electromagnetic_ophelie_cli.h` | normalization / operator-audit CLI |
| 36 | `.../stage2/electromagnetic_ophelie_self_induction.h` | Picard A_ind relax |
| 37 | `.../team7/electromagnetic_ophelie_team7_probe.h` | probe reference |
| 38 | `.../team7/electromagnetic_ophelie_team7_coil_path_source.h` | filament source |
| 39 | `.../diagnostics/aphi_team7_native_geometry_config.h` | TEAM7 geometry |
| 40 | `.../particle_generation_TEAM7/particle_generation_team7.cpp` | particle generation |
| 41 | `.../test_3d_ophelie_team7_complex_edge_flux/README.md` | CLI description |

(`...` = corresponding path under `/home/yyc/SPHinXsysSYCL/tests/extra_source_and_tests/`)

---

## D. One-click packaging command

```bash
cd /home/yyc/SPHinXsysSYCL

cp -f build/output/team7_probe_one-way_p4_audit.csv \
 discussion_bundle/team7_l2_outputs/ 2>/dev/null || true

zip -r discussion_bundle/team7_l2_gpt_upload_round3.zip \
 docs/ophelie/archive/discussion/OPHELIE_GPT_DISCUSSION_TEAM7_L2_EDGE_FLUX_ROUND3.md \
 docs/ophelie/archive/discussion/OPHELIE_GPT_DISCUSSION_TEAM7_L2_EDGE_FLUX.md \
 docs/ophelie/archive/discussion/TEAM7_L2_UPLOAD_MANIFEST_ROUND3.md \
 discussion_bundle/TEAM7_L2_GPT_FIRST_PROMPT_ROUND3.txt \
 docs/ophelie/archive/discussion/TEAM7_L2_UPLOAD_MANIFEST.md \
 discussion_bundle/TEAM7_L2_GPT_FIRST_PROMPT_ROUND2.txt \
 discussion_bundle/README.md \
 discussion_bundle/team7_l2_outputs/ \
 OPHELIE_TEAM7_HIGH_SIGMA_EDGE_FLUX_BOTTLENECK_AND_NEXT_PLAN.md \
 tests/extra_source_and_tests/3d_examples/test_3d_ophelie_team7_complex_edge_flux/TEAM7_VALIDATION_LOG.md \
 tests/extra_source_and_tests/3d_examples/test_3d_ophelie_team7_complex_edge_flux/README.md \
 tests/extra_source_and_tests/3d_examples/test_3d_ophelie_team7_complex_edge_flux/test_3d_ophelie_team7_complex_edge_flux.cpp \
 tests/extra_source_and_tests/3d_examples/test_3d_ophelie_high_sigma_edge_flux_scaling/test_3d_ophelie_high_sigma_edge_flux_scaling.cpp \
 tests/extra_source_and_tests/3d_examples/test_3d_ophelie_high_sigma_edge_flux_scaling/CMakeLists.txt \
 tests/extra_source_and_tests/3d_examples/test_3d_ophelie_high_sigma_edge_flux_scaling/HIGH_SIGMA_EDGE_FLUX_P1A_LOG.md \
 tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/diagnostics/electromagnetic_ophelie_edge_flux_operator_audit.h \
 tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/team7/electromagnetic_ophelie_team7_validation.h \
 tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/team7/electromagnetic_ophelie_team7_probe.h \
 tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/team7/electromagnetic_ophelie_team7_coil_path_source.h \
 tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/diagnostics/electromagnetic_ophelie_aind_diagnostic.h \
 tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/diagnostics/aphi_team7_native_geometry_config.h \
 tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_edge_flux.h \
 tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_parameters.h \
 tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_phi_component.h \
 tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/diagnostics/electromagnetic_ophelie_phi_rhs_diagnostics.h \
 tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_phi_gmres.h \
 tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_cli.h \
 tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/stage2/electromagnetic_ophelie_self_induction.h \
 tests/extra_source_and_tests/3d_examples/particle_generation_TEAM7/particle_generation_team7.cpp

ls -lh discussion_bundle/team7_l2_gpt_upload_round3.zip
```

---

## E. First prompt for ChatGPT

See [`TEAM7_L2_GPT_FIRST_PROMPT_ROUND3.txt`](TEAM7_L2_GPT_FIRST_PROMPT_ROUND3.txt)

---

*Manifest version: Round 3, 2026-06-11*
