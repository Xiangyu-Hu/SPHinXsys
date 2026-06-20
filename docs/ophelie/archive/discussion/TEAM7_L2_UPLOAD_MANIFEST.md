# TEAM7 L2 — ChatGPT upload file manifest (Round 2)

> **Primary discussion document (required; upload first)**: 
> [`/home/yyc/SPHinXsysSYCL/docs/ophelie/archive/discussion/OPHELIE_GPT_DISCUSSION_TEAM7_L2_EDGE_FLUX.md`](OPHELIE_GPT_DISCUSSION_TEAM7_L2_EDGE_FLUX.md)

> **One-click zip (recommended)**: 
> `/tmp/team7_l2_gpt_upload_round2.zip` (see §D)

---

## A. Required — discussion and milestones (5)

| # | Absolute path | Description |
|---|----------|------|
| 0 | `.../OPHELIE_GPT_DISCUSSION_TEAM7_ROUND3_CLOSURE.md` | **Round 3 closure + GPT decision request (2026-06-12, upload first)** |
| 1 | `.../OPHELIE_GPT_DISCUSSION_TEAM7_L2_EDGE_FLUX.md` | Round 2 historical discussion pack |
| 2 | `.../test_3d_ophelie_team7_complex_edge_flux/TEAM7_VALIDATION_LOG.md` | Dedicated validation milestone log |
| 3 | `OPHELIE_TEAM7_L2_VALIDATION_DEBUG_AND_CURSOR_PLAN.md` | Cursor/GPT execution plan |
| 4 | `TEAM7_L2_GPT_FIRST_PROMPT_ROUND2.txt` | First prompt (copy-paste) |

---

## B. Required — run CSV / logs (`discussion_bundle/team7_l2_outputs/`)

### B.1 Reference data

| # | File | Description |
|---|------|------|
| 5 | `TEAM7_Bz_A1_B1_reference_mT.csv` | Bz reference (17 probes) |
| 6 | `TEAM7_Jey_A1_surface_reference_Am2.csv` | Jy reference (50 Hz) |
| 7 | `TEAM7_probe_definitions.csv` | Probe definitions |

### B.2 Filament 50 Hz baseline (P1/P3b)

| # | File | Description |
|---|------|------|
| 8 | `team7_probe_one-way_f50_fil_asign_p1_p90_minus_src1.000_jraw.csv` | Probe decomposition |
| 9 | `team7_bz_one-way_f50_fil_*_coil_only.csv` / `*_total.csv` | phase0 comparison |
| 10 | `team7_jey_probe_one-way_f50_fil_*.csv` | Jy probe |
| 11 | `team7_j_vs_b_probe_split_one-way_f50_fil_*.csv` | **J/B split** |
| 12 | `team7_omega_scaling_one-way_f50_fil_*.csv` | **P3c ω sweep** |
| 13 | `team7_probe_phase90_xbands_one-way_f50_fil_*.csv` | x-band phase90 |
| 14 | `team7_plate_depth_profile_one-way_f50_fil_*.csv` | z-layered J/B |
| 15 | `team7_validation_history_one-way_f50_fil_*.txt` | append log |
| 16 | `team7_summary_one-way_f50_fil_*.txt` | single-line summary |

### B.3 P5 feedback / P3 200 Hz / P2b Picard

| # | File | Description |
|---|------|------|
| 17 | `team7_probe_feedback1.csv` | feedback=1 probe |
| 18 | `team7_probe_one-way_f200_fil_*.csv` | 200 Hz probe |
| 19 | `team7_picard_r0.05.csv` / `r0.15` / `r0.35` | Picard per iteration |
| 20 | `team7_j_vs_b_probe_split_picard_f50_fil_*.csv` | Picard J/B |

### B.4 Historical volume / a_sign experiments (context)

| # | File | Description |
|---|------|------|
| 21 | `team7_probe_one-way.csv` | volume baseline `a_sign=-1` |
| 22 | `team7_probe_one-way_jscale_auto.csv` | `j_post_scale=auto` |
| 23 | `team7_probe_phase90_xbands_one-way.csv` | x-band history |
| 24 | `team7_plate_depth_profile_one-way.csv` | depth profile history |

### B.5 Full stdout logs (`logs/`)

| # | File | Description |
|---|------|------|
| 25 | `logs/team7_run_filament_one_way_src1.log` | filament 50 Hz one-way |
| 26 | `logs/team7_run_omega_sweep.log` | P3c sweep |
| 27 | `logs/team7_run_feedback1.log` | P5 feedback |
| 28 | `logs/team7_run_f200_one_way.log` | 200 Hz |
| 29 | `logs/team7_run_picard_r0.05.log` | Picard relax=0.05 |
| 30 | `logs/team7_run_picard_r0.15.log` | Picard relax=0.15 |
| 31 | `logs/team7_run_picard_r0.35.log` | Picard relax=0.35 |
| 32 | `team7_run_baseline_a_sign_plus1.log` | historical `a_sign=+1` |
| 33 | `team7_run_imag_a_sign_minus1.log` | historical `a_sign=-1` |
| 34 | `team7_run_a_sign_m1_jscale_auto.log` | historical auto jscale |

### B.6 Round 3 push (2026-06-12)

| # | File | Description |
|---|------|------|
| 35 | `high_sigma_edge_flux_scaling_P1a_v4_lenz.csv` | P1a + P4 Lenz columns |
| 36 | `team7_aind_lenz_audit_one-way_f50_fil_*.csv` | TEAM7 P4 single row |
| 37 | `team7_edge_partition_audit_one-way_f50_fil_*.csv` | P2 partition + Jn/Jt + **P4.1 moment** |
| 38 | `team7_edge_conductance_audit_one-way_f50_fil_*.csv` | P4.1 global moment eigenvalues |
| 39 | `HIGH_SIGMA_EDGE_FLUX_P1A_LOG.md` | P1a detailed log |
| 40 | `OPHELIE_TEAM7_ROUND3_GAUGE_BOUNDARY_NEXT_CURSOR_PLAN.md` | Round 3 master plan |

---

## C. Recommended upload — core source (12)

| # | Absolute path | Description |
|---|----------|------|
| 35 | `.../team7/electromagnetic_ophelie_team7_validation.h` | validation, ω sweep, Picard |
| 36 | `.../team7/electromagnetic_ophelie_team7_probe.h` | Bz/Jey reference |
| 37 | `.../team7/electromagnetic_ophelie_team7_coil_path_source.h` | filament source |
| 38 | `.../diagnostics/electromagnetic_ophelie_aind_diagnostic.h` | one-way + restore |
| 39 | `.../electromagnetic_ophelie_edge_flux.h` | edge-flux kernel |
| 40 | `.../electromagnetic_ophelie_parameters.h` | parameters |
| 41 | `.../electromagnetic_ophelie_cli.h` | CLI |
| 42 | `.../diagnostics/aphi_team7_native_geometry_config.h` | TEAM7 geometry |
| 43 | `.../test_3d_ophelie_team7_complex_edge_flux.cpp` | test entry |
| 44 | `.../test_3d_ophelie_team7_complex_edge_flux/README.md` | CLI description |
| 45 | `.../reference_data/team7/tools/import_team7_jey_comsol_table2.py` | Jey 200 Hz import |
| 46 | `.../particle_generation_TEAM7/particle_generation_team7.cpp` | particle generation |

(`...` = corresponding path under `/home/yyc/SPHinXsysSYCL/tests/extra_source_and_tests/`)

---

## D. Optional — background documents

| # | File | Description |
|---|------|------|
| 47 | `OPHELIE_PHI_DEVELOPMENT_LOG.md` | master dev log |
| 48 | `OPHELIE_GPT_DISCUSSION_SINCE_COMPLEX_EDGE_FLUX_PLAN.md` | French complex edge-flux |
| 49 | `OPHELIE_COMPLEX_EDGE_FLUX_SOLVER_PLAN.md` | original plan |
| 50 | `build/output/PlateBody_ite_0000000000.vtp` | ParaView field (large) |

---

## E. One-click packaging command

```bash
cd /home/yyc/SPHinXsysSYCL

zip -r /tmp/team7_l2_gpt_upload_round2.zip \
 docs/ophelie/archive/discussion/OPHELIE_GPT_DISCUSSION_TEAM7_L2_EDGE_FLUX.md \
 docs/ophelie/archive/discussion/TEAM7_L2_UPLOAD_MANIFEST.md \
 discussion_bundle/TEAM7_L2_GPT_FIRST_PROMPT_ROUND2.txt \
 discussion_bundle/README.md \
 discussion_bundle/team7_l2_outputs/ \
 OPHELIE_TEAM7_L2_VALIDATION_DEBUG_AND_CURSOR_PLAN.md \
 tests/extra_source_and_tests/3d_examples/test_3d_ophelie_team7_complex_edge_flux/TEAM7_VALIDATION_LOG.md \
 tests/extra_source_and_tests/3d_examples/test_3d_ophelie_team7_complex_edge_flux/README.md \
 tests/extra_source_and_tests/3d_examples/test_3d_ophelie_team7_complex_edge_flux/test_3d_ophelie_team7_complex_edge_flux.cpp \
 tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/team7/electromagnetic_ophelie_team7_validation.h \
 tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/team7/electromagnetic_ophelie_team7_probe.h \
 tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/team7/electromagnetic_ophelie_team7_coil_path_source.h \
 tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/diagnostics/electromagnetic_ophelie_aind_diagnostic.h \
 tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/diagnostics/aphi_team7_native_geometry_config.h \
 tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_edge_flux.h \
 tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_parameters.h \
 tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_cli.h \
 tests/extra_source_and_tests/3d_examples/reference_data/team7/tools/import_team7_jey_comsol_table2.py \
 tests/extra_source_and_tests/3d_examples/particle_generation_TEAM7/particle_generation_team7.cpp

ls -lh /tmp/team7_l2_gpt_upload_round2.zip
```

---

## F. First prompt for ChatGPT

See [`TEAM7_L2_GPT_FIRST_PROMPT_ROUND2.txt`](TEAM7_L2_GPT_FIRST_PROMPT_ROUND2.txt)

---

## G. Reproduce run data (filament L2 baseline)

```bash
cd /home/yyc/SPHinXsysSYCL/build

# 50 Hz baseline + J/B split
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_team7_complex_edge_flux/bin/test_3d_ophelie_team7_complex_edge_flux \
 --reload=1 --team7-level=one-way --native-dp-mm=3 --ophelie-edge-flux-complex=1 \
 --coil-source-model=filament-racetrack --team7-coil-source-scale=1.0 --team7-frequency=50 \
 --team7-reference-dir=$HOME/SPHinXsysSYCL/tests/extra_source_and_tests/3d_examples/reference_data/team7

# P3c ω sweep
 --team7-omega-sweep=1

# P5 feedback
 --ophelie-aind-one-way-feedback=1

# P2b Picard
 --team7-level=picard --self-induction-relax=0.05 --self-induction-max-iter=6
```

---

*Manifest version: Round 2, 2026-06-11*
