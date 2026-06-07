# Stage 10.17 Pre-Commit Cleanup and CLI Fix Plan

> Project: SPHinXsys SYCL electromagnetic A-phi solver  
> Scope: TEAM7 native geometry, matrix-free Laplace A-phi CK solver, pre-commit cleanup  
> Status: before milestone commit  
> Main goal: fix known configuration/preset inconsistency, clean obsolete files safely, and make the repository state suitable for a GitHub milestone commit.

---

## 0. Background and Current Decision

We will **not** change the main numerical route in this cleanup phase.

The confirmed mainline is:

```text
Matrix-free Laplace-form A-phi solver
+ SPHinXsys SYCL CK style
+ Contact multi-body GMRES
+ TEAM7 native geometry / SI unit setup
+ impressed-current path source
```

This cleanup phase is not about improving the solver accuracy, big-air multi-resolution robustness, or Joule-thermal coupling. Those tasks should be continued **after** this pre-commit cleanup is finished and committed.

The current priority is:

1. Fix the known CLI preset parsing bug.
2. Make configuration, reload metadata, scripts, reports, and documentation mutually consistent.
3. Clean or archive obsolete electromagnetic development paths.
4. Prepare a clean milestone commit.

---

## 1. Scope of This Cleanup Phase

### 1.1 Tasks included

This phase includes:

```text
A. Configuration and record consistency audit
B. CLI preset parsing bug fix
C. Reload metadata audit
D. Documentation update
E. File cleanup plan
F. Minimal pre-commit gate tests
G. Commit preparation
```

### 1.2 Tasks explicitly excluded

Do **not** work on the following in this phase:

```text
- Big-air multi-resolution optimization
- Laplace operator local-h / pair-h regularization experiment
- Adaptive Bz interpolation probe implementation
- GMRES/preconditioner tuning
- New Joule-thermal coupling development
- Contact A-penalty research
- divA/curlB advanced diagnostics
- Full TEAM7 <10% profile error target
```

Those are next-stage tasks after the cleanup commit.

---

## 2. Critical Bug: `--team7-case=` CLI Preset Parsing

### 2.1 Problem

In:

```text
tests/extra_source_and_tests/extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_team7_native_geometry_config.h
```

there is likely an off-by-one parsing bug:

```cpp
else if (arg.rfind("--team7-case=", 0) == 0)
{
    if (!applyTeam7NativeCasePreset(arg.substr(12), config))
    {
        std::cerr << "warning: unknown --team7-case preset: " << arg.substr(13) << "\n";
    }
}
```

The prefix string is:

```text
--team7-case=
```

Its length is **13**, not 12.

Current behavior likely passes a string such as:

```text
=uniform-legacy-dp6
=multires-small-dp6-l2
```

instead of:

```text
uniform-legacy-dp6
multires-small-dp6-l2
```

This can silently prevent presets from being applied correctly.

### 2.2 Required fix

Replace all magic-number `substr(N)` parsing with prefix-size based parsing.

Recommended style:

```cpp
const std::string prefix = "--team7-case=";
if (arg.rfind(prefix, 0) == 0)
{
    const std::string preset = arg.substr(prefix.size());
    if (!applyTeam7NativeCasePreset(preset, config))
    {
        std::cerr << "warning: unknown --team7-case preset: " << preset << "\n";
    }
}
```

Apply the same principle to every CLI option using string prefixes:

```text
--team7-case=
--team7-reload-case=
--team7-air-preset=
--team7-dp-mm=
--team7-local-refinement-levels=
--team7-discretization=
--coil-source-model=
--racetrack-inset-mm=
--racetrack-z-inset-mm=
--target-power=
```

Do not use hard-coded substring indices when the prefix string can be used directly.

### 2.3 Required grep checks

Run:

```bash
rg "substr\([0-9]+\)" tests/extra_source_and_tests/extra_src/shared/electromagnetic_dynamics
rg "rfind\(.*0\)" tests/extra_source_and_tests/extra_src/shared/electromagnetic_dynamics
rg "team7-case|team7-reload-case|team7-air-preset|local-refinement|discretization" tests/extra_source_and_tests
```

For every CLI parser match, verify that:

```text
parsed_value = arg.substr(prefix.size())
```

rather than a hard-coded numeric index.

### 2.4 Acceptance criteria

The fix is accepted only if the program clearly prints or records the resolved preset information, for example:

```text
TEAM7 native config:
  team7_case                 = uniform-legacy-dp6
  reload_case_id             = team7_native_uniform_legacy_dp6
  air_preset                 = legacy
  discretization             = uniform
  dp_reference               = 0.006 m
  local_refinement_levels    = 0
```

and for MR:

```text
TEAM7 native config:
  team7_case                 = multires-small-dp6-l1
  reload_case_id             = team7_native_multires_small_dp6_l1
  air_preset                 = small
  discretization             = multires
  dp_reference               = 0.006 m
  local_refinement_levels    = 1
```

If a preset name is unknown, the executable should emit a warning and preferably fail in audit tests.

---

## 3. Configuration and Record Consistency

### 3.1 Why this is required

The current records mention several case labels, including:

```text
small uniform
small MR L1
legacy uniform
legacy MR L2
multires-small-dp6-l2
uniform-legacy-dp6
```

Some of these labels are ambiguous. For the next commit, every case should be described by the same fields:

```text
air_preset
dp_reference
discretization
local_refinement_levels
reload_case_id
source_model
target_power
```

Avoid relying on informal names such as "legacy MR" or "small l2" without the full configuration.

### 3.2 Standard naming convention

Use the following naming convention:

```text
team7_native_<discretization>_<air_preset>_dp<dp_mm>_l<level>
```

Examples:

```text
team7_native_uniform_small_dp6_l0
team7_native_uniform_legacy_dp6_l0
team7_native_multires_small_dp6_l1
team7_native_multires_legacy_dp6_l1
team7_native_multires_legacy_dp6_l2
team7_native_multires_big_dp6_l2
```

For uniform cases, `l0` is acceptable even if no adaptive refinement is used.

### 3.3 Required fields in `team7_native_geometry.txt`

The reload metadata file should contain at least:

```text
case_id
air_preset
discretization
dp_reference
local_refinement_levels
air_lower_bound
air_upper_bound
coil_lower_bound
coil_upper_bound
plate_lower_bound
plate_upper_bound
number_of_air_particles
number_of_coil_particles
number_of_plate_particles
smoothing_length_ratio_min_air
smoothing_length_ratio_max_air
smoothing_length_ratio_mean_air
smoothing_length_ratio_min_coil
smoothing_length_ratio_max_coil
smoothing_length_ratio_mean_coil
smoothing_length_ratio_min_plate
smoothing_length_ratio_max_plate
smoothing_length_ratio_mean_plate
```

If the current file format is kept as a simple text file, use stable key-value pairs such as:

```text
case_id = team7_native_multires_small_dp6_l1
air_preset = small
discretization = multires
dp_reference = 0.006
local_refinement_levels = 1
```

### 3.4 Required runtime audit

When a TEAM7 native test loads reload particles, it should print:

```text
Requested reload case:
Loaded reload metadata:
Resolved TEAM7 config:
Particle counts:
SmoothingLengthRatio statistics:
```

If the requested reload case and the metadata case do not match, the test should fail.

Recommended behavior:

```cpp
if (requested_case_id != metadata_case_id)
{
    std::cerr << "ERROR: requested TEAM7 reload case does not match loaded metadata.\n";
    std::cerr << "  requested = " << requested_case_id << "\n";
    std::cerr << "  loaded    = " << metadata_case_id << "\n";
    return 1;
}
```

### 3.5 Required script consistency check

Audit these scripts:

```text
tests/extra_source_and_tests/3d_examples/run_team7_native_si_smoke.sh
tests/extra_source_and_tests/3d_examples/run_team7_uniform_legacy_bz_vtp.sh
tests/extra_source_and_tests/3d_examples/run_team7_highres_bz.sh
tests/extra_source_and_tests/3d_examples/sync_team7_native_si_reload.sh
```

For each script, record:

```text
which executable it calls
which reload case it uses
which dp it assumes
which air preset it assumes
whether it uses uniform or MR
whether it expects Bz reference output
```

Each script should print the resolved config before solving.

---

## 4. Documentation Update Required Before Commit

### 4.1 Main milestone document

Create or update:

```text
tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_STAGE10_17_19_PRECOMMIT_CLEANUP_RECORD.md
```

This document should contain:

```text
1. What was fixed
2. Which CLI parsing bugs were found
3. Which presets were verified
4. Which reload metadata fields were added
5. Which files were cleaned or archived
6. Which gate tests were run
7. Which known problems remain
```

### 4.2 Update existing records

Update:

```text
tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_TEAM7_SMALL_AIR_MULTIRES_PLAN.md
tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_STAGE10_17_TEAM7_NATIVE_RECORD.md
tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/STAGE10_APHI_TEST_CASE_INDEX.md
tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/STAGE10_APHI_FILE_MAP_FOR_DISCUSSION.md
tests/extra_source_and_tests/extra_src/shared/electromagnetic_dynamics/README.md
```

### 4.3 Required wording corrections

Do not claim:

```text
MR has been validated.
Big-air MR is ready.
TEAM7 quantitative validation is fully closed.
```

Use:

```text
MR infrastructure has been implemented.
Small-air MR currently requires consistency and stability audit.
Legacy uniform air improves Bz profile accuracy.
Big-air MR should be attempted only after the MR reload/operator/probe audit.
TEAM7 native SI validation has reached engineering-pass level but not final quantitative benchmark level.
```

### 4.4 Suggested English summary for GitHub

Use this wording in the milestone document:

```text
Stage 10.17 established the SPHinXsys SYCL CK matrix-free Laplace A-phi mainline for the TEAM7 native geometry. The SI-unit geometry, path-based impressed coil source, Contact GMRES solve, and Bz reference probing are now connected in a single production-oriented workflow. Uniform legacy-air runs indicate that enlarging the air domain improves the Bz profile error, confirming domain truncation as a relevant error source. The multi-resolution air infrastructure is implemented but not yet validated; current MR cases require reload metadata, CLI preset, operator, and probe consistency audits before big-air MR can be used as the production route.
```

---

## 5. File Cleanup Plan

### 5.1 Cleanup principle

Do not delete files merely because they are old. Delete or archive files only if they satisfy all three conditions:

```text
1. They are not included by any active test or production file.
2. They belong to an abandoned route.
3. Their historical information is already preserved in Git history or a short archive README.
```

### 5.2 Must keep

Keep the current production CK A-phi route:

```text
tests/extra_source_and_tests/extra_src/shared/electromagnetic_dynamics/
```

Especially keep:

```text
aphi_matrix_free_operator_ck.*
aphi_laplace_ck.*
aphi_grad_phi_coupling_ck.*
aphi_div_sigma_a_coupling_ck.*
aphi_phi_gauge_penalty_ck.*
aphi_reaction_ck.*
aphi_joule_heating_ck.*
aphi_gmres_solver_ck.*
aphi_multibody_contact_gmres_ck.*
aphi_matrix_free_solve_ck.*
aphi_block_vector_ops_ck.*
aphi_block_jacobi_preconditioner_ck.*
aphi_variables_ck.*
aphi_field_names_ck.h
diagnostics/aphi_team7_native_*.h
test_helpers/*
benchmark/*
```

Keep TEAM7 native production cases:

```text
tests/extra_source_and_tests/3d_examples/particle_generation_em/
tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_team7_native_geometry_audit/
tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_team7_native_source_rhs_audit/
tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_team7_native_vacuum_source_b_sanity/
tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_team7_native_geometry_bz_reference_probe/
tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_team7_native_geometry_source_driven_smoke/
tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_team7_native_joule_post_smoke/
```

Keep core gate tests:

```text
test_3d_aphi_ck_fused_apply_vs_debug_lhs
test_3d_aphi_ck_contact_apply_vs_monolithic
test_3d_aphi_ck_contact_fused_apply_equivalence
test_3d_aphi_ck_gmres_contact_interface_flux_matched_manufactured
test_3d_aphi_ck_joule_uniform_field_analytic_verification
```

### 5.3 Safe archive/delete candidates

The following paths are candidates for deletion or archival, but only after `rg` confirms no active dependency:

```text
tests/extra_source_and_tests/3d_examples/legacy_aphi_archive/
tests/extra_source_and_tests/extra_src/shared/legacy_aphi_archive/
tests/extra_source_and_tests/extra_src/shared/aphi_sphinxsys/
```

Reasons:

```text
- legacy_aphi_archive contains abandoned curl-curl / Hessian / Eigen sparse / old matrix-free prototypes.
- aphi_sphinxsys is not the current CK-style production route.
- Current TEAM7 native results are produced by electromagnetic_dynamics, not by these paths.
```

Before deleting, run:

```bash
rg "legacy_aphi_archive" tests/extra_source_and_tests
rg "aphi_sphinxsys" tests/extra_source_and_tests
rg "#include .*aphi_sphinxsys|#include .*legacy_aphi" tests/extra_source_and_tests
```

If no active dependency exists, either remove them or move them outside the active source tree.

### 5.4 Research files: move, do not delete yet

Do not delete immediately:

```text
A-divergence penalty studies
Contact A-penalty studies
GMRES/preconditioner sensitivity diagnostics
BiCGStab historical comparisons
curlB/divA diagnostic tests
alternate_krylov/
```

Recommended destination:

```text
tests/extra_source_and_tests/extra_src/shared/electromagnetic_dynamics/experimental/
tests/extra_source_and_tests/extra_src/shared/electromagnetic_dynamics/diagnostics/research/
tests/extra_source_and_tests/3d_examples/aphi_ck_diagnostics/
```

Rationale:

```text
These files are not production, but they document why some routes were not adopted.
They may still be useful for paper writing, numerical-method discussion, or future debugging.
```

### 5.5 Empty file and generated file scan

Run:

```bash
find tests/extra_source_and_tests -empty -print
find tests/extra_source_and_tests -name "*.vtp" -o -name "*.vtu" -o -name "*.pvd"
find tests/extra_source_and_tests -name "build" -type d
find tests/extra_source_and_tests -name "output" -type d
find tests/extra_source_and_tests -name "restart" -type d
find tests/extra_source_and_tests -name "*.csv" -size +5M
find tests/extra_source_and_tests -name "*.xml" -size +5M
```

Do not commit generated visualization output unless it is deliberately used as small reference data.

### 5.6 CMake cleanup

After deleting or moving files, run:

```bash
rg "legacy_aphi_archive|aphi_sphinxsys|test_3d_em_|curl_curl|hessian|frequency" tests/extra_source_and_tests -g "CMakeLists.txt"
```

Remove stale CMake entries.

Then configure from a clean build directory.

---

## 6. Suggested Repository Organization After Cleanup

Do not perform a large structural refactor in this commit unless necessary. But the target organization should be documented.

Recommended long-term structure:

```text
tests/extra_source_and_tests/extra_src/shared/electromagnetic_dynamics/
├── core/
├── operators/
├── solvers/
├── team7_native/
├── diagnostics/
├── test_helpers/
├── benchmark/
└── experimental/
```

For this pre-commit cleanup, a lighter version is acceptable:

```text
electromagnetic_dynamics/
├── production files at current level
├── diagnostics/
│   ├── team7 native helper files
│   └── research diagnostics
├── test_helpers/
├── benchmark/
├── alternate_krylov/
└── README.md
```

Do not move too many files at once if that causes large include-path churn.

---

## 7. Minimal Pre-Commit Gate Tests

After fixing CLI parsing and cleanup, run at least:

```bash
# TEAM7 native geometry and source
./test_3d_aphi_ck_team7_native_geometry_audit
./test_3d_aphi_ck_team7_native_source_rhs_audit
./test_3d_aphi_ck_team7_native_vacuum_source_b_sanity
./test_3d_aphi_ck_team7_native_geometry_bz_reference_probe

# Core CK operator gates
./test_3d_aphi_ck_fused_apply_vs_debug_lhs
./test_3d_aphi_ck_contact_apply_vs_monolithic
./test_3d_aphi_ck_contact_fused_apply_equivalence

# Contact GMRES and Joule postprocess
./test_3d_aphi_ck_gmres_contact_interface_flux_matched_manufactured
./test_3d_aphi_ck_joule_uniform_field_analytic_verification
```

If executable names differ in the local build directory, run their corresponding CTest targets.

Recommended CTest filter:

```bash
ctest -R "aphi_ck|team7_native" --output-on-failure
```

If the full `aphi_ck` suite is too large, run the minimal gate list above and record exactly which tests were run.

---

## 8. Commit Plan

### 8.1 Before commit

Run:

```bash
git status --short
git diff --stat
git diff -- tests/extra_source_and_tests/extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_team7_native_geometry_config.h
```

Check that:

```text
- CLI parsing bug is fixed.
- Metadata audit is added or updated.
- Documentation is updated.
- Legacy files are deleted/moved only if no active include remains.
- No build/output/VTP/restart files are accidentally staged.
```

### 8.2 Recommended staging

Stage source and docs first:

```bash
git add tests/extra_source_and_tests/extra_src/shared/electromagnetic_dynamics
git add tests/extra_source_and_tests/3d_examples/particle_generation_em
git add tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_team7_native_geometry_audit
git add tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_team7_native_source_rhs_audit
git add tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_team7_native_vacuum_source_b_sanity
git add tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_team7_native_geometry_bz_reference_probe
git add tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_team7_native_geometry_source_driven_smoke
git add tests/extra_source_and_tests/3d_examples/test_3d_aphi_ck_team7_native_joule_post_smoke
git add tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/STAGE10_APHI_TEST_CASE_INDEX.md
git add tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/STAGE10_APHI_FILE_MAP_FOR_DISCUSSION.md
git add tests/extra_source_and_tests/3d_examples/docs/electromagnetic_aphi/CURSOR_APHI_STAGE10_17_19_PRECOMMIT_CLEANUP_RECORD.md
```

Stage deletions only after review:

```bash
git add -u tests/extra_source_and_tests/3d_examples/legacy_aphi_archive
git add -u tests/extra_source_and_tests/extra_src/shared/legacy_aphi_archive
git add -u tests/extra_source_and_tests/extra_src/shared/aphi_sphinxsys
```

### 8.3 Recommended commit split

Prefer two commits:

```text
Commit 1:
Stage 10.17: fix TEAM7 native preset parsing and metadata audit

Commit 2:
Archive obsolete electromagnetic A-phi prototype paths
```

If the current worktree is too tangled, one combined milestone commit is acceptable:

```text
Stage 10.17: TEAM7 native A-phi cleanup and legacy archive
```

### 8.4 Suggested commit message

```text
Stage 10.17: fix TEAM7 preset parsing and cleanup EM A-phi files

- Fix TEAM7 native CLI preset parsing by replacing hard-coded substring indices with prefix-size based parsing.
- Add or update TEAM7 native configuration/reload metadata audit.
- Standardize TEAM7 native case naming for uniform and multi-resolution air domains.
- Update Stage 10.17 documentation and test-case index.
- Keep the SPHinXsys SYCL CK matrix-free Laplace A-phi solver as the production route.
- Archive obsolete non-CK, curl-curl, Hessian, and sparse-matrix prototype paths.
- Preserve research diagnostics for A-divergence, Contact A-penalty, and Krylov sensitivity studies.
```

---

## 9. Known Issues That Should Remain Open After This Commit

Do not try to solve these in this cleanup phase. Record them as known issues:

```text
1. TEAM7 Bz profile is engineering-pass but not final quantitative validation.
   Current legacy uniform real L2 is about 13%, imaginary L2 about 21%, depending on case and probe setup.

2. Peak Bz x-location still shows an offset relative to the reference.

3. Multi-resolution air is implemented but not validated.
   Small MR and legacy MR require reload/operator/probe audits.

4. Current nearest-particle Bz probe may not be suitable for MR comparison.
   Adaptive kernel interpolation probe should be added next.

5. The pairwise Laplace regularization currently uses global reference smoothing length.
   A local-h or pair-h diagnostic mode should be tested next for MR.

6. Big-air MR should not be treated as production until small/legacy MR L1 is stable.

7. Joule heat coupling is available as a pipeline but should not be the focus of this cleanup commit.
```

---

## 10. Next Stage After This Commit

After this cleanup is committed, continue with:

```text
Stage 10.18 / 10.19:
MR reload/operator/probe stabilization
```

Recommended order:

```text
1. Add MR reload audit test.
2. Add adaptive Laplace constant/linear MMS tests.
3. Add finite first-matvec and positive-Jacobi diagnostics.
4. Add adaptive kernel-interpolated Bz probe.
5. Run small uniform dp6 vs small MR L1.
6. Run legacy uniform dp6 vs legacy MR L1.
7. Only then try big-air MR.
```

The optimization of the solver and multi-resolution big-air strategy should begin only after the repository has a clean and reproducible Stage 10.17 commit.

---

## 11. Short Instruction for Cursor

```text
Please do not optimize the solver in this step. First fix the TEAM7 native CLI preset parser, make the reload metadata and records consistent, update the documentation, and clean/archive obsolete electromagnetic prototype files. Keep the matrix-free Laplace A-phi CK route as the production route. Do not claim MR validation yet. After the cleanup passes the minimal TEAM7/native A-phi gate tests, prepare a milestone commit.
```
