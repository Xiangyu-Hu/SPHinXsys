# OPHELIE TEAM7 — complete material pack for ChatGPT discussion

> repository root：`/home/yyc/SPHinXsysSYCL`
> builddirectory：`/home/yyc/SPHinXsysSYCL/build`
> externalReference：[COMSOL — Solving TEAM Problem 7](https://www.comsol.com/blogs/solving-team-problem-7-acdc-module/)

---

## 1. what we are attempting（one sentence）

implement in **SPHinXsysSYCL** **OPHELIE-style particle method** induction heating：apply known coil current density **J**，use ****discrete Biot–Savart** for aluminum plate **A/B**，then use **σ constitutive relation + optional scalar φ fix** for **E, induction J, Joule heat**；GeometryonIn progressfrom **1 m analytic scaffold** pivot to **TEAM7 native STL**（`particle_generation_em`）， and and `TEAM7-reference` **COMSOL/muFEM curves** for comparison——**curve comparison code not yet implemented**。

---

## 2. vs reference comparison: current actual setup (must read)

### 2.1 comparison target（Literature / repositoryReference）

| Reference | path | content |
|------|------|------|
| probe definitions | [`docs/TEAM7-reference/reference_data/team7/TEAM7_probe_definitions.csv`](../../docs/TEAM7-reference/reference_data/team7/TEAM7_probe_definitions.csv) | A1–B1: Bz line x∈[0,288] mm, y=72, z=34 mm, 50/200 Hz; Jy line z≈18.99 mm |
| Bz Reference (mT) | [`TEAM7_Bz_A1_B1_reference_mT.csv`](../../docs/TEAM7-reference/reference_data/team7/TEAM7_Bz_A1_B1_reference_mT.csv) | 17 x points, 50/200 Hz phase components |
| Bz Reference (Gauss) | [`TEAM7_Bz_A1_B1_reference_Gauss.csv`](../../docs/TEAM7-reference/reference_data/team7/TEAM7_Bz_A1_B1_reference_Gauss.csv) | same, units Gauss |
| A2–B2 line | [`TEAM7_Bz_A2_B2_reference_mT.csv`](../../docs/TEAM7-reference/reference_data/team7/TEAM7_Bz_A2_B2_reference_mT.csv) | y=144 mm |
| COMSOL Note | [COMSOL Blog TEAM7](https://www.comsol.com/blogs/solving-team-problem-7-acdc-module/) | 2742 turns, 1 A/turn, aluminum plate eddy currents，50/200 Hz |

**Important：current OPHELIE code does not read the CSVs above or output along-probe interpolation error tables。** comparison is **planned P0**, not a completed result.

### 2.2 geometry and units（native STL mode，`--native-stl`）

| item | Setup |
|----|------|
| STL source | [`particle_generation_em/data/coil.stl`](../../tests/extra_source_and_tests/3d_examples/particle_generation_em/data/coil.stl)、[`plate.stl`](../../tests/extra_source_and_tests/3d_examples/particle_generation_em/data/plate.stl)（mm coordinates） |
| import scaling | `stl_scale_to_meter_ = 1e-3`（mm → m） |
| particle spacing | default `dp_mm = 6` → `dp = 0.006 m`（`--native-dp-mm=` can change） |
| computational domain (small air, default) | air box mm：[-50,350]×[-50,350]×[-50,200]，scaled to m |
| EM parameter and body | only **PlateBody**（conductive plate）、**CoilSourceBody**（Coil）；**no air particles**parameterEM |
| Relax（optional） | default **coil + plate + air** SYCL-CK body-fitted；Reload only write **Coil + Plate** two bodies |

implementation see：[`electromagnetic_ophelie_team7_native_geometry.h`](../../tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_team7_native_geometry.h)

### 2.3 where particles come from（initial geometry artifact, not EM initial field）

three mutually exclusive modes（Companion to [`mr_free_stream_around_cylinder`](../../tests/tests_sycl/2d_examples/test_2d_free_stream_around_cylinder_mr_sycl/mr_free_stream_around_cylinder.cpp) same idea）：

| mode | command row | Particle |
|------|--------|------|
| lattice | default / `--skip-relax` | `generateParticles<Lattice>`，no relax yet |
| Relax | `--relax=1` | lattice → SYCL relax (1000 steps default)→ `Reload.xml` |
| Reload | `--reload=1` | read from `./reload/Reload.xml` **Position + VolumetricMeasure** |

your most recent successful run：

```bash
cd /home/yyc/SPHinXsysSYCL/build
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_team7/bin/test_3d_ophelie_team7 \
 --native-stl --reload=1 --no-phi --sigma=1e4 --ophelie-smoke --state_recording=0
```

Log：`n_plate=7700`, `n_coil=7712`, `source=reload`, `passed=1`。

### 2.4 EM initial/boundary conditions — how excitation is applied（not COMSOL boundary conditions）

**not solving coil voltage or full Maxwell; explicit source + plate conductive constitutive relation.**

#### Step A — plate conductivity

- operator ：`AssignOphelieGlassSigmaCK`
- each plate particle：`sigma[i] = sigma_glass`（CLI `--sigma=`，you used `1e4`；native defaultphysical Al `3.54e7` unless user overrides）

#### Step B — coil excitation current density J

- operator ：`InitializeOphelieCoilSourceCK`
- each coil particle i：
 - `r_xy = position_i - coil_center`（native：`coil_center` = coil particle centroid）
 - `e_theta = (-dy/r, dx/r, 0)`
 - **`JSrcReal = J0 * e_theta`**，`JSrcImag = 0`
- **J0**（native TEAM7）：
 - `N_turns = 2742`, `I_per_turn = 1 A`
 - `J0 = N * I / A_cross`（`A_cross` ≈ coil particle bounding-box radial × axial scale，**not strict turn cross-section**）
 - your log：`J0 ≈ 146975 A/m²`

code：[`electromagnetic_ophelie_source.hpp`](../../tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_source.hpp)

harmonic approximation：real part stored in `JSrcReal`；frequency `f=50 Hz` → `ω=2πf`；subsequently **E imaginary part** passed `-ω·A` represents 90° phase。

#### Step C — Coil → plate：Biot–Savart（discrete summation）

- operator ：`ComputeOphelieCoilToGlassBiotSavartCK`
- for each plate particle i, over coil particles j：
 - `r = x_i - x_j`
 - current element `jv = JSrcReal_j * V_j`（V = particle volume）
 - `A += (μ₀/4π) * jv / |r|`
 - `B += (μ₀/4π) * jv × r / |r|³`
- write to：`ACoilReal`, `BCoilReal`（imaginary part as 0）
- then `CombineOphelieCoilAndInducedVector PotentialCK` → **`ASrcReal`, `BSrcReal`**

code：[`electromagnetic_ophelie_biot_savart.hpp`](../../tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_biot_savart.hpp)

#### Step D — Level 0 (stops here with `--no-phi`)

- operator ：`ComputeOphelieEJQFromASrcNoPhiCK`
- each plate particle：
 - **`EImag = -ω * ASrcReal`**
 - **`JImag = σ * EImag`**
 - **`JouleHeat = ½ (JReal·EReal + JImag·EImag)`**（harmonic power）

code：[`electromagnetic_ophelie_postprocess.hpp`](../../tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_postprocess.hpp)

#### Step E — φ fix（on by default; you used `--no-phi` to disable）

1. `setupOpheliePhiImagRhsProblem`： by `ASrc`、`σ` build `PhiRhsImag`
2. `solvePhiImag`（PCG/GMRES/Jacobi）：SPH discretized Laplace + gauge → `PhiImag`
3. `ComputeOphelieScalarPhiGradientCK` → `grad_phi_imag`
4. `ComputeOphelieEJQWithPhiCK`：
 - **`EImag = -grad_phi_imag - ω * ASrcReal`**
 - **`JImag = σ * EImag`**

code：[`electromagnetic_ophelie_phi.hpp`](../../tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_phi.hpp)、[`electromagnetic_ophelie_phi_gmres.h`](../../tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_phi_gmres.h)

native 7700 plate particles: φ often does not converge（`rel_res_linf` can reach O(10²)）。

#### Step F — PowerScaling（option B, **problematic when comparing to COMSOL 1 A/turn**）

- operator ：`ScaleOphelieElectromagneticFieldsCK`
- `P_raw = Σ (JouleHeat * V)`（volume-weighted sum）
- `power_scale = P_target / P_raw`（`P_target = 50 kW` default）
- **`field_scale = sqrt(power_scale)`**
- Scaling：**A, B, E, J, Phi** × `field_scale`；**JouleHeat** × `power_scale`

your log：`P_raw≈8e-4 W` → `field_scale≈7895` → print's `max_BSrc≈9.74` etc. is ****post-scaling** quantities。

code：[`electromagnetic_ophelie_diagnostics.h`](../../tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_diagnostics.h) in `computeOpheliePowerScalingFactors`

**still no `--no-power-scaling` switch.**

### 2.5 what results are actually collected

#### terminal one-line summary（`test_3d_ophelie_team7` end）

| Quantity | meaning |
|----|------|
| `n_plate`, `n_coil` | Particle count |
| `max_ASrc`, `max_BSrc` | plate source field peak（post-scaling） |
| `max_EImag`, `max_JImag` | plate E/J peak values（post-scaling） |
| `max_ACoil` | coil contribution to A on plate, peak |
| `max_PhiImag`, `phi_rel_res` | φ solution (if φ enabled) |
| `P_raw`, `P_scaled`, `field_scale`, `I_eff` | power scaling diagnostics |
| `divJ_L0`, `divJ_phi`, `divJ_red` | SPH discrete div(J) relative quantity |
| `passed` | heuristic pass（`--ophelie-smoke` relaxed φ） |

#### VTP (with `--state_recording=1`, under `build/output/`)

written variables see [`test_3d_ophelie_team7.cpp`](../../tests/extra_source_and_tests/3d_examples/test_3d_ophelie_team7/test_3d_ophelie_team7.cpp) end `BodyStatesRecordingToVtp`：

- Coil：`JSrcReal`
- plate：`ASrcReal`, `BSrcReal`, `EImag`, `JImag`, `JouleHeat`, `sigma`；if φ enabled also has `PhiImag`, `grad_phi_imag`

#### Relax VTP（`--relax=1`）

- `PlateBody_ite_*`, `CoilSourceBody_ite_*` etc. （`output/`）

#### not yet collected

- along TEAM7 probe lines **Bz(x)、Jy(x)** and CSV 'serror
- 200 Hz second frequency（parameter can set `frequency`， but no comparison script yet）
- Coil/plateinterior 's COMSOL consistentphaseapproximate（COMSOL use `sign(real)*abs`）

---

## 3. main executables and test cases（file links）

### 3.1 main integration test (TEAM7)

| File | Note |
|------|------|
| [`test_3d_ophelie_team7.cpp`](../../tests/extra_source_and_tests/3d_examples/test_3d_ophelie_team7/test_3d_ophelie_team7.cpp) | **main executable**：analytic / native STL、relax、reload、EM、passed |
| [`test_3d_ophelie_team7/README.md`](../../tests/extra_source_and_tests/3d_examples/test_3d_ophelie_team7/README.md) | RunNote、command example |
| [`test_3d_ophelie_team7/reload/README.md`](../../tests/extra_source_and_tests/3d_examples/test_3d_ophelie_team7/reload/README.md) | Reload.xml approximate |
| [`test_3d_ophelie_team7/CMakeLists.txt`](../../tests/extra_source_and_tests/3d_examples/test_3d_ophelie_team7/CMakeLists.txt) | build; copy STL to `bin/input/` |

executable (after build)：

```text
/home/yyc/SPHinXsysSYCL/build/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_team7/bin/test_3d_ophelie_team7
```

### 3.2 particle generation + relax prototype (you already tested)

| File | Note |
|------|------|
| [`particle_generation_em.cpp`](../../tests/extra_source_and_tests/3d_examples/particle_generation_em/particle_generation_em.cpp) | TEAM7 STL lattice + air + SYCL relax 1000 steps; body names **Coil/Plate/Air** (different from OPHELIE) |
| [`particle_generation_em/data/`](../../tests/extra_source_and_tests/3d_examples/particle_generation_em/data/) | coil.stl, plate.stl |

### 3.3 box OPHELIE regression

| File | Note |
|------|------|
| [`test_3d_ophelie.cpp`](../../tests/extra_source_and_tests/3d_examples/test_3d_ophelie/test_3d_ophelie.cpp) | box glass + coil；`--ophelie-compare-level0` |
| [`test_3d_ophelie/README.md`](../../tests/extra_source_and_tests/3d_examples/test_3d_ophelie/README.md) | |

### 3.4 operator-level tests (already passed, unrelated to TEAM7 curves)

| Test | path | Validation |
|------|------|------|
| Biot B direction | [`test_3d_ophelie_biot_savart_direction.cpp`](../../tests/extra_source_and_tests/3d_examples/test_3d_ophelie_biot_savart_direction/test_3d_ophelie_biot_savart_direction.cpp) | Bz>0 |
| Coil J direction/Scaling | [`test_3d_ophelie_coil_source_direction.cpp`](../../tests/extra_source_and_tests/3d_examples/test_3d_ophelie_coil_source_direction/test_3d_ophelie_coil_source_direction.cpp) | e_θ，J0×2 |
| Level0 uniform A | [`test_3d_ophelie_level0_uniform_A.cpp`](../../tests/extra_source_and_tests/3d_examples/test_3d_ophelie_level0_uniform_A/test_3d_ophelie_level0_uniform_A.cpp) | E=-ωA, J=σE |
| I Scaling | [`test_3d_ophelie_scaling.cpp`](../../tests/extra_source_and_tests/3d_examples/test_3d_ophelie_scaling/test_3d_ophelie_scaling.cpp) | J0×2 → A,B,E,J×2, P×4 |

### 3.5 A–φ branch TEAM7 reference (read-only comparison, not yet hooked into OPHELIE build)

| File | Note |
|------|------|
| [`docs/TEAM7-reference/team7/aphi_team7_native_stl_geometry_helpers.h`](../../docs/TEAM7-reference/team7/aphi_team7_native_stl_geometry_helpers.h) | another TEAM7 STL + relax pipeline |
| [`test_3d_aphi_ck_team7_native_geometry_audit.cpp`](../../docs/TEAM7-reference/test_3d_aphi_ck_team7_native_geometry_audit/test_3d_aphi_ck_team7_native_geometry_audit.cpp) | geometry audit |

### 3.6 SYCL relax official reference

| File | Note |
|------|------|
| [`particle_relaxation_single_resolution_3D.cpp`](../../tests/tests_sycl/3d_examples/test_3d_particle_relaxation_single_resolution_sycl/particle_relaxation_single_resolution_3D.cpp) | official SYCL single-resolution relax |
| [`mr_free_stream_around_cylinder.cpp`](../../tests/tests_sycl/2d_examples/test_2d_free_stream_around_cylinder_mr_sycl/mr_free_stream_around_cylinder.cpp) | **relax + reload same main** example |

---

## 4. OPHELIE source file modules（`extra_src/shared/electromagnetic_ophelie/`）

| File | responsibility |
|------|------|
| [`electromagnetic_ophelie.h`](../../tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie.h) | umbrella include |
| [`electromagnetic_ophelie_parameters.h`](../../tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_parameters.h) | ω、σ、φ solver、target_P etc. |
| [`electromagnetic_ophelie_cli.h`](../../tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_cli.h) | `--no-phi`, `--native-stl`, `--ophelie-smoke`, `--sigma=` etc. |
| [`electromagnetic_ophelie_source.h/.hpp`](../../tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_source.h) | coil J initialization |
| [`electromagnetic_ophelie_biot_savart.h/.hpp`](../../tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_biot_savart.h) | Biot–Savart |
| [`electromagnetic_ophelie_postprocess.h/.hpp`](../../tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_postprocess.h) | Level0 E/J/Q |
| [`electromagnetic_ophelie_phi.h/.hpp`](../../tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_phi.h) | φ Laplace, gradient, E/J with φ |
| [`electromagnetic_ophelie_phi_gmres.h`](../../tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_phi_gmres.h) | PCG/GMRES（Note L2 vs L∞ residualLog） |
| [`electromagnetic_ophelie_laplace.h`](../../tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_laplace.h) | SPH discretized Laplace |
| [`electromagnetic_ophelie_diagnostics.h`](../../tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_diagnostics.h) | divJ、PowerScaling factor |
| [`electromagnetic_ophelie_team7_geometry.h`](../../tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_team7_geometry.h) | 1 m analytic TEAM7 scaffold |
| [`electromagnetic_ophelie_team7_native_geometry.h`](../../tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_team7_native_geometry.h) | **native STL TEAM7** + native relax |
| [`electromagnetic_ophelie_relaxation.h`](../../tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_relaxation.h) | single-body SYCL relax (analytic use) |
| [`electromagnetic_ophelie_self_induction.h`](../../tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_self_induction.h) | experimental self-induction（**not used for TEAM7 acceptance**） |
| [`electromagnetic_ophelie_field_names.h`](../../tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_field_names.h) | VTP variable names |
| [`electromagnetic_ophelie_observables.h`](../../tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_observables.h) | host reduction, metrics |

build entry：[`tests/extra_source_and_tests/extra_src/CMakeLists.txt`](../../tests/extra_source_and_tests/extra_src/CMakeLists.txt) → `extra_sources_3d` library。

---

## 5. documents and review records

| File | Note |
|------|------|
| [`OPHELIE_CODE_REVIEW_NEXT_STEPS.md`](../../OPHELIE_CODE_REVIEW_NEXT_STEPS.md) | ChatGPT review manifest, tests to add, TEAM7 order recommended |
| [`docs/ophelie/TEAM7_GPT_DISCUSSION_PACKAGE.md`](TEAM7_GPT_DISCUSSION_PACKAGE.md) | **this document** |
| [`FRENCH_OPHELIE_PARTICLE_DEVELOPMENT_PLAN.md`](../../FRENCH_OPHELIE_PARTICLE_DEVELOPMENT_PLAN.md) | French plan (if still present) |
| [`french_ophelie_particle_plan_and_starter/`](../../french_ophelie_particle_plan_and_starter/) | early starter code (not hooked into main build) |
| [`biot_reference_selection_for_cursor//`](../../biot_reference_selection_for_cursor//) | A–φ CK Reference（read-only, not hooked into OPHELIE） |

---

## 6. recorded run results (you can paste to GPT)

### 6.1 native + reload + no-phi（2025-06，success）

```
native-stl n_plate=7700 n_coil=7712 dp=0.006 source=reload
J0=146975 sigma=10000
phi_correction=0
max_BSrc=9.74111 max_EImag=116.419 max_JImag=1.16419e+06
P_raw=0.000802152 P_scaled=50000 field_scale=7895.08
divJ_L0=0.181924
passed=1
```

### 6.2 native + skip-relax + φ（φ poor）

- `phi_rel_res_linf` ~ 17–22，PCG hits 6000 steps
- Log `rel_res_l2_vol` ~ 5e-4 and `rel_res_linf` ~ 20 **not the same metric**

### 6.3 analytic 1 m scaffold

- `--skip-relax --no-phi`：``passed=1` easier
- inconsistent with COMSOL geometry

---

## 7. Model gap vs COMSOL TEAM7（discussion points）

| COMSOL / TEAM7 Benchmark | current OPHELIE attempt |
|---------------------|------------------|
| full-field EM / eddy-current solve | Coil **known J**，plate **σE** + optional φ |
| 2742-turn fine geometry and excitation | STL Particle + **e_θ** + bounding-box **J0** |
| probe Bz, Jy **curve** validation | ****not yet implemented** along-line interpolation |
| 1 A/turn under absolute field | **50 kW PowerScaling** distorts A/B/E/J amplitude |
| air-domain field | no air EM particles |
| φ / gauge | SPH Laplace + gauge penalty，poor convergence |

---

## 8. recommended open issues for ChatGPT discussion (manifest)

1. **acceptance first line**：use `TEAM7_Bz_A1_B1_reference_mT.csv` @ 50 Hz as criterion？what tolerance？
2. **PowerScaling**：when comparing to COMSOL, must disable `target_P=50kW`（need to add CLI）？
3. **J0 and N×I**：bounding-box area vs native turn cross-section，is acceptance error source？
4. **Biot particle summation**：7700×7712 enough？need multipole or grid convergence study？
5. **should φ enter TEAM7 v1**：still first **no-phi + Bz line**？
6. **`e_θ` axis**：handling when STL coil geometry is not fully symmetric？
7. **phase**：COMSOL `sign(real)*abs` vs we `EImag=-ωA` how to align with reference CSV phase0/90？
8. **divJ≈0.18** use as plate conductor weak-solution metric？
9. **A–φ branch** `TEAM7-reference` relationship: align same STL/dp before comparing？
10. **Next stepsImplementationPriority**：probe CSV comparison vs disable scaling vs φ/GMRES？

---

## 9. Recommended attachment order for ChatGPT

1.this document（or sections 2, 7, 8）
2. `test_3d_ophelie_team7.cpp`（main flow）
3. `electromagnetic_ophelie_source.hpp` + `electromagnetic_ophelie_biot_savart.hpp`（excitation + B）
4. `electromagnetic_ophelie_postprocess.hpp` + `electromagnetic_ophelie_phi.hpp`（E/J/φ）
5. `TEAM7_probe_definitions.csv` + `TEAM7_Bz_A1_B1_reference_mT.csv`（comparison target）
6. terminal log section（section 6）
7. `OPHELIE_CODE_REVIEW_NEXT_STEPS.md`（if discussing roadmap）

---

## 10. common commands (complete paths, no `...`)

```bash
cd /home/yyc/SPHinXsysSYCL/build
ninja test_3d_ophelie_team7

# Relax → Reload
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_team7/bin/test_3d_ophelie_team7 \
 --native-stl --relax=1 --state_recording=0 --relax-log-every=100

# EM（your already-successful run）
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_team7/bin/test_3d_ophelie_team7 \
 --native-stl --reload=1 --no-phi --sigma=1e4 --ophelie-smoke --state_recording=0

# VTP field
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_team7/bin/test_3d_ophelie_team7 \
 --native-stl --reload=1 --no-phi --sigma=1e4 --state_recording=1
```

---

*documentversion： and repository native STL + reload fixafterStatusconsistent。*
