# SPHinXsys-SPDB Capability Audit (SPDB)

**Date:** 2026-10-15
**Repo:** `/home/ubuntu/software/SPHinXsys-SPDB`
**Baseline:** `4e686ee35` on branch `feature/spdb-thermo-viscoplastic`
**Scope:** Geometry, solid particles, stress update, contact/damping, thermal, CK/SYCL solids, skeleton cases, CMake registration.

**Project constraints applied throughout:**
- Names containing `plastic` do **not** prove usability for high-temperature metal creep.
- Dynamic contact cases do **not** automatically verify quasi-static contact.
- Prefer stable existing solid formulations for P0/P1 elastic baselines.

---

## Summary matrix

| # | Capability | Real files / classes | CPU / SYCL | Verification status | Notes / gaps for SPDB |
|---|------------|----------------------|------------|---------------------|------------------------|
| 1 | Cylinder / box / boolean / custom Shape | See §1 | Geometry on host; reused by SYCL cases | Unit + example evidence | Cylinder+subtract practiced via mesh cylinders; parse cylinder+subtract not co-demoed |
| 2 | Solid particles from shapes | `ParticleGenerator<…, Lattice/Reload>` | Host generation; SYCL cases same API | Taylor bar CPU/SYCL | Sufficient for P0 lattice solids |
| 3 | Stress update + elastic/plastic | `Decomposed*` / `StructureIntegration*`; `NeoHookean` / `HardeningPlasticSolid` | Elastic + hardening plastic: CPU+CK; viscous plastic: CPU only | Taylor / twisting / cream-drop | **No high-T creep**; plastic ≠ creep |
| 4 | Contact / damping / relaxation | `ContactForce*`, `RepulsionForceCK`, `DampingWithRandomChoice` | CPU contact verified; CK contact present | Soft-body contact, self-contact, sliding, Taylor wall | Dynamic repulsion + velocity damping ≠ verified QS contact |
| 5 | Thermal / T-dependent solid | `IsotropicDiffusion*`, `DiffusionRelaxationCK` | Diffusion CPU+SYCL | Diffusion examples (often named temperature) | **No** `E(T)` / `σ_y(T)` / thermo-mechanical coupling |
| 6 | CK device material params (solids) | Nested `ConstituteKernel` + `DelegatedData` | Elastic + `HardeningPlasticSolid` | SYCL Taylor / beams | ViscousPlastic / nonlinear hardening / HJC lack CK kernels |
| 7 | Skeleton 3D cases | See §7 | CPU primary; SYCL Taylor twin | Regression / gtest | Soft-body contact best for bonding press; Taylor for plastic+wall |
| 8 | CMake 3D registration | `SUBDIRLIST` auto-discovery | CPU `3d_examples`; SYCL `tests_sycl/3d_examples` | Build system | New folder + `CMakeLists.txt` linking `sphinxsys_3d` |

---

## 1. Geometry

### 1.1 Capability table

| Capability | Real files / classes | CPU / SYCL | Verification | Notes / gaps for SPDB |
|------------|----------------------|------------|--------------|------------------------|
| 3D cylinder (analytic) | `GeometricCylinder`, `GeometricShapeCylinder` — `src/for_3D_build/geometries/geometric_element_and_shape_3d.h` | Host Shape; SYCL examples reuse | Unit: `tests/unit_tests_src/for_3D_build/geometrices/test_3d_geometric_cylinder/` | Local axis = **x**; orient with `Transform` / `FromTwoVectors(UnitX, axis)` |
| 3D cylinder (mesh) | `TriangleMeshShapeCylinder` — `src/for_3D_build/geometries/triangle_mesh_shape.h` | Host | Taylor bar; SYCL mixed Poiseuille | Preferred when boolean subtract needed |
| 3D cylinder (SDF) | `SDFCylinder` — `src/for_3D_build/geometries/sdf_primitive.h` | Host SDF | SYCL `test_3d_sdf_shape_sycl` | Optional for smooth fillets later |
| Box / cube (analytic) | `GeometricBox` — `src/shared/geometries/geometric_element.h`; `GeometricShapeBox` — `src/shared/geometries/geometric_shape.h` | Host | Unit geometric shape; many wall examples | Blind square slot: subtract box from cylinder face |
| Box (mesh) | `TriangleMeshShapeBrick` — `triangle_mesh_shape.h` | Host | STFB / walls | Same role as analytic box |
| Boolean add / subtract | `GeometricOps::{add,sub,…}` · `BinaryShapes` · `ComplexShape` — `src/shared/geometries/base_geometry.h`, `complex_geometry.h` | Host boolean | Unit `test_complex_shape`; STFB walls | **Only add+sub** on Binary/Complex; intersect/sym_diff not implemented there |
| SDF boolean subtract | `SDFSubtraction`, `SDFSmoothSubtraction` — `src/for_3D_build/geometries/sdf_operation.h` | Host | SDF examples | Useful for smooth cavity edges later |
| Custom Shape / distance / normal | Base `Shape`: pure virtual `checkContain`, `findClosestPoint`, `findBounds`; derived `findSignedDistance` (interior **negative**), `findNormalDirection` (outward) — `base_geometry.h` | Host | `test_complex_shape` normals | Subclass `Shape` or compose via `ComplexShape` |
| Aggregator include | `src/for_3D_build/geometries/all_geometries.h` | — | — | 3D geometry entry |

### 1.2 Example usage

```cpp
// Analytic cylinder (axis via Transform; local cylinder axis = UnitX)
Transform T(Rotation3d(Eigen::Quaternion<Real>::FromTwoVectors(Vec3d::UnitX(), Vec3d::UnitZ())),
            Vec3d(0, 0, 0.01));
GeometricShapeCylinder cyl(T, /*radius*/0.02, /*halflength*/0.01);

// Box / cube half-sizes
GeometricShapeBox slot(Transform(slot_center), Vec3d(0.005, 0.005, 0.005), "BlindSlot");

// Boolean difference (recommended for outer − cavity)
ComplexShape lower("LowerCylinderWithSlot");
lower.add<GeometricShapeCylinder>(T, 0.02, 0.01);
lower.subtract<GeometricShapeBox>(Transform(slot_center), Vec3d(0.005, 0.005, 0.005));

// Mesh path used in existing cylinder+subtract demos
auto wall = makeShared<ComplexShape>("Wall");
wall->add<TriangleMeshShapeCylinder>(axis, R_out, L/2, res, translation, "Outer");
wall->subtract<TriangleMeshShapeCylinder>(axis, R_in, L/2, res, translation, "Inner");
```

**SPDB gap:** No in-tree demo of `GeometricShapeCylinder` + `subtract` together; practice evidence for cylinder subtract is `TriangleMeshShapeCylinder` in `tests/tests_sycl/3d_examples/test_3d_mixed_poiseuille_flow_sycl/mixed_poiseuille_flow.cpp`. Analytic path should still work via `ComplexShape` API but needs SPDB case verification.

---

## 2. Particle generation for solids from shapes

| Capability | Real files / classes | CPU / SYCL | Verification | Notes / gaps for SPDB |
|------------|----------------------|------------|--------------|------------------------|
| Lattice fill from Shape | `ParticleGenerator<BaseParticles, Lattice>`, `GeneratingMethod<Lattice>` — `src/shared/particle_generator/particle_generator_lattice.h` | Host; used in SYCL solid cases | Taylor bar, soft-body contact, twisting column | Uses body `initial_shape_` / `checkContain` |
| Reload from file | `ParticleGenerator<BaseParticles, Reload>` — `base_particle_generator.h` | Host | Taylor bar `--relax` / `--reload` | Good after particle relaxation |
| Unstructured mesh | `ParticleGenerator<BaseParticles, UnstructuredMesh>` — `particle_generator_mesh.h` | Host | Mesh-driven cases | Not required for P0 cylinders |
| Lattice blockers / inserts | OrientedBox blockers; `GeometricShapeBox` / `GeometricShapeCylinder` inserts in lattice generator | Host | API present in header | Optional refinement inserts |

```cpp
body.defineAdaptationRatios(1.15, 1.0);
body.defineMatterMaterial<NeoHookeanSolid>(rho0, E, nu);
body.generateParticles<BaseParticles, Lattice>();
// or after relaxation:
body.generateParticles<BaseParticles, Reload>(body.Name());
```

---

## 3. Solid stress update, elastic / plastic / creep

### 3.1 Stress update (prefer stable)

| Capability | Real files / classes | CPU / SYCL | Verification | Notes / gaps for SPDB |
|------------|----------------------|------------|--------------|------------------------|
| Decomposed elastic 1st half | `solid_dynamics::DecomposedIntegration1stHalf` — `src/shared/particle_dynamics/solid_dynamics/elastic_dynamics.h` | CPU | Twisting / shell / beams | **Preferred CPU elastic** (hourglass-resistant) |
| Decomposed plastic 1st half | `solid_dynamics::DecomposedPlasticIntegration1stHalf` — `inelastic_dynamics.h` | CPU | Taylor bar | Preferred CPU plastic |
| Classic PK2 / Cauchy / Kirchhoff | `Integration1stHalfPK2`, `Integration1stHalf`, … + `Integration2ndHalf` — `elastic_dynamics.h` | CPU | Many solids | Alternative formulations |
| CK structure integration | `solid_dynamics::StructureIntegration1stHalf<…, MaterialType, …>`, `StructureIntegration2ndHalf`, `AcousticTimeStepCK` — `src/shared/shared_ck/particle_dynamics/solid_dynamics/structure_dynamics.h`; aggregate `all_solid_dynamics_ck.h` | CK/SYCL | SYCL Taylor, twisting column, oscillating beam | **Preferred device path** |
| Aggregates | `all_solid_dynamics.h`, `all_solid_dynamics_ck.h` | — | — | Include hubs |

### 3.2 Elastic materials

| Capability | Real files / classes | CPU / SYCL | Verification | Notes / gaps for SPDB |
|------------|----------------------|------------|--------------|------------------------|
| Linear elastic | `LinearElasticSolid` — `src/shared/materials/elastic_solid.h` (+ `.hpp`) | CPU + `ConstituteKernel` | Shells / beams | Good synthetic P0 material |
| SVK | `SaintVenantKirchhoffSolid` | CPU + CK | SYCL oscillating / anisotropic beam | Finite strain elastic |
| Neo-Hookean | `NeoHookeanSolid` | CPU + CK | CPU+SYCL twisting column; soft-body contact | Stable soft-solid baseline |
| Incompressible NH | `NeoHookeanSolidIncompressible` | CPU (comment: with `DecomposedIntegration1stHalf`) | Limited | No independent CK path called out |
| Orthotropic / muscle / FENE | `OrthotropicSolid`, `Muscle`, … | Mostly CPU | Specialty | Not SPDB metal path |

### 3.3 Plastic / “viscoplastic” (not creep)

| Capability | Real files / classes | CPU / SYCL | Verification | Notes / gaps for SPDB |
|------------|----------------------|------------|--------------|------------------------|
| Plastic base | `PlasticSolid` : `NeoHookeanSolid` — `src/shared/materials/inelastic_solid.h` | — | — | Abstract; **not creep** |
| J2 linear hardening | `HardeningPlasticSolid` (+ nested `ConstituteKernel`) | CPU + CK/SYCL | CPU+SYCL Taylor bar | Only solid plastic with full CK evidence |
| Nonlinear hardening | `NonLinearHardeningPlasticSolid` | CPU only (virtual `ElasticLeftCauchy`) | `test_2d_stretching` | No CK kernel |
| Herschel–Bulkley viscous plastic | `ViscousPlasticSolid` | **CPU only** | `test_2d_viscous_cream_drop` (paste) | Rate-dependent yield; **≠ high-T metal creep**; no CK |
| HJC | `HJCSolid` + `HJCIntegration1stHalf` | CPU | `test_3d_taylor_bar_hjc` | Concrete impact; irrelevant to SPDB creep |
| Continuum J2 / DP | `J2Plasticity`, `PlasticContinuum` — `general_continuum.h` | CK for continuum UL | Granular / continuum | Not solid bonding constitutive |

### 3.4 Creep / high-T viscoplasticity

| Capability | Real files / classes | CPU / SYCL | Verification | Notes / gaps for SPDB |
|------------|----------------------|------------|--------------|------------------------|
| Norton / Arrhenius / creep | **None found** | — | — | **Must implement** for 20Cr13 hold; do not reuse cream-drop viscous plastic as creep |

---

## 4. Contact / surface forces / damping / relaxation

| Capability | Real files / classes | CPU / SYCL | Verification | Notes / gaps for SPDB |
|------------|----------------------|------------|--------------|------------------------|
| Contact relations | `SurfaceContactRelation` — `src/shared/body_relations/contact_body_relation.h`; `SelfSurfaceContactRelation` — `inner_body_relation.h` | Host topology | Soft-body / self-contact / Taylor | Required wiring for contact |
| Contact density factor | `ContactFactorSummation`, `SelfContactFactorSummation` — `…/contact_dynamics/repulsion_factor_summation.h` | CPU | Soft-body contact | Precedes force |
| Solid–solid repulsion | `ContactForce` (`RepulsionForce<Contact<>>`) — `…/contact_dynamics/contact_repulsion.h` | CPU | `test_3d_muscle_soft_body_contact` | Dynamic penalty-style contact |
| Self-contact | `SelfContactForce` | CPU | `test_3d_self_contact` | Coil / folding |
| Wall contact | `ContactForceFromWall`, `ContactForceToWall` | CPU | Sliding; Taylor custom wall force | Fixtures / platens |
| Wall friction | `PairwiseFrictionFromWall` — `contact_friction.h` | CPU | Shell collision | Optional |
| CK repulsion | `RepulsionForceCK` — `src/shared/shared_ck/…/contact_dynamics/repulsion_force.h` | CK/SYCL API | Present in tree; fewer dedicated 3D QS demos | Device contact path exists |
| Velocity damping (QS approach) | `DampingWithRandomChoice` + `DampingPairwiseInner<Vec3d, FixedDampingRate>` — `…/dissipation_dynamics/particle_dynamics_dissipation.h` | CPU | Soft-body / self-contact comments “quasi-static” | Damps dynamics toward equilibrium; **does not certify QS contact accuracy** |
| Particle distribution relaxation | `relax_dynamics::RelaxationStepInner`, `RandomizeParticlePosition` — `…/relax_dynamics/` | CPU (+ SYCL relaxation examples) | Taylor `--relax`; SYCL particle_relaxation | Initial packing, not contact |
| Loading (not contact) | `ForceInBodyRegion`, springs, surface pressure — `loading_dynamics.h` | CPU | Various | Apply bonding pressure loads |

```cpp
SurfaceContactRelation contact_ab(body_a, {&body_b});
InteractionDynamics<solid_dynamics::ContactFactorSummation> factor(contact_ab);
InteractionWithUpdate<solid_dynamics::ContactForce> force(contact_ab);
DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec3d, FixedDampingRate>>>
    damp(0.2, inner_relation, "Velocity", physical_viscosity);
```

**SPDB gap:** Hold-pressure bonding needs separate quasi-static / long-hold contact verification; existing cases are short dynamic contact with optional velocity damping.

---

## 5. Thermal diffusion / temperature-dependent material

| Capability | Real files / classes | CPU / SYCL | Verification | Notes / gaps for SPDB |
|------------|----------------------|------------|--------------|------------------------|
| Isotropic / directional diffusion | `IsotropicDiffusion`, `LocalIsotropicDiffusion`, `DirectionalDiffusion` — `src/shared/materials/diffusion_reaction.h` | CPU | Optimization + diffusion examples | Can represent temperature as species |
| CK diffusion relaxation | `DiffusionRelaxationCK` — `src/shared/shared_ck/…/diffusion_reaction_dynamics/diffusion_dynamics_ck.h` | CK/SYCL | `tests/tests_sycl/2d_examples/test_2d_diffusion_NeumannBC_sycl/` etc. | Thermal field OK as diffusion |
| T-dependent solid moduli / yield | **None** | — | — | No `E(T)`, `σ_y(T)`, thermal expansion, or thermo-mechanical constitutive |
| Electro-thermal / Joule heat | **None in solid continuum path** | — | — | Phase-2 SPDB (pulse current); out of P0 |

**SPDB gap:** Diffusion module is **decoupled** from solid stress update. Coupling + T-dependent metal laws are new work.

---

## 6. CK / SYCL device material parameter support (solids)

| Capability | Real files / classes | CPU / SYCL | Verification | Notes / gaps for SPDB |
|------------|----------------------|------------|--------------|------------------------|
| Elastic kernel params | `LinearElasticSolid` / `SVK` / `NeoHookeanSolid`::`ConstituteKernel` — `elastic_solid.hpp` | Scalars copied into kernel | SYCL beams / twisting | Ready for device elastic |
| Hardening plastic device state | `HardeningPlasticSolid::ConstituteKernel` with `DelegatedData` for plastic strain / hardening — `inelastic_solid.hpp`; variables via `sphinxsys_variable.h` | CK/SYCL | `test_3d_taylor_bar_sycl` | Only verified solid plastic on device |
| ViscousPlastic / NonLinearHardening / HJC on device | No solid `ConstituteKernel` for these | CPU only | Cream-drop / HJC CPU | Cannot assume SYCL for creep-like viscous plastic |
| Structure dynamics CK | `StructureIntegration1stHalf` uses material `ConstituteKernel` | SYCL | Taylor / twisting | Prefer this over inventing new CPU-only loops for GPU path |

---

## 7. Recommended skeleton examples

| Priority | Path | Demonstrates | CMakeLists | Why for SPDB |
|----------|------|--------------|------------|--------------|
| **1** | `tests/3d_examples/test_3d_muscle_soft_body_contact/` (`muscle_soft_body_contact.cpp`) | Two NeoHookean solids, bidirectional `SurfaceContactRelation` + `ContactForce`, velocity damping, no Simbody | `…/test_3d_muscle_soft_body_contact/CMakeLists.txt` | Closest to upper/lower workpiece press-fit contact |
| **2** | `tests/3d_examples/test_3d_taylor_bar/` (`taylor_bar.h/.cpp`) | `TriangleMeshShapeCylinder`, `HardeningPlasticSolid`, Lattice/Reload, wall contact, LABEL `solid_dynamics; contact` | `…/test_3d_taylor_bar/CMakeLists.txt` | Cylinder solid + plastic + wall; SYCL twin available |
| **3** | `tests/3d_examples/test_3d_self_contact/` (`3d_self_contact.cpp`) | Relax/reload pipeline, self + wall contact, QS damping comments | `…/test_3d_self_contact/CMakeLists.txt` | Geometry/packing workflow for complex cavities |
| Alt | `tests/3d_examples/test_3d_sliding/` | `ContactForceFromWall` + gtest | case `CMakeLists.txt` | Rigid support / platen contact |
| Alt SYCL | `tests/tests_sycl/3d_examples/test_3d_taylor_bar_sycl/` | CK plastic + contact labels | `…/CMakeLists.txt` | Device baseline after CPU elastic case |

### Concrete recommended skeleton file paths

1. **Soft-body contact (primary contact skeleton)**
   - Source: `/home/ubuntu/software/SPHinXsys-SPDB/tests/3d_examples/test_3d_muscle_soft_body_contact/muscle_soft_body_contact.cpp`
   - CMake: `/home/ubuntu/software/SPHinXsys-SPDB/tests/3d_examples/test_3d_muscle_soft_body_contact/CMakeLists.txt`

2. **Taylor bar (cylinder + plastic + wall)**
   - Source: `/home/ubuntu/software/SPHinXsys-SPDB/tests/3d_examples/test_3d_taylor_bar/taylor_bar.cpp`
   - Header: `/home/ubuntu/software/SPHinXsys-SPDB/tests/3d_examples/test_3d_taylor_bar/taylor_bar.h`
   - CMake: `/home/ubuntu/software/SPHinXsys-SPDB/tests/3d_examples/test_3d_taylor_bar/CMakeLists.txt`
   - SYCL twin CMake: `/home/ubuntu/software/SPHinXsys-SPDB/tests/tests_sycl/3d_examples/test_3d_taylor_bar_sycl/CMakeLists.txt`

3. **Self-contact / relax-reload (geometry pipeline)**
   - Source: `/home/ubuntu/software/SPHinXsys-SPDB/tests/3d_examples/test_3d_self_contact/3d_self_contact.cpp`
   - CMake: `/home/ubuntu/software/SPHinXsys-SPDB/tests/3d_examples/test_3d_self_contact/CMakeLists.txt`

---

## 8. CMake registration pattern (new 3D case)

### 8.1 Discovery chain (CPU)

```
CMakeLists.txt (root)
  → add_subdirectory(tests)                    # when SPHINXSYS_BUILD_TESTS
tests/CMakeLists.txt
  → ADD_SUBDIRECTORY(3d_examples)              # SPHINXSYS_3D && SPHINXSYS_BUILD_3D_EXAMPLES
                                               # skipped when SPHINXSYS_USE_SYCL
tests/3d_examples/CMakeLists.txt               # SUBDIRLIST auto-discovery
  → foreach subdir with CMakeLists.txt: add_subdirectory(${subdir})
```

Parent auto-discovery (exact pattern):

```cmake
# tests/3d_examples/CMakeLists.txt
SUBDIRLIST(SUBDIRS ${CMAKE_CURRENT_SOURCE_DIR})
foreach(subdir ${SUBDIRS})
    if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/CMakeLists.txt)
        add_subdirectory(${subdir})
    endif()
endforeach()
```

**No manual name entry** in the parent list is required for CPU `3d_examples`.

### 8.2 SYCL twin location

When `SPHINXSYS_USE_SYCL` is on, `tests/CMakeLists.txt` does **not** add `3d_examples`; instead use:

```
tests/tests_sycl/3d_examples/CMakeLists.txt   # same SUBDIRLIST pattern
```

### 8.3 Optional user examples

`SPHINXSYS_BUILD_USER_EXAMPLES` (default **OFF**) → `tests/user_examples/` with the same `SUBDIRLIST` auto-discovery.

### 8.4 Minimal new case `CMakeLists.txt` (copy from soft-body contact)

Create directory: `tests/extra_source_and_tests/3d_examples/test_3d_spdb_<name>/`
Reusable sources go in `tests/extra_source_and_tests/extra_src/shared` or `for_3D_build`, and the case links `extra_sources_3d`.

```cmake
STRING(REGEX REPLACE ".*/(.*)" "\\1" CURRENT_FOLDER ${CMAKE_CURRENT_SOURCE_DIR})
PROJECT("${CURRENT_FOLDER}")

SET(EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}/bin/")
SET(BUILD_INPUT_PATH "${EXECUTABLE_OUTPUT_PATH}/input")
file(MAKE_DIRECTORY ${BUILD_INPUT_PATH})
# optional: file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/regression_test_tool/ DESTINATION ${BUILD_INPUT_PATH})

aux_source_directory(. DIR_SRCS)
add_executable(${PROJECT_NAME} ${DIR_SRCS})
target_link_libraries(${PROJECT_NAME} extra_sources_3d)
set_target_properties(${PROJECT_NAME} PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${EXECUTABLE_OUTPUT_PATH}")

add_test(NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME}
    WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH})
```

Relax/reload variant (Taylor / self-contact style): add `--relax` / `--reload` tests and `DEPENDS` between them; optionally `LABELS "solid_dynamics; contact"`.

---

## SPDB readiness verdict (P0-relevant)

| Need | Status |
|------|--------|
| Parametric cylinders + box blind slot (boolean) | **Available** via `ComplexShape` + cylinder/box shapes |
| Lattice solid particles | **Available** |
| Stable elastic compression baseline | **Available** (`NeoHookean` / `LinearElastic` + `DecomposedIntegration1stHalf` or CK `StructureIntegration1stHalf`) |
| Solid–solid contact + damping toward equilibrium | **Available** as **dynamic** contact; QS accuracy **unverified** |
| High-T metal creep / viscoplastic hold | **Missing** — do not treat `ViscousPlasticSolid` as creep |
| Temperature-dependent constitutive + thermo-mechanical coupling | **Missing** |
| Pulse electro-thermal | **Missing** (later phase) |
| CK device plastic | **Partial** — `HardeningPlasticSolid` only |

**Recommended path:** new case under `tests/extra_source_and_tests/3d_examples/test_3d_spdb_*`, linking `extra_sources_3d`. Reusable code goes in `extra_src/shared` or `extra_src/for_3D_build`. Geometry from `GeometricShapeCylinder` / `GeometricShapeBox` + `ComplexShape::subtract`; synthetic elastic material (explicitly not 950 °C 20Cr13); coarse lattice; VTK + volume checks; CPU first.
