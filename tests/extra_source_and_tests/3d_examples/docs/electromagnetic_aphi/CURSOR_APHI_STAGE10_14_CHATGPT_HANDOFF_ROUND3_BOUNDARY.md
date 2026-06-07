# Stage 10.14 Round 3 Handoff — Boundary Strategy / Dummy vs Ghost (For ChatGPT)

> **User concern**: SphinxSys typically **does not use ghost particles for boundaries**; need to discuss with ChatGPT whether this round's boundary attempts deviate from framework convention.  
> **Cursor clarification**: This round P5 **implemented dummy support shell (Strategy C)**, **did not implement** addendum ghost/mirror particles (Strategy D). The repo also has an **earlier** `aphi_ghost_buffer_diva_diagnostic` (analytic ghost points supplement divA), a **different route** from P5 boundary diagnostics.

---

## 1. What This Round (Round 3) Actually Did

| Item | Status | Notes |
|----|------|------|
| P4/P4b Contact source-driven baseline | done | true two-body/three-body; coil excitation on left_body |
| P5 `boundary_width_scale` 1/2/3 | done | only enlarges outer air padding, does not change TEAM7 inner box |
| P5 **`dummy_shell`** | done | see §2; **not** SPH ghost boundary particle mechanism |
| P5 ghost/mirror (Strategy D) | **not done** | addendum explicitly "after B/C" |
| Contact interface spikes | done | warning only, `C_spike=10`, no fail gate |
| P2 thermal coupling | done | same-body spatial Joule + order-of-magnitude energy threshold; SYCL strict 5% not closed |
| P7 probe CSV | done in prior round | `aphi_probe_metric_helpers.h` |

---

## 2. Three Boundary Lines Compared (Discussion Focus)

### 2.1 Strategy B — Enlarged Air Padding (In Use, Uncontroversial)

- **Approach**: `AphiSourceDrivenEmSolveSpec::boundary_width_scale` enlarges `3*dp` outer bbox; TEAM7 physical box `[0,L]×[0,H]×[0,W]` unchanged.
- **Files**: `aphi_source_driven_em_solve_helpers.h`, `aphi_boundary_support_policy_diagnostic_helpers.h`
- **Observation**: conductor Joule ~0.0805 W, scale 1→3 relative change ~1e-4 order; baseline requires GMRES convergence.

### 2.2 Strategy C — Passive Dummy Support Shell (New This Round, Needs Discussion)

- **Approach**:
  - `AphiLhsTestBody(..., boundary_width, dummy_shell_width)` further expands bbox by `dummy_shell_width` (usually `2*dp`).
  - Still uses **regular Lattice real particles** filling extended box; **not** a separate ghost layer data structure.
  - `AssignZeroSigmaOutsidePhysicalBoxCK`: particles outside physical box get `sigma=0` (passive support, not in conductor/air TEAM7 region integration).
- **Naming**: code calls it `dummy_shell` / `dummy_shell_width_scale`, addendum calls it "passive dummy support shell".
- **Relation to SphinxSys**: essentially **more padding particles + zero material**, no framework-level ghost/mirror BC API invoked.
- **Test**: `test_3d_aphi_ck_boundary_support_policy_diagnostic` row 4 `policy=dummy_shell`.
- **Open discussion**: Under "no ghost for boundaries", is this **σ=0 shell particle** approach acceptable as Stage 10.14 diagnostic? Or keep Strategy B only?

### 2.3 Strategy D + Legacy Ghost divA Diagnostic (Not on P5 Production Path)

| Route | File | Purpose |
|------|------|------|
| addendum **ghost/mirror** | not implemented | explicit mirror/ghost support, different from dummy_shell |
| **`aphi_ghost_buffer_diva_diagnostic_helpers.h`** | existing research code | place **analytic ghost sample points** outside physical box to improve divA; `test_3d_aphi_ck_ghost_buffer_diva_diagnostic` |
| **`buildAnalyticGhostSamplePoints`** | same | not SPHBody ghost particles, host-side supplemental point set |

**Important**: do not conflate `ghost_buffer_diva` with P5 `dummy_shell`; discuss separately with user.

---

## 3. Key Code Locations (dummy_shell)

```text
aphi_lhs_test_helpers.h
  AphiLhsTestBody(..., dummy_shell_width=0)  // extend BoundingBox

aphi_source_driven_em_solve_helpers.h
  AphiSourceDrivenEmSolveSpec::dummy_shell_width_scale
  AssignZeroSigmaOutsidePhysicalBoxCK
  isInsideTeam7PhysicalBox()
  runSourceDrivenEmSolve() assigns dummy sigma before pipeline

aphi_boundary_support_policy_diagnostic_helpers.h
  enum AphiBoundarySupportPolicy::DummyShell
  runBoundaryPolicyRow(..., dummy_shell_width_scale=2.0)
```

---

## 4. Local Test Results (Summary)

```text
test_3d_aphi_ck_boundary_support_policy_diagnostic passed=1
  baseline / enlarged_air 2x,3x / dummy_shell — all converged=1, conductor_Joule~0.0805

test_3d_aphi_ck_contact_source_driven_heating_baseline passed=1
test_3d_aphi_ck_em_joule_thermal_one_way passed=1
```

---

## 5. Suggested Upload File List

### 5.1 Plans and Records (Required Reading)

| File | Path |
|------|------|
| Main plan | [`CURSOR_APHI_STAGE10_14_REAL_EM_HEATING_SOLVER_PLAN.md`](../../../../CURSOR_APHI_STAGE10_14_REAL_EM_HEATING_SOLVER_PLAN.md) |
| Contact/boundary addendum | [`CURSOR_APHI_STAGE10_14_CONTACT_BOUNDARY_ADDENDUM.md`](../../../../CURSOR_APHI_STAGE10_14_CONTACT_BOUNDARY_ADDENDUM.md) |
| Stage record | [`CURSOR_APHI_STAGE10_14_REAL_EM_HEATING_SOLVER_RECORD.md`](CURSOR_APHI_STAGE10_14_REAL_EM_HEATING_SOLVER_RECORD.md) |
| Round 2 handoff | [`CURSOR_APHI_STAGE10_14_CHATGPT_HANDOFF_ROUND2.md`](CURSOR_APHI_STAGE10_14_CHATGPT_HANDOFF_ROUND2.md) |
| **This round boundary handoff** | [`CURSOR_APHI_STAGE10_14_CHATGPT_HANDOFF_ROUND3_BOUNDARY.md`](CURSOR_APHI_STAGE10_14_CHATGPT_HANDOFF_ROUND3_BOUNDARY.md) |

### 5.2 Helpers Modified/Added This Round (Boundary + Related)

| File | Changes This Round |
|------|----------|
| [`aphi_boundary_support_policy_diagnostic_helpers.h`](../extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_boundary_support_policy_diagnostic_helpers.h) | +DummyShell policy |
| [`aphi_source_driven_em_solve_helpers.h`](../extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_source_driven_em_solve_helpers.h) | +dummy_shell_width_scale, AssignZeroSigmaOutsidePhysicalBoxCK |
| [`aphi_lhs_test_helpers.h`](../extra_src/shared/electromagnetic_dynamics/test_helpers/aphi_lhs_test_helpers.h) | AphiLhsTestBody extended bbox |
| [`aphi_em_observable_helpers.h`](../extra_src/shared/electromagnetic_dynamics/test_helpers/aphi_em_observable_helpers.h) | interface spike metrics |
| [`aphi_contact_source_driven_heating_baseline_helpers.h`](../extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_contact_source_driven_heating_baseline_helpers.h) | spike printing |
| [`aphi_em_joule_thermal_coupling_helpers.h`](../extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_em_joule_thermal_coupling_helpers.h) | same-body thermal coupling |
| [`aphi_ghost_buffer_diva_diagnostic_helpers.h`](../extra_src/shared/electromagnetic_dynamics/diagnostics/aphi_ghost_buffer_diva_diagnostic_helpers.h) | **unchanged**; for comparing another "ghost" research line |

### 5.3 Test Cases (This Round Related)

| Test | Path |
|------|------|
| P5 boundary diagnostic | [`test_3d_aphi_ck_boundary_support_policy_diagnostic/`](test_3d_aphi_ck_boundary_support_policy_diagnostic/) |
| P1 source-driven | [`test_3d_aphi_ck_source_driven_em_solve/`](test_3d_aphi_ck_source_driven_em_solve/) |
| P2 thermal coupling | [`test_3d_aphi_ck_em_joule_thermal_one_way/`](test_3d_aphi_ck_em_joule_thermal_one_way/) |
| P4 | [`test_3d_aphi_ck_contact_source_driven_heating_baseline/`](test_3d_aphi_ck_contact_source_driven_heating_baseline/) |
| P4b | [`test_3d_aphi_ck_three_body_contact_source_driven_baseline/`](test_3d_aphi_ck_three_body_contact_source_driven_baseline/) |
| (comparison) ghost divA | [`test_3d_aphi_ck_ghost_buffer_diva_diagnostic/`](test_3d_aphi_ck_ghost_buffer_diva_diagnostic/) |

### 5.4 Do Not Upload

- Unchanged large blocks such as full `aphi_block_jacobi_preconditioner_ck.hpp`
- Build artifacts `build/`

---

## 6. Questions for ChatGPT (Boundary Focus)

1. **Framework consistency**: Under SphinxSys **no ghost particles for boundaries**, should Stage 10.14 **forbid** `dummy_shell` (σ=0 expanded lattice particles) and keep only `boundary_width_scale` (Strategy B)?
2. **dummy_shell definition**: Current implementation is "extend bbox + sigma=0 outside box", **not** ghost/mirror BC. Is this acceptable **kernel support padding**, or still "non-physical particles" to remove?
3. **Relation to `aphi_ghost_buffer_diva`**: Should analytic ghost **point sets** for divA research stay **decoupled** from TEAM7 heating solver? Should P5 never reference that path?
4. **Strategy D**: Should ghost/mirror particles be **explicitly deferred to 10.15+** in Stage 10.14, or must there be a comparison before TEAM7 sign-off?
5. **Integration and observables**: dummy_shell particles already have `sigma=0`; conductor Joule integration still uses `team7ParticleInRegion` inside inner box. Need explicit exclusion of outside-box particles, or is current enough?
6. **P2 thermal + boundary**: Do boundary strategies affect "same-body Joule→thermal" acceptance? Can order-of-magnitude energy threshold serve as Stage 10.14 closure condition?
7. **Stage 10.14 closure**: Is minimum set = P1–P5 (with or without dummy_shell) + P4 + probe infrastructure, with standard TEAM7 reference comparison in 10.15?

---

## 7. One-Line Prompt for ChatGPT

```text
Stage 10.14 A-phi: P1–P7 baselines pass. Boundary work added Strategy B (padding scale)
and Strategy C (dummy_shell = extra lattice particles with sigma=0 outside [0,L]^3 box),
NOT SphinxSys ghost BC and NOT ghost/mirror Strategy D. Separate legacy ghost_buffer_diva
research exists. Please advise: keep dummy_shell for TEAM7 scaffold, or padding-only;
whether ghost_buffer_diva should stay decoupled; and Stage 10.14 closure criteria.
Files: CONTACT_BOUNDARY_ADDENDUM §5, RECORD, ROUND3_BOUNDARY handoff, P5 test + helpers listed above.
```
