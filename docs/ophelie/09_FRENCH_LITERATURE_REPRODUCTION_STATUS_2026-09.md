# French Literature Reproduction — Current Status (2026-09)

> **Purpose:** Snapshot of branch `feature/EM-Biot-Savart` for roadmap discussion (ChatGPT / advisor).  
> **Not a claim of completion.** Midterm goal remains an auditable SPH–EM–thermal–stirring *prototype*, not strict Jacoutot quantitative reproduction.  
> **Primary paper:** Jacoutot et al., *Chem. Eng. Sci.* **63** (2008) 2391–2401 — *Numerical modeling of coupled phenomena in a mechanically stirred molten-glass bath heated by induction*.  
> **Related:** earlier Jacoutot / OPHELIE natural-convection work; thesis σ(T) tables under `docs/ophelie/reference/jacoutot_sigma_t/`.  
> **Repo snapshot:** `89d300b64` (2026-09-12).  
> **Companion audits:** [`OPHELIE_CURRENT_CODE_AUDIT.md`](OPHELIE_CURRENT_CODE_AUDIT.md), [`01_MASTER_DEVELOPMENT_TIMELINE.md`](01_MASTER_DEVELOPMENT_TIMELINE.md), [`06_RH200_GLASS_EM_STIRRING.md`](06_RH200_GLASS_EM_STIRRING.md), stirring README `tests/.../test_3d_ophelie_french_stirring_em/README.md`.

---

## 1. One-sentence verdict

**Glass-domain complex edge-flux EM + power calibration + one-shot Eulerian Q → SPH stirring is runnable; full Jacoutot 2008 CES reproduction (Fig. 4 iterative EM↔thermo coupling, cold-crucible carriers, skull, field-figure sign-off) is not done.**

Estimated maturity (discussion numbers only):

| Layer | ~Completion | Meaning |
|-------|-------------|---------|
| EM kernel (glass, φ + edge-flux, P_recon calibrate) | **70–80%** | Solves and books power; pointwise field quality still open |
| Alignment with *this* stirred paper | **40–50%** | Closest case: `french_stirring_em` (60 kW, 10 rpm) |
| Paper-grade quantitative reproduction | **20–30%** | Missing iterative coupling, A_metal, skull, Fig. benchmarks |

---

## 2. What the paper actually requires

### 2.1 Physics / numerics (CES 2008)

| Item | Paper |
|------|--------|
| Process | Cold-crucible vitriﬁcation, direct induction, **mechanical stirrer** |
| Supply side | ~400 kW generator @ **300 kHz** (not equal to glass absorbed power) |
| Stirred case absorbed power | **~60 kW** in glass; stirrer **10 rpm** |
| EM model | **OPHELIE** (integral / surface methods); `A` from Biot–Savart; `j = −σ(∇V + … A)`; `Q ∝ \|j\|²/σ` |
| Cold crucible | Segmented cooled walls; surface currents; skull layer of frozen glass |
| Properties | Strong **T-dependence**: μ, cp, k, **σ** (Fig. 2) |
| Coupling (Fig. 4) | Supervisor iterates **FLUENT 3D thermo-hydro** ↔ **OPHELIE 2D axisymmetric EM**: azimuthal-average σ → EM → rotate `Q(r,z)` back to 3D; under-relax σ |
| Thermo-hydro | Laminar 3D unsteady Newtonian; **Boussinesq**; wall / free-surface heat losses |
| Macro indicators | Power number **Np**, Nusselt **Nu** (adapted for glass) |
| Results horizon | Coupling to ~**1000 s**; skull shape after convergence |

### 2.2 What “reproduction” should mean (proposal for discussion)

Recommend splitting acceptance into three gates (do not mix them):

1. **Gate A — EM operator trust:** closed-form / reduced cases (cylinder axial B, skin-depth, French r–z probes); power book + local J/E quality.  
2. **Gate B — Literature coupling architecture:** something that *behaves like* Fig. 4 (periodic EM update from σ(T), Q in lab frame), even if still SPH and still reduced geometry.  
3. **Gate C — Quantitative paper match:** Fig. 5–8 style `Q` / T / velocity / skull; Np, Nu; long-time heat balance without ad-hoc loss scaling.

Current branch is strongest on **Gate A infrastructure** and a **partial Gate B demo**; **Gate C is open**.

---

## 3. What we have on this branch

### 3.1 Production EM path

```text
Multiloop filament coil
  → Biot–Savart A_coil (/ B_coil)
  → complex edge-flux (φ_r, φ_i + edge EMF drop)
  → E/J reconstruction → Q = ½ σ (|E_r|² + |E_i|²)
  → glass volume P_recon
  → coil current scale: I ← I √(P_target / P_raw)   [literature mode: no post-hoc field_scale]
```

Key docs: `03_EDGE_FLUX_PRODUCTION.md`, `FRENCH_LITERATURE_MODE.md`, `08_EDGE_FLUX_UNDIRECTED_POWER_CLOSURE.md`  
(calibration must use **`P_recon`**, not graph Laplace energy).

### 3.2 Main executables (French / RH200 track)

| Executable | Role |
|------------|------|
| `test_3d_ophelie_french_reduced` | Reduced cylinder EM regression; `--literature-mode` @ ~50 kW |
| `test_3d_ophelie_french_natural_em` | Natural geometry EM + power gate |
| `test_3d_ophelie_french_complex_joule_to_heat_one_way` | One-way EM → heat; σ(T) options |
| `test_3d_ophelie_french_natural_convection_frozen_q` | Frozen particle Q + WCSPH + Boussinesq (no stirring transport of Q) |
| `test_3d_ophelie_french_stirring_em` | **Closest to CES 2008 stirred case:** EM (+ optional A_ind Picard) → Euler Q grid → SPH + paddle @ 10 rpm, target **60 kW** |
| `test_3d_ophelie_rh200_glass_em_stirring` | Midterm **client geometry** demo (STL); em-grid; not strict French CAD |
| `test_3d_ophelie_analytic_induction_simple` | Cylinder / thin plate / skin-depth / σ–f sweep / French probe template |
| `test_3d_ophelie_french_aind_diagnostic` / `*_self_induction_picard*` | A_glass / Picard (still experimental acceptance) |
| TEAM7 cases | High-σ side-track; **not** French delivery gate |

### 3.3 `french_stirring_em` — closest literature demo

Documented intent (see case README):

- Phase A: complex edge-flux (+ self-induction Picard by default) → calibrate to **60 kW** → CIC deposit Q on **fixed Euler grid**.  
- Phase B: WCSPH melt + kinematic paddle + Boussinesq + French Robin/radiation BC; resample Q each step (lab-frame source).  
- Geometry / ops: EREBUS-like melt R≈0.25 m; coil 7 turns, 300 kHz; **10 rpm**; T₀≈1473 K; σ law from thesis Fig. IV.10 at reference T.

**Documented simplifications (already acknowledged in README):**

- No skull → raw literature heat-loss coeffs over-cool; production uses `--balance-heat-loss` scale (~0.4) so loss ≈ 60 kW at T₀.  
- No SPH thermal conduction in long stir (same as frozen-Q production choice).  
- Paddle adiabatic.  
- EM not re-solved every thermo cycle (not Fig. 4).  
- Rotor VTP Position NaN fixed on this branch (`89d300b64`, host paddle kinematics).

Recent run (server output, 2026-08-30 UTC): ~600 s physical / ~100 rev completed in one wall-clock budget; monitor + VTP available; subsampled soft-copyright pack from **new** `output/` (Rotor Position finite).

### 3.4 Analytic / validation side (Gate A)

| Case | Status |
|------|--------|
| Cylinder + axial B₀ | Smoke/sweep **pass** at relaxed gates (J_rel ~0.3) |
| Thin plate + axial B₀ | **Pass** |
| Plate skin-depth B_x | Gate-pass only with **relaxed** J_rel &lt; 0.90 (physics quality ~0.87 — open-boundary / resolution / complex A_imag issues) |
| French r–z probe template | 30-point CSV committed; solve at stable σ/f; external COMSOL/OPHELIE ref **TBD** |
| TEAM7 L2 / no-flux | Side-track; incomplete |

---

## 4. Gap matrix vs Jacoutot CES 2008

| Paper requirement | Current status | Gap severity |
|-------------------|----------------|--------------|
| Glass volume EM of OPHELIE *form* | Edge-flux + A_coil (+ experimental A_ind) | Medium (method differs: SPH volume vs OPHELIE surface/integral) |
| A_total = A_coil + A_glass + A_metal | A_coil main; A_glass experimental; **no A_metal** | **High** for cold crucible |
| Segmented cold-crucible / skin carriers | Not implemented | **High** |
| Fig. 4 iterative 3D↔2D EM–thermo coupling | One-shot (or Picard) EM → frozen Euler Q | **High** (architecture) |
| Azimuthal σ average + Q(r,z) rotate | Not implemented | **High** |
| Full μ(T), cp(T), k(T), σ(T) in stir | σ(T) tables exist; stir mostly const props @ 1473 K | Medium |
| Skull layer | Absent; heat balance patched by loss scaling | **High** for long-time T |
| Field / T figure quantitative match | No systematic Fig. 5–8 campaign | **High** for “reproduction” claim |
| Np, Nu reporting | Not paper-style | Medium |
| Power book ~50–60 kW | Yes (P_recon + current calibrate) | Low (book ≠ field shape) |
| 10 rpm stirred SPH with lab-frame Q | Yes (`french_stirring_em`) | Low–medium |

---

## 5. Explicit non-claims (do not oversell)

1. **`literature_passed=1`** on reduced/natural EM ≠ full OPHELIE or CES 2008 reproduction.  
2. **P_scaled / P_recon ≈ 50–60 kW** ≠ correct spatial Joule map or skin-depth profile.  
3. **RH200** is midterm delivery geometry, **not** French cold-crucible CAD.  
4. **`french_stirring_em` with `--balance-heat-loss`** is an engineering heat-budget fix, **not** skull physics.  
5. **A_ind Picard “on”** ≠ validated self-induction against OPHELIE reference.  
6. **Skin-depth case gate-pass** ≠ skin-depth physics quality closed.

---

## 6. Candidate development directions (for ChatGPT debate)

Use this section as a menu; pick **one primary north star** for the next 4–8 weeks.

### Direction A — Trust the EM field first (Gate A)

**Goal:** Make J/E believable before more coupling.

- Tighten skin-depth plate (finer dp, boundary closure, complex A_imag).  
- French probe CSV vs external OPHELIE/COMSOL once geometry ready.  
- Optional: TEAM7 as high-σ regression (do not block French).  
- Pointwise J·E residual audit (Task B).

**Pros:** Reduces risk of coupling garbage-in/garbage-out.  
**Cons:** Less visible “stirring paper” progress; client may want figures sooner.

### Direction B — Weak Fig. 4 coupling (Gate B light)

**Goal:** Re-solve EM every N thermo steps from σ(T) (or azimuthally averaged σ), refresh Euler Q; keep SPH 3D.

**Pros:** Directly addresses largest architectural gap vs CES 2008.  
**Cons:** Costly; needs stable σ under-relaxation; still no A_metal/skull.

### Direction C — Thermal fidelity (skull / BC)

**Goal:** Replace `--balance-heat-loss` with skull shell or effective conductivity layer; free-surface crust option.

**Pros:** Long-time T meaningful vs paper ~1000 s.  
**Cons:** Hard multiphase/solidification; may still heat with wrong Q shape if Gate A weak.

### Direction D — Cold-crucible carriers (A_metal)

**Goal:** Segmented wall / bottom thin conductors in EM.

**Pros:** Closer to true OPHELIE process EM.  
**Cons:** Large EM surface/volume modeling effort; may be out of midterm scope.

### Direction E — Client / demo polish (RH200 + soft-copyright figures)

**Goal:** Stable movies, energy CSVs, Np-like indicators on RH200 or french_stirring; not paper L2 match.

**Pros:** Delivery and IP figures.  
**Cons:** Does not close reproduction claim.

### Suggested default package (discussion strawman)

1. **Primary:** Direction **A** (skin + French probes) until interior J_rel is in the same ballpark as cylinder axial-B.  
2. **Parallel light:** Direction **B** prototype (EM refresh every N steps on french_stirring or natural) with fixed σ law and power re-calibrate.  
3. **Defer:** D until A+B show stable Q(r,z).  
4. **Keep:** E only for client deliverables, not as science gate.

---

## 7. Questions to ask ChatGPT / advisor

1. For SPH midterm, is **Gate B-light** (periodic 3D EM refresh) enough, or must we emulate paper’s **2D axisymmetric OPHELIE** projection exactly?  
2. Priority: **field validation (A)** vs **stirring figures (E)** vs **iterative coupling (B)** given remaining calendar?  
3. Is **skull** mandatory for any CES-like claim, or is calibrated effective BC acceptable with explicit disclaimer?  
4. Should **A_metal** be in-scope for this PhD/project phase, or declared out-of-scope with coil+glass-only EM?  
5. Accept skin-depth as **diagnostic with relaxed gate**, or treat J_rel≈0.87 as a **blocking** EM defect?  
6. How to define a single **“French reproduction milestone”** sentence that is honest for midterm vs final thesis?

---

## 8. Key paths (quick index)

```text
docs/ophelie/09_FRENCH_LITERATURE_REPRODUCTION_STATUS_2026-09.md   ← this file
docs/ophelie/OPHELIE_CURRENT_CODE_AUDIT.md
docs/ophelie/01_MASTER_DEVELOPMENT_TIMELINE.md
docs/ophelie/04_FRENCH_REDUCED_AND_THERMAL.md
docs/ophelie/06_RH200_GLASS_EM_STIRRING.md
docs/ophelie/FRENCH_LITERATURE_MODE.md
docs/ophelie/08_ANALYTIC_INDUCTION_SIMPLE.md          (if present)
tests/.../test_3d_ophelie_french_stirring_em/README.md
tests/.../test_3d_ophelie_analytic_induction_simple/
Paper PDF (repo root): Numerical modeling of coupled phenomena in a mechanically stirred...
```

---

## 9. Change log

| Date | Note |
|------|------|
| 2026-09-12 | Initial status pack for ChatGPT / advisor roadmap discussion on `feature/EM-Biot-Savart`. |
| 2026-09-14 | Round-1 coupling implementation started: see [`10_THERMO_EM_COUPLING_ROUND1.md`](10_THERMO_EM_COUPLING_ROUND1.md). |
| 2026-09-14 | Post-brief close-out (natural hydro + 50 s, stirring 5 rev, spatial stats): [`11_THERMO_EM_COUPLING_AFTER_CHATGPT_TASK.md`](11_THERMO_EM_COUPLING_AFTER_CHATGPT_TASK.md). |
