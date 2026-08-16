# Master Development Timeline

> Consolidated from: `OPHELIE_PHI_DEVELOPMENT_LOG.md`, edge-flux run records, French three-way summaries, TEAM7 logs, RH200 status plans.  
> Originals: [`archive/plans/`](archive/plans/), [`archive/discussion/`](archive/discussion/)

---

## 1. Project goal

Build an **OPHELIE-inspired** complex phasor EM solver on SPHinXsys SYCL, then couple **Joule heating** to **SPH flow/thermal** for RH200 glass stirring. Midterm delivery is a **coupling prototype with auditable energy books**, not strict French cold-crucible quantitative reproduction.

---

## 2. Milestone timeline

| Date | Milestone | Outcome |
|------|-----------|---------|
| 2026-06 early | φ RHS lifecycle fix; MMS linear/discrete self-consistency | Sprint 1–2 passed |
| 2026-06 | One-sided Neumann boundary; slab/cylinder MMS | Sprint B first version |
| 2026-06-02 | Edge-flux Stage 1: pairwise SPH current equation | Sign MMS passed |
| 2026-06-08 | **Power calibration fix**: graph energy ≠ physical Joule power → use `P_recon` | French H production **passed @ 50 kW** |
| 2026-08-08 | French natural glass relax + EM (`R=0.25`, `H=0.185`) target-power 50 kW | `passed=1`, `I_peak≈170 A/loop`; legacy `P_graph/P_recon~1e5` still open |
| 2026-08-11 | Stage 3.0–3.4: A_glass Picard → Q→T; σ(T); diffusion; III.12 sensitivity | φ residual floor ~2.1e-4 accepted via dual gate |
| 2026-08-11 | Stage 3.4 closeout + Stage 3.5 plan: Natural σ_nat + 50 kW + A_glass | See `OPHELIE_STAGE3_4_CLOSEOUT_AND_FRENCH_NATURAL_NEXT_PLAN_2026-08-11.md` |
| 2026-08-11 | Stage 3.5 passed: 50 kW + A_glass (`P=50000`, `phi≈1.78e-4`) | Then Stage 4.0 Natural Robin/radiation BC |
| 2026-08-11 | Stage 4.0 Natural thermal BC smoke passed (loss channels >0; transient overcooling expected) | Stage 4.1 frozen-Q WCSPH+Boussinesq test added |
| 2026-08-08 | **Stage 1 Phase A**: volume-consistent graph power + undirected `j>i` audit | See `08_EDGE_FLUX_UNDIRECTED_POWER_CLOSURE.md`; calibration still `P_recon` |
| 2026-06-08 | Complex edge-flux: φ_r/φ_i, E_r/E_i, J_r/J_i | Production path |
| 2026-06 | French thermal one-way: EM → JouleHeat → explicit heating | Energy closure tests |
| 2026-06 | TEAM7 L2 Round 2–3: probe L2, high-σ scaling, operator audit | Side-track; not delivery gate |
| 2026-06 | TEAM7 P5: no-flux boundary normal / tangent-LS diagnostics | Partial closure |
| 2026-06 | RH200 Step 0–1: stirring + fake Joule → energy CSV | Baseline thermal |
| 2026-06 | RH200 Step 2: OPHELIE EM solve once on glass reload | P_scaled ≈ 50 kW |
| 2026-06 | RH200 Step 3b: **em-grid** fixed Euler `Q(x)` + per-step sampling | Fixes Lagrangian tag error |
| 2026-06 | Audit Task A+C: EXCITATION_TO_JOULE + HEATING_RATE CSV | Transparent scaling |
| 2026-06 | French-like material presets (Jacoutot @1473 K) | demo / stable / full μ=4 |
| 2026-06 | RH200 VTP: JReal/JImag on glass particles; em-grid live J mapping | ParaView animation |

---

## 3. Current production configuration

```text
Complex edge-flux (--ophelie-current-form=edge-flux, edge_flux_complex=true)
  → Biot–Savart A_coil (one-way, no A_ind Picard in RH200 midterm)
  → φ_r, φ_i (PCG/GMRES)
  → edge EMF drop → edge flux → E/J reconstruction
  → JouleHeatEdgeReconComplex
  → coil current calibration: I_new = I_old * sqrt(P_target / P_raw)
```

RH200 flow coupling:

```text
--joule-mode=em-grid (recommended)
  EM solve once → CIC deposit to Euler grid (full_oil.stl domain)
  each step: sample Q(x_i(t)), JReal/JImag(x_i(t)) on glass particles
  explicit Joule → T; glass/rotor diffusion; Simbody rotor FSI
```

---

## 4. Key fixes (do not regress)

1. **φ RHS overwrite bug** — Neumann/finalize must enter GMRES RHS pipeline.
2. **Edge-flux power scaling** — never calibrate on graph Laplace energy; use **`P_recon`** body integral.
3. **em-fixed vs em-grid** — frozen particle Q follows material tags (wrong); em-grid uses fixed lab-frame `Q(x)`.
4. **Energy book ρ₀** — Joule heating and thermal audit use reference density, not instantaneous ρ.
5. **Separate scale factors** — `em_power_scale_factor` (current²) vs `grid_rescale_factor` (CIC sampling).

---

## 5. Validation matrix (representative)

| Test | Validates |
|------|-----------|
| `test_3d_ophelie_edge_flux_sign` | Constant A → edge drop ≈ 0 |
| `test_3d_ophelie_edge_flux_power_uniform_field` | P_recon / P_exact = 1 |
| `test_3d_ophelie_edge_flux_scaling` | I×2 ⇒ P×4 |
| `test_3d_ophelie_phi_neumann_slab/cylinder` | Neumann boundary MMS |
| `test_3d_ophelie_phi_compatible_operator_mms` | G_c/D_c discrete consistency |
| `test_3d_ophelie_french_reduced --literature-mode` | H production @ 50 kW |
| `test_3d_ophelie_french_complex_joule_to_heat_one_way` | Thermal energy closure |
| `test_3d_ophelie_rh200_glass_em_stirring` | Full RH200 midterm case |

---

## 6. Open items (not blocking midterm demo)

- **A_ind** self-induction Picard — Stage 3.1 / 3.2b engineering pass (`phi_tol=2e-4`)
- **σ(T)** — thesis III.12 law wired (`--sigma-law=thesis-iii12`); Table-1 `16@1473` ≠ representative point
- Cold-crucible water cooling / free-surface losses
- **Stage 4** natural convection / Boussinesq / 50 kW full outer coupling
- TEAM7 strict L2 benchmark closure
- **Task B**: pointwise J·E formula residual audit
- French CAD / 400 kW supply-side model

Stage 3 progress snapshot (2026-08-11): [`archive/discussion/OPHELIE_STAGE3_PROGRESS_SNAPSHOT_2026-08-11.md`](archive/discussion/OPHELIE_STAGE3_PROGRESS_SNAPSHOT_2026-08-11.md)

---

## 7. Archive file manifest (repo-root plans moved here)

All former root-level plan files are in [`archive/plans/`](archive/plans/). Filename prefix indicates topic:

| Prefix | Count | Topic |
|--------|-------|-------|
| `OPHELIE_PHI_*` | 6 | φ operator, boundary, discretization |
| `OPHELIE_EDGE_FLUX_*` | 4 | Edge-flux solver and power fix |
| `OPHELIE_COMPLEX_EDGE_FLUX_*` | 2 | Complex mode + TEAM7 cross-plan |
| `OPHELIE_FRENCH_*` | 5 | French reduced geometry and midterm |
| `OPHELIE_TEAM7_*` | 8 | TEAM7 benchmark side-track |
| `RH200_*` | 3 | RH200 coupling and audit |
| `OPHELIE_CODE_REVIEW_*`, `OPHELIE_GUIDE_*`, etc. | 2+ | Reviews and audits |

Use `ls archive/plans/` for the full list.
