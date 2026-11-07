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
| 2026-09-21 | Stage 1 thermal host+MMS (Doc 42): double Q/Robin identities; insulated GPU energy_rel ~1e-7; multi-face BC | Production area ~14% short; coupled H from float T·V still open |
| 2026-09-21 | Ledger H → double Σ ρ cp (T−T0) V; stirring 2 s validation_passed=1 | energy_residual_H_rel 1.02 → 1.3e-4; area gate still open |
| 2026-09-21 | `--boundary-area-normalize`: SPH face A → geometric A (h unchanged) | melt A ratios → 1.000; T0 loss matches 0-D |
| 2026-09-21 | Freeze normalize **on**; nonuniform-T flux host+reload gates | production rel_side/bot ~1e-4; Stage 1 R03 closed for uniform+linear T |
| 2026-09-21 | Stage 2 slice: analytic α=(1−i)/δ (R09); re-compare SPH/GetDP azi-avg | Doc 43 |
| 2026-09-21 | Stage 2 R05: coil-only rz dump vs GetDP azi-avg (P ratio **1.009**) | three-mode table in Doc 43 |
| 2026-09-21 | Stage 2 R06–R07: domain enlarge gate pass; common-support L2 | axial ΔP~1.6%; midplane corr Q≈0.90 |
| 2026-09-21 | Stage 2 R08: coil-only axis-B audit; remove 23%/5 mm hard-code | ~20% residual; enlarge +1.2%; js0 kept |
| 2026-09-21 | Stage 3 r0: 1-rev stirred T baseline (Doc 44) | validation_passed=1; E_rel≈1.2%; T_std=4.3 K |
| 2026-09-21 | Stage 3 sensitivity: EM interval + area-normalize off | interval robust; normalize Δ(loss)≈−8% |
| 2026-09-22 | Stage 3 5-rev obs + skull-δ=10 mm diagnostic; GPT pack | Doc 44 obs_5rev; δ lever ~−34% loss |
| 2026-09-22 | Doc 45 review: F01–F07; not a validated thermal baseline | Batch A/B started |
| 2026-09-22 | Doc 46 Batch A/B closeout: F01–F07 code+evidence; 2 s smoke + cadence twin | se20/se100 `validation_passed=1`; F05 `nonuniform_ok=1`; Batch C not started |
| 2026-09-22 | Doc 47: accept A/B with qualifications; cadence twin withdrawn; enter Batch C | F04 strict parse first; one-point coil-only B |
| 2026-09-22 | Doc 48: §8 one-point B — SPH≈filament; GetDP \|B\|~−20% + sign | artifacts under `em_rz_getdp/work/batch_c_*` |
| 2026-09-22 | Doc 49: GetDP B shape error (not scale); SimpleTable OK; no js0 retune | midplane ratio falls toward coil |
| 2026-09-23 | Doc 50: root cause — `coil_loops[]` + `Js=I/A`; axis \|B\| ~1% vs filament | retire legacy 2π js0 compensation |
| 2026-09-23 | Doc 51: eddy J/Q re-compare under Doc 50 GetDP; signed-B note | P_sph/P_getdp **0.694** (was 1.009); midplane corr_Q **0.966** |
| 2026-09-23 | Doc 52: `--csv-every=` + F02 genuine-cadence protocol | GPU pair pending (no driver) |
| 2026-09-23 | Doc 53: F03 exclusive-endpoint gate; F07 Mass/ρ/V; EM P-gap note | no blanket slack; timestamps still OPEN |
| 2026-09-23 | Doc 54: F03 nominal-tick freshness + F01 production-seq host | GPU smoke still pending |
| 2026-09-23 | Doc 55: GPU F02 genuine cadence + freshness/F07 on 3090 | `f02_cadence_ok=true`; CSV 30 vs 6 |
| 2026-09-23 | Doc 56: midplane \|J\| scale ≈√P; corrected 1-rev GPU | `validation_passed=1` @ 1 rev |
| 2026-09-24 | Doc 57: SPH Bcoil vs GetDP coil-only; P-gap≈(B_drive)² | midplane \|B\| ratio ≈0.84 |
| 2026-10-06 | Doc105/106 P0: wall dens compatibility | AcousticDensMode 0/1; P2b C_ref=0; A1 step+hold PASS |
| 2026-10-06 | Doc105/107 P1: closed-box pressure convention | base p_h only; strat a_net_inner~2e-4 vs buoyancy 0.098 |
| 2026-10-06 | Doc105/108 P2: French FS hydro audit | dens overwrite+zones; quiet hydro OPEN (no P2b copy) |
| 2026-10-06 | Doc105/109 P2b: FS atm density | p_total~0 on surface band; a_net still OPEN |
| 2026-10-06 | Doc105/110 P2c: wall zones + budgets | no_extra mech-edge; keep U_rms~7e-3; quiet OPEN |
| 2026-10-06 | Doc105/111 P2d: French AcousticDensMode | mode1+a_ref=−g no gain; default mode0 |
| 2026-10-06 | Doc105/112 P2e: ff vs wall force split | interior/FS residual = fluid–fluid; quiet OPEN |
| 2026-10-06 | Doc105/113 P2f: packing + ρ0 freeze | no_extra=∇p_h floor; dens_sum packing; quiet OPEN |
| 2026-10-06 | Doc105/114 P2g: deep vs near-FS interior | deep a_net=0.39 with freeze; no_extra deep=2.63 |
| 2026-10-06 | Doc105/115 P2h: settle-only dens + packing audit | pack 0.8=support; settle_dens overpack p_eos; freeze deep=0.37 |
| 2026-10-06 | Doc105/116 P2i: FreeSurface overpack cap | no_extra+cap deep a_net=0.22; full-field OPEN |
| 2026-10-07 | Doc105/117 P2j: cap quiet-keep + surface budget | keep deep/int quiet; a_surf~2.6 U_rms~5e-3 OPEN |
| 2026-10-07 | Doc105/118 P2k: surface ff residual on cap | ff/wall~3e4; fsatm worse; deep quiet held |
| 2026-10-07 | Doc105/119 P2l: surface incomplete + zonal U | a_surf_cmp=0.14 quiet; U floor=side uz~1e-2 |
| 2026-10-07 | Doc105/120 P2m: incomplete-vel-zero keep | U 5.2e-3→4.8e-3; side wall floor confirmed |
| 2026-10-07 | Doc105/121 P2n: incomplete force-cancel live | REJECTED; keep a_surf/U worse |
| 2026-10-07 | Doc105/122 P2o: side-band U floor | side0 U 5.2e-3→1.3e-3; complete-side dominates |
| 2026-10-07 | Doc105/123 P2p: residual U floor source | side0 keep8 U→1.07e-3; dt0.5/ivel0 weak |
| 2026-10-07 | Doc105/124 P2q: bottom/corner force + bot0 | wall/ff_bot~1.06; side0+bot0+ivel0 U=7.4e-4 |
| 2026-10-07 | Doc105/125 P2r: FS incomplete pack-scale | REJECTED; a_surf/U worse; pause ops design |
| 2026-10-15 | Doc145/146: CEP EM audit + parameter contract | Trad=303; σ conflict OPEN; no long runs |
| 2026-10-15 | Doc145/146: L0/L1 finite-section coil + offline weak-Q | A5/B7 I/Bz checks pass; F0/F1 GPU-blocked |
| 2026-10-15 | Doc145/147: E2 refine + E3 aind + F2/F3 | δ/R triad pass@dp025; F2/F3 short OK w/ caveats |
| 2026-10-15 | Doc145/148: L1 French EM wire + Q/mapping + next_round | L0/L1 ΔP~8%; NR scripts refuse long runs |
| 2026-10-15 | Doc145/149: Doc61 L1 audit + dense Q | audit ratio=1; matched Q map ok (nr20×nz12) |
| 2026-10-15 | Doc145/150: round summary + pack | tar.gz ~628KB; long runs still OFF |
| 2026-10-15 | Doc151 plan: L1 visual + mapping + short coupling | freeze L2; same-chain C1–C3 |
| 2026-10-15 | Doc151/152: P0/P1 geometry+visual | 8.4%→0.026% after turn-center unify |
| 2026-10-17 | Doc151/153: P2/P3/P4 mapping+C1 | Q shape L2~10.8% borderline; C1 passed 2s |
| 2026-10-17 | Doc151/154: C2/C3 sigma feedback | subprocess T→σ→L1 Picard; C2/C3 passed 2s |
| 2026-10-17 | Doc151/155: round summary + pack | tar.gz ~1.5MB; long runs still OFF |
| 2026-10-17 | Doc156 plan: spatial-σ consistent coupling | freeze L1; no volume_mean / end-flush as evidence |
| 2026-10-17 | Doc156/157–159: spatial σ + schedule + direct Q | D0–D2 / DIRECT passed; P7 partial |
| 2026-10-17 | Doc156/160: round summary + pack | tar.gz ~5.8MB; long runs still OFF |
| 2026-10-17 | Doc161/162: circulation viz + D2_DIRECT + mapping/thermal prep | early circulation signs; long runs still OFF |
| 2026-10-17 | Doc161/163: round summary + pack | tar.gz; long runs still OFF |
| 2026-10-18 | Doc177–179: CLL bridge + CK-match implicit viscosity + thermal ledger | classic_nb≠0; warm Umax~−8%; cold energy~0.013% |
| 2026-10-18 | Doc180: viscosity analytic refs (CK explicit + CPU implicit) + full EM–thermal–flow cold 2.3 s | P1/P2 PASS; u0+0.7/1.4/2.1; energy~6e-6; freeze thermal-dt=0.005; pack `ophelie_doc180_viscosity_coupling_2026-10-18.tar.gz` |
| 2026-10-18 | Doc183/184: field export + Python circulation diagnostics + A/B/C short matrix | flow yes, closed cell no at 2.3 s cold T0; pack `ophelie_doc183_python_circulation_2026-10-18.tar.gz` |
| 2026-10-19 | Doc186: bulk-drift classifier correction + uniform buoyancy U0/U1/U2 | classifier fixed; U0 quiet; U1/U2 persistent drift; downgrade circulation detectability |
| 2026-11-07 | Doc212: shared complex LegacyFlux EMF; cavity CFL; French finalize; F5/F6 | phase+edge PASS; field L2~5–6% PARTIAL; DivΣA/pairwise diagnostic-only |
| 2026-11-07 | Doc213 plan: LegacyFlux EMF skin/field closure (G0 baselines) | production stays shared EMF; no RHS swap |

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
