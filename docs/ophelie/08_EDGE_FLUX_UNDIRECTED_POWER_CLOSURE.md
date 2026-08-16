# Stage 1: Edge-flux undirected / volume-consistent power closure

> Status: **in progress** (2026-08-08)  
> Audience: project notes for advisor / ChatGPT discussion  
> Related: B1 in `archive/discussion/OPHELIE_GPT_DISCUSSION_BLOCKERS.md`, French plan Stage 1 (§5–§7), `03_EDGE_FLUX_PRODUCTION.md`

---

## 1. Why this work now

French natural glass EM (`test_3d_ophelie_french_natural_em`) already:

- solves complex edge-flux with `edge_res_red ≈ 9×10³`
- calibrates **50 kW** on `P_recon` with `passed=1`
- `I_peak/loop ≈ 169.7 A` at `σ=16 S/m`, `f=282 kHz`, `dp=0.01`, `n=10474`

But production still logs:

```text
P_graph_over_recon ≈ 1.05×10⁵
```

Same order as the 2026-06-08 blocker **B1** (`≈1.15×10⁵`). So this is **not a new bug**; it was deferred when calibration switched to `P_recon`.

New French full-coupling plan asks to close Stage 1 **before** `A_glass` / `σ(T)` / thermal coupling: undirected edge + edge-consistent Joule power with `P_edge ≈ P_J_recon`.

---

## 2. Two ChatGPT positions (must not confuse)

| When | Guidance | Practical effect |
|------|----------|------------------|
| 2026-06 (power-fix plan) | Graph `½ C e²` is **Laplace discrete energy**, not physical Joule; calibrate on **`P_recon`** | French H `passed=1`; B1 left open |
| 2026-08 (French full-coupling plan) | Build **undirected** edges, `P_ij=½ G\|ℰ\|²`, deposit half/half; require `P_edge` vs `P_J_recon` within ~5% | Re-opens B1 as Stage 1 blocker |

These are not contradictory if interpreted carefully:

1. June: *current* graph accumulator was unsafe for 50 kW.
2. August: a **properly defined** edge-consistent power *should* match reconstructed continuum power.

Question for advisor / ChatGPT: confirm that August `G` is allowed to be the **SPH pair weight already used in φ**, or must be a separate geometric conductance `σ A_eff/L`.

---

## 3. What we found in the present code

File: `electromagnetic_ophelie_edge_flux.h` → `ComputeOphelieEdgeFluxJouleHeatCK`

Legacy (pre-fix) semantics:

```text
power_i = Σ_j 0.25 * C_ij * e_ij²     # treated as Watts
Q_i     = power_i / V_i               # W/m³
P_graph = Σ_i power_i                 # Σ Watts?
```

But `C_ij` is the **same pairwise Laplace weight** as the φ operator:

```text
C_ij = σ̄_ij * (-2 r (dW_ij V_j) / (r² + ε h²))
```

so `L_i = Σ_j C_ij (φ_i − φ_j)` approximates a **pointwise** `∇·(σ∇φ)` residual, not a resistor-network Siemens conductance.

Dimensional / scaling argument (French natural, `V̄ ~ 3×10⁻⁶ m³`):

```text
P_graph / P_recon ~ 10⁵  ~  O(1 / V̄)
```

So the dominant error is **missing particle volume when promoting pair energy to Watts**, not merely “counting each SPH neighbor twice”.

Note: the factor `0.25` already anticipates **double visitation** (`i→j` and `j→i`) relative to undirected `½ G e²` with arithmetic-mean `G`. Undirected bookkeeping alone cannot explain a `10⁵` gap.

---

## 4. SPH-compatible closure (Phase A — implemented now)

Keep SYCL `Inner<>` neighbor loops (no external mesh / CSR). Do **not** change φ LHS/RHS in Phase A.

### 4.1 Volume-consistent graph power

```text
q_i     = Σ_j 0.25 * C_ij * e_ij²     # interpret as W/m³
Q_graph = q_i
P_i     = q_i * V_i                   # Watt share on particle i
P_edge  = Σ_i P_i
```

Optional symmetrized weight used only in the power kernel (solver unchanged):

```text
V̄_ij = 0.5 (V_i + V_j)
C^sym_ij = σ̄_ij * (-2 r (dW_ij V̄_ij) / (r² + ε h²))
```

### 4.2 Undirected host/device audit (diagnostic)

For neighbors with `j > i` only:

```text
α_ij = σ̄ * (-2 r dW / (r²+εh²))
P_undirected += 0.5 * α_ij * V_i * V_j * e_ij²
```

Expect `P_undirected ≈ P_edge` after Phase A (same weak-form pair energy, counted once).

### 4.3 Production policy (unchanged until gates pass)

- **Calibration / primary Q still `P_recon` / `JouleHeatEdgeRecon*`**
- Graph fields remain diagnostic until soft gate holds:

```text
0.5 ≤ P_edge / P_recon ≤ 2.0     (soft)
target after recon polish: relative difference < 5%, then < 1%
```

### 4.4 Phase B (not in this change)

- Force `G_ij = G_ji` inside **residual + RHS** (true single edge system `Bᵀ G B`)
- Optionally promote edge-consistent `Q` to production heating field
- Only then reconsider calibrating on `P_edge` instead of `P_recon`

---

## 5. Baseline numbers to quote in discussion

### Fixed-current probe (`I_peak/loop = 1 A`)

From earlier French natural run:

```text
P_recon ≈ 1.74 W
P_graph (legacy) ≈ 1.83×10⁵ W
ratio ≈ 1.05×10⁵
edge_res_red ≈ 9.4×10³
```

### Target-power 50 kW

```text
I_peak/loop ≈ 169.655 A
P_recon = 50000 W
P_graph (legacy) ≈ 5.264×10⁹ W
ratio ≈ 105280
Q_outer/Q_center ≈ 17.3
passed = 1
```

After Phase A, re-run the same commands and fill:

```text
P_edge_vol_consistent = ?
P_undirected = ?
P_edge / P_recon = ?
P_undirected / P_edge = ?
```

### Uniform-field MMS (2026-08-08)

**Run #1 (before complex-chain fix):** `complex_potential_real` failed with `P_graph=0`.

**Run #2 (after `|e_i|²+|e_r|²` graph + post-both-chain timing):** `summary passed=1`

```text
potential_field:              P_recon/P_exact=1  P_graph/P_exact=0.78285  P_undirected/P_graph≈1  P_legacy/P_recon=12232
induction_field:              same
complex_potential_real:       same as potential_field                                              passed=1
complex_induction_imag_chain: P_recon/P_exact≈1  P_graph/P_exact=0.78285  P_undirected/P_graph=1   passed=1
complex_induction_real_chain: same as imag_chain                                                   passed=1
```

Interpretation:

- Volume fix works: legacy ratio `12232` on this box (`dp=0.04`, `n=512`) matches old B1 `~1/V` inflation.
- Undirected audit matches volume-consistent graph (`ratio=1`) on all five cases.
- Complex imag/real chains now agree on graph power (potential on `phi_real` and induction on either A chain).
- Remaining `P_graph/P_exact≈0.78` is SPH discrete pair energy vs continuum `½σ|E|²` (soft diagnostic gap); **recon is exact** and stays the calibration quantity.
### French natural 50 kW (2026-08-08, after Phase A)

Reload `n=10474`, `dp=0.01`, `σ=16`, `f=282 kHz`, target-power:

```text
# unit-current probe
P_recon=1.73715 W
P_graph=0.617243 W
P_undirected=0.617246 W
P_graph_over_recon=0.355319
P_undirected_over_graph=1
P_legacy_unweighted_over_recon=105280   # old B1 still visible as diagnostic

# after I² calibration to 50 kW
I_peak/loop=169.655 A
P_recon=50000 W
P_graph=17766 W
P_undirected=17766 W
P_graph_over_recon=0.35532
edge_res_red=9418
phi_rel≈1.14e-4
passed=1
```

Interpretation:

- **B1 closed in the old sense**: ratio fell from `~1.05e5` to `0.355` (same volume fix).
- Undirected identity holds on the production geometry (`P_undirected/P_graph=1`).
- Soft gate `[0.5, 2]` for `P_graph/P_recon` is **not** met yet (`0.355`); larger discrete bias than uniform-box MMS (`0.78`), likely irregular/relaxed neighbors + non-uniform induction field.
- Production calibration unchanged and healthy: `P_recon=50000`, `I` same as pre-fix run.
- Stage 1 remaining: decide whether `0.355` is acceptable diagnostic bias, tighten recon/edge consistency toward `<5%`, and/or Phase B symmetrize residual `G`.

### French natural 50 kW (2026-08-11, Stage1.5 closeout: reload/dp + coil geometry + weak-current audit)

Reload `n=10474`, `dp=0.015`, `R=0.25 m`, `H=0.185 m`, `f=282 kHz`, `σ=16`, `coil_R=0.285 m`, `coil_stack_total_height=230 mm` (7 匝等间距，stack 与 glass center 对齐 via `z_min/z_max`):

```text
[ophelie][stage1.5] reload audit:
  declared_dp=0.015 dp_eff_from_meanVol=0.015 volume_rel_error=0.026839
  dp_rel_error=2.29e-08
[ophelie][stage1.5] neighbor-count: p50=30 mean=28.3359 min=12 max=39

unit-current probe (I=1 A/loop):
  P_graph_edge=1.58762
  P_total_recon=1.762
  P_graph_over_recon=0.901036
  P_undirected/P_graph=1
  edge_flux_weak_current_antisym (imag): q_antisym_rel_l2=5.23e-08
  edge_flux_q_spatial: Q_max/mean=2.47068 (soft_gate=1)

after I² calibration to 50 kW:
  P_graph_edge=45052
  P_total_recon=50000
  P_graph_over_recon=0.90104
  passed=1
```

Interpretation:

- **Stage1.5 gate pass**: reload/dp consistency ok (`volume_rel_error=2.68% < 3%`, `dp_rel_error~2e-8`).
- weak-current antisymmetry now extremely small (`q_antisym_rel_l2 ~ 5e-8`)，这满足了 Stage1.5 决策中对 weak-current audit 的预期。
- `P_graph_over_recon ≈ 0.901` 已在软门 `[0.5, 2]` 内，因此不再需要为 Stage2 阻塞而继续打磨 Stage1 的功率比。


Commands (user build tree):

```bash
# from build/.../test_3d_ophelie_french_natural_em/bin
./test_3d_ophelie_french_natural_em \
  --reload=1 --reload-dir="$RELAX_RELOAD" \
  --excitation-mode=target-power --target-power=50000

./test_3d_ophelie_edge_flux_power_uniform_field
```

---

## 6. Questions for advisor / ChatGPT

1. Is Phase A volume weighting the intended continuum identity for SPH pointwise Laplace `C_ij`, or should `G` be redesigned as geometric `σ A_eff/L` separate from φ weights?
2. After `P_edge/P_recon → O(1)`, should production Joule heat switch from LS recon to edge deposit, or keep recon as primary and edge as audit forever?
3. Must Phase B symmetrize residual/RHS before claiming Stage 1 closed, or is power-only closure enough to proceed to `A_glass`?
4. Complex phasor: confirm `P = Σ ¼ C (|e_r|²+|e_i|²)` (with volume weights) vs treating imag/real chains separately then summing.
5. Soft gate: keep calibrating on `P_recon` until `<5%`, or switch earlier?

---

## 7. Acceptance checklist (Stage 1 power)

- [x] Uniform-field test: `P_edge/P_exact ∈ [0.5,2]`, `P_recon/P_exact ∈ [0.5,2]`, identity `|J|²/(2σ)` vs `Re(J·E*)` still ~0
- [x] French natural 50 kW (2026-08-11, Stage1.5 closeout): `P_graph/P_recon ≈ 0.901` (was `0.355`); soft `[0.5,2]` green; undirected match OK
- [x] Antisymmetry / undirected count documented (`j>i` once)
- [x] Timeline + this note updated with post-fix numbers
- [ ] No `A_glass` / `σ(T)` merge until soft gate green (per French plan) — **discuss with advisor**: allow proceed on `P_recon` with `0.355` diagnostic, or require `<5%` first

---

## 8. Change log

| Date | Change |
|------|--------|
| 2026-08-08 | Opened Stage 1 power-closure track; diagnosed B1 as volume semantics + deferred undirected audit; Phase A code path landed |
| 2026-08-08 | Uniform-field: `P_undirected/P_graph=1`, `P_graph/P_exact≈0.78`, legacy ratio `~1e4`; fixed complex graph to use both phasor chains |
| 2026-08-08 | Uniform-field **summary passed=1** (5/5); ready for French natural 50 kW re-check of `P_graph_over_recon` |
| 2026-08-08 | French natural 50 kW: `P_graph_over_recon` **105280 → 0.355**; `P_undirected/P_graph=1`; `passed=1`; soft `[0.5,2]` still open |

---

## 9. Copy-paste brief for advisor / ChatGPT

```text
Context: SPHinXsys SYCL OPHELIE-like complex edge-flux on French natural glass
(R=0.25 m, H=0.185 m, σ=16, f=282 kHz). Production calibrates 50 kW on P_recon;
legacy P_graph/P_recon ≈ 1.05e5 (blocker B1 from 2026-06).

Diagnosis: C_ij is SPH pairwise Laplace weight (pointwise ∇·σ∇), not network Siemens.
Legacy code treated Σ_j 0.25 C e² as Watts and set Q=that/V → ~1/V inflation.
Factor 0.25 already accounts for double neighbor visitation; undirected recount alone
cannot explain 1e5.

Phase A (implemented, needs rebuild+rerun):
  Q_graph_i = Σ 0.25 C_ij e²          (W/m³)
  P_graph   = Σ Q_graph_i V_i
  P_undirected = Σ_{j>i} 0.5 α V_i V_j e²
Calibration still on P_recon. φ LHS/RHS unchanged.

Please confirm:
1) Is volume weighting the right continuum identity for this C_ij?
2) After P_graph/P_recon → O(1), switch production Q to edge deposit or keep recon?
3) Must we symmetrize residual/RHS (Phase B) before A_glass?
```

