# test_3d_ophelie_french_complex_joule_to_heat_one_way

French reload **end-to-end**: complex edge-flux EM → `JouleHeatEdgeReconComplex` → one-way thermal (no EM↔thermal feedback).

Default path is **coil-only**. Optional `--self-induction` enables Stage 3.1:

```text
A_coil + A_glass Picard → Q → T
```

Current French Natural geometry for this track:

```text
R=0.25 m, H=0.185 m, dp=0.015 m, n≈10474
coil R=0.285 m, 7 turns, z∈[-0.0225, 0.2075]
f=282 kHz, σ=16 S/m
```

## Run (cwd = `build/`)

**Always run GPU cases in a normal terminal (not Cursor sandbox).**

```bash
cmake --build . --target test_3d_ophelie_french_complex_joule_to_heat_one_way -j$(nproc)

RELAX_RELOAD=~/SPHinXsysSYCL/build/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_natural_glass_relax/bin/reload

# Stage 3.0 — coil-only EM → frozen Q → heating (energy closure gate)
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_complex_joule_to_heat_one_way/bin/test_3d_ophelie_french_complex_joule_to_heat_one_way \
  --reload-dir="$RELAX_RELOAD" \
  --dp=0.015 --glass-radius=0.25 --glass-height=0.185 \
  --frequency=282000 --sigma=16 \
  --coil-radius=0.285 --coil-num-loops=7 --coil-segments-per-loop=64 \
  --coil-z-min=-0.0225 --coil-z-max=0.2075 \
  --ophelie-current-form=edge-flux --ophelie-edge-flux-complex=1 \
  --thermal-dt=0.1 --thermal-steps=1

# Stage 3.1 — A_glass Picard then one-way thermal
# engineering phi gate: residual floor is ~1.78e-4 on this mesh/PCG setting
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_complex_joule_to_heat_one_way/bin/test_3d_ophelie_french_complex_joule_to_heat_one_way \
  --reload-dir="$RELAX_RELOAD" \
  --dp=0.015 --glass-radius=0.25 --glass-height=0.185 \
  --frequency=282000 --sigma=16 \
  --coil-radius=0.285 --coil-num-loops=7 --coil-segments-per-loop=64 \
  --coil-z-min=-0.0225 --coil-z-max=0.2075 \
  --ophelie-current-form=edge-flux --ophelie-edge-flux-complex=1 \
  --self-induction --self-induction-max-iter=25 --self-induction-relax=0.3 \
  --self-induction-tol=1e-3 --self-induction-phi-tol=2e-4 \
  --thermal-dt=0.1 --thermal-steps=1

# Stage 3.2 — provisional σ(T) outer loop (Arrhenius stand-in; not thesis III.12)
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_complex_joule_to_heat_one_way/bin/test_3d_ophelie_french_complex_joule_to_heat_one_way \
  --reload-dir="$RELAX_RELOAD" \
  --dp=0.015 --glass-radius=0.25 --glass-height=0.185 \
  --frequency=282000 --sigma=16 \
  --coil-radius=0.285 --coil-num-loops=7 --coil-segments-per-loop=64 \
  --coil-z-min=-0.0225 --coil-z-max=0.2075 \
  --ophelie-edge-flux-complex=1 \
  --self-induction --self-induction-max-iter=25 --self-induction-relax=0.3 \
  --self-induction-tol=1e-3 --self-induction-phi-tol=2e-4 \
  --sigma-t --sigma-t-max-iter=8 --sigma-t-relax=0.3 --sigma-t-tol=0.05 \
  --sigma-t-seed-delta-k=100 \
  --thermal-dt=0.1 --thermal-steps=1

# Stage 3.4 — thesis Fig. III.12 σ(T) (log10 σ = 3.7921 − 3179.8/T); expect σ(1473K)≈43, not Table-1 16
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_complex_joule_to_heat_one_way/bin/test_3d_ophelie_french_complex_joule_to_heat_one_way \
  --reload-dir="$RELAX_RELOAD" \
  --dp=0.015 --glass-radius=0.25 --glass-height=0.185 \
  --frequency=282000 --sigma=16 \
  --coil-radius=0.285 --coil-num-loops=7 --coil-segments-per-loop=64 \
  --coil-z-min=-0.0225 --coil-z-max=0.2075 \
  --ophelie-edge-flux-complex=1 \
  --sigma-t --sigma-law=thesis-iii12 \
  --sigma-t-max-iter=8 --sigma-t-relax=0.3 --sigma-t-tol=0.05 \
  --sigma-t-seed-delta-k=100 \
  --thermal-dt=0.1 --thermal-steps=1

# Stage 3.4b — thesis III.12 + A_glass Picard
# Dual-layer φ: preferred 2e-4 (Picard / preferred_ok); engineering hard 2.5e-4 (phi_ok)
# residual ~2.1e-4 → phi_gate_warn=1, passed=1 (not a physical blocker)
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_complex_joule_to_heat_one_way/bin/test_3d_ophelie_french_complex_joule_to_heat_one_way \
  --reload-dir="$RELAX_RELOAD" \
  --dp=0.015 --glass-radius=0.25 --glass-height=0.185 \
  --frequency=282000 --sigma=16 \
  --coil-radius=0.285 --coil-num-loops=7 --coil-segments-per-loop=64 \
  --coil-z-min=-0.0225 --coil-z-max=0.2075 \
  --ophelie-edge-flux-complex=1 \
  --self-induction --self-induction-max-iter=25 --self-induction-relax=0.3 \
  --self-induction-tol=1e-3 --self-induction-phi-tol=2e-4 \
  --sigma-t --sigma-law=thesis-iii12 \
  --sigma-t-max-iter=8 --sigma-t-relax=0.3 --sigma-t-tol=0.05 \
  --sigma-t-seed-delta-k=100 \
  --thermal-dt=0.1 --thermal-steps=1

# Optional Stage 3.4 diagnostic — tighten inner φ PCG 10× (check audit floor vs Krylov)
# ./... --phi-pcg-tol=5e-5 --phi-pcg-max-iter=20000 ...

# Stage 3.5 — 50 kW target-power + A_glass (calibrate I, re-Picard; |P−50kW|/50kW < 1e-3)
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_complex_joule_to_heat_one_way/bin/test_3d_ophelie_french_complex_joule_to_heat_one_way \
  --reload-dir="$RELAX_RELOAD" \
  --dp=0.015 --glass-radius=0.25 --glass-height=0.185 \
  --frequency=282000 --sigma=16 \
  --coil-radius=0.285 --coil-num-loops=7 --coil-segments-per-loop=64 \
  --coil-z-min=-0.0225 --coil-z-max=0.2075 \
  --ophelie-edge-flux-complex=1 \
  --self-induction --self-induction-max-iter=25 --self-induction-relax=0.3 \
  --self-induction-tol=1e-3 --self-induction-phi-tol=2e-4 \
  --target-power=50000 \
  --thermal-dt=0.1 --thermal-steps=1

# Production Natural σ(T) law (Table-1 anchor 16 @1473 + III.12 shape, clip [1e-6,30])
# --sigma-t --sigma-law=french-natural ...

# Stage 3.3a — A_glass + isotropic diffusion + cold-wall Dirichlet (no σ(T))
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_complex_joule_to_heat_one_way/bin/test_3d_ophelie_french_complex_joule_to_heat_one_way \
  --reload-dir="$RELAX_RELOAD" \
  --dp=0.015 --glass-radius=0.25 --glass-height=0.185 \
  --frequency=282000 --sigma=16 \
  --coil-radius=0.285 --coil-num-loops=7 --coil-segments-per-loop=64 \
  --coil-z-min=-0.0225 --coil-z-max=0.2075 \
  --ophelie-edge-flux-complex=1 \
  --self-induction --self-induction-max-iter=25 --self-induction-relax=0.3 \
  --self-induction-tol=1e-3 --self-induction-phi-tol=2e-4 \
  --use-literature-thermal=1 --thermal-diffusion=1 \
  --thermal-dt=0.1 --thermal-steps=30

# Stage 3.3b — σ(T) converge, then frozen-Q diffusion + cold wall (auditable energy/Dirichlet)
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_complex_joule_to_heat_one_way/bin/test_3d_ophelie_french_complex_joule_to_heat_one_way \
  --reload-dir="$RELAX_RELOAD" \
  --dp=0.015 --glass-radius=0.25 --glass-height=0.185 \
  --frequency=282000 --sigma=16 \
  --coil-radius=0.285 --coil-num-loops=7 --coil-segments-per-loop=64 \
  --coil-z-min=-0.0225 --coil-z-max=0.2075 \
  --ophelie-edge-flux-complex=1 \
  --self-induction --self-induction-max-iter=25 --self-induction-relax=0.3 \
  --self-induction-tol=1e-3 --self-induction-phi-tol=2e-4 \
  --sigma-t --sigma-t-max-iter=8 --sigma-t-relax=0.3 --sigma-t-tol=0.05 \
  --sigma-t-seed-delta-k=100 \
  --use-literature-thermal=1 --thermal-diffusion=1 \
  --thermal-dt=0.1 --thermal-steps=30

# Stage 4.1 — isotropic diffusion + cold-crucible Dirichlet shell (regression BC)
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_complex_joule_to_heat_one_way/bin/test_3d_ophelie_french_complex_joule_to_heat_one_way \
  --reload-dir="$RELAX_RELOAD" --dp=0.015 \
  --ophelie-edge-flux-complex=1 --thermal-diffusion=1 --thermal-bc=dirichlet-regression

# Stage 4.0 — French Natural production thermal BC (Robin side/bottom + free-surface conv/rad)
# h_side=300 selected within journal sensitivity range [100–400], not a unique prescribed value
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_complex_joule_to_heat_one_way/bin/test_3d_ophelie_french_complex_joule_to_heat_one_way \
  --reload-dir="$RELAX_RELOAD" \
  --dp=0.015 --glass-radius=0.25 --glass-height=0.185 \
  --frequency=282000 --sigma=16 \
  --coil-radius=0.285 --coil-num-loops=7 --coil-segments-per-loop=64 \
  --coil-z-min=-0.0225 --coil-z-max=0.2075 \
  --ophelie-edge-flux-complex=1 \
  --self-induction --self-induction-max-iter=25 --self-induction-relax=0.3 \
  --self-induction-tol=1e-3 --self-induction-phi-tol=2e-4 \
  --target-power=50000 \
  --thermal-bc=french-natural --use-literature-thermal=1 \
  --thermal-dt=0.1 --thermal-steps=5

# Stage 4.2 — Jacoutot Table 1 thermal props + Temperature VTP (ParaView)
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_complex_joule_to_heat_one_way/bin/test_3d_ophelie_french_complex_joule_to_heat_one_way \
  --reload-dir="$RELAX_RELOAD" --dp=0.015 \
  --ophelie-edge-flux-complex=1 --thermal-diffusion=1 \
  --use-literature-thermal=1 --thermal-steps=30 --thermal-state-recording=1 --thermal-record-interval=5
```

Prerequisite: `Reload.xml` from `test_3d_ophelie_french_natural_glass_relax` at matching `dp`.

## Thermal material presets

| Preset | CLI | ρ [kg/m³] | cp [J/(kg·K)] | k [W/(m·K)] | T₀ [K] |
|--------|-----|-----------|---------------|-------------|--------|
| **reduced** (default, regression) | `--thermal-material=reduced` | 2500 | 1200 | 1 | 300 |
| **literature** | `--use-literature-thermal=1` or `--thermal-material=literature` | 2750 | 1150 | 4 | 1473 |

Literature values: Jacoutot et al., *Chem. Eng. Process.* 47 (2008) Table 1 @ 1473 K.

## Acceptance

| Mode | Gates |
|------|--------|
| **default** (`thermal_diffusion=0`) | EM OK; vol-weighted ΔT / energy balance / `P·Δt` **< 5%**; mismatch vol **< 1%** |
| **`--self-induction`** | above + `J_rel < j_tol`; dual-layer φ: preferred=`phi-tol`, hard=`max(phi-tol,2.5e-4)`; `A_ind/A_coil > 0` |
| **`--target-power=`** | probe EM → calibrate `I` → re-solve; `|P−P_target|/P_target < 1e-3` |
| **`--sigma-t`** | above + σ spatial non-uniform + outer σ/P converged |
| **`--thermal-diffusion=1`** | EM OK; `max_T > T₀`; Dirichlet regression compliance **> 90%**; **no spurious energy creation** (`E_thermal ≤ E_joule`) |
| **`--thermal-bc=french-natural`** | diffusion + Robin side/bottom + free-surface conv/rad; side/bottom/free losses **> 0**; `E_thermal ≤ E_joule` |
| **`--sigma-t --thermal-diffusion=1`** | σ gates on coupling; then diffusion gates on **frozen-Q** post step (T reset to T₀ first) |

With diffusion, pointwise Joule closure vs frozen-Q formula is **not** required (heat redistributes / exits cold walls).

**Note on σ(T) sources:** CEP/CES **Figure 2** is the coupling flowchart, **not** σ(T). Traceable curves are Jacoutot thesis **Fig. III.12 / IV.10** (`--sigma-law=thesis-iii12|thesis-iv10`, sensitivity). **Production Natural:** `--sigma-law=french-natural` = Table-1 `σ(1473K)=16` × III.12 shape, clipped to `[1e-6, 30]`. See `docs/ophelie/reference/jacoutot_sigma_t/README.md` and Stage 3.4 plan.
## Related

- MMS diffusion + cold box: `test_3d_ophelie_thermal_diffusion_mms`
- Picard-only EM: `test_3d_ophelie_french_self_induction_picard`
- Implementation: `electromagnetic_ophelie_thermal_diffusion_one_way.h`, `electromagnetic_ophelie_joule_to_heat_one_way.h`, `electromagnetic_ophelie_aind_diagnostic.h`
