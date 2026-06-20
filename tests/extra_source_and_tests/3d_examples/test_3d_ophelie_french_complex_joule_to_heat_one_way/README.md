# test_3d_ophelie_french_complex_joule_to_heat_one_way

French reload **end-to-end**: complex edge-flux EM → `JouleHeatEdgeReconComplex` → one-way thermal (no EM feedback).

## Run (cwd = `build/`)

```bash
cmake --build . --target test_3d_ophelie_french_complex_joule_to_heat_one_way -j$(nproc)

# Stage 4.0 — frozen Q, explicit heating only (energy closure gate)
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_complex_joule_to_heat_one_way/bin/test_3d_ophelie_french_complex_joule_to_heat_one_way \
  --reload=1 --ophelie-current-form=edge-flux --ophelie-edge-flux-complex=1

# Stage 4.1 — isotropic diffusion + cold-crucible Dirichlet shell
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_complex_joule_to_heat_one_way/bin/test_3d_ophelie_french_complex_joule_to_heat_one_way \
  --reload=1 --ophelie-edge-flux-complex=1 --thermal-diffusion=1

# Stage 4.2 — Jacoutot Table 1 thermal props + Temperature VTP (ParaView)
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_complex_joule_to_heat_one_way/bin/test_3d_ophelie_french_complex_joule_to_heat_one_way \
  --reload=1 --ophelie-edge-flux-complex=1 --thermal-diffusion=1 \
  --use-literature-thermal=1 --thermal-steps=30 --thermal-state-recording=1 --thermal-record-interval=5

# optional overrides: --thermal-dt=1.0 --thermal-steps=3 --thermal-t0=300 --rho=2500 --cp=1200 --k=1
# material presets: --thermal-material=reduced|literature
```

Prerequisite: a folder containing `Reload.xml` with **~20500 GlassBody** particles (not the tiny 3-particle demo reload).

Generate once from `build/`:

```bash
./tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_reduced/bin/test_3d_ophelie_french_reduced \
  --glass-radius=0.325 --glass-height=0.50 --dp=0.02 \
  --relax=1 --state_recording=0 --relax-log-every=50
# writes build/.../bin/reload/Reload.xml
```

Use `--reload-dir=` only when sharing reload from another build target (e.g. `build/reload`).

**Note:** `<french_reload>` in docs is a placeholder — use a real path, e.g. `--reload-dir=$HOME/SPHinXsysSYCL/build/tests/extra_source_and_tests/3d_examples/test_3d_ophelie_french_reduced/bin/reload`.

VTP fields (`output/GlassBody_ite_*.vtp`): `Temperature`, `OphelieThermalDeltaT`, `JouleHeat`, `JouleHeatEdgeReconComplex` (complex edge-flux); with diffusion also `OphelieThermalLaplaceT`, `OphelieThermalConductivity`, `OphelieThermalBoundaryMask`.

## Thermal material presets

| Preset | CLI | ρ [kg/m³] | cp [J/(kg·K)] | k [W/(m·K)] | T₀ [K] |
|--------|-----|-----------|---------------|-------------|--------|
| **reduced** (default, regression) | `--thermal-material=reduced` | 2500 | 1200 | 1 | 300 |
| **literature** | `--use-literature-thermal=1` or `--thermal-material=literature` | 2750 | 1150 | 4 | 1473 |

Literature values: Jacoutot et al., *Chem. Eng. Process.* 47 (2008) Table 1 @ 1473 K. Regression gates use **reduced**; literature is for longer runs / ParaView comparison.

## Acceptance

| Mode | Gates |
|------|--------|
| **default** (`thermal_diffusion=0`) | EM OK; vol-weighted ΔT / energy balance / `P·Δt` **< 5%**; mismatch vol **< 1%** |
| **`--thermal-diffusion=1`** | EM OK; `max_T > T₀`; boundary Dirichlet compliance **> 90%**; **no spurious energy creation** (`E_thermal ≤ E_joule`) |

With diffusion, pointwise Joule closure vs frozen-Q formula is **not** required (heat redistributes / exits cold walls).

## Related

- MMS diffusion + cold box: `test_3d_ophelie_thermal_diffusion_mms`
- Implementation: `electromagnetic_ophelie_thermal_diffusion_one_way.h`, `electromagnetic_ophelie_joule_to_heat_one_way.h`, `electromagnetic_ophelie_french_thermal_material.h`, `electromagnetic_ophelie_thermal_vtp.h`
