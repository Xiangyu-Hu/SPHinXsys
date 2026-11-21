# Client case audit (Doc229)

- repo_root: `/home/ubuntu/SPHinXsys`
- output: `/home/ubuntu/SPHinXsys/docs/ophelie/reference/client_case_audit_2026-11-20`

- formula_check: PASS (`client_formula_check.csv`)
  - 1473 K: mu=5.874 (ref 5.87400286), sigma=0.280681 (ref 0.28068105)
  - **retired**: pow(10, 4.5-8500/T) client sigma reading
- geometry NC/stir manifests: PASS/WARN see CSV
- NC nominal probe solver_m = [4.5166e-06, -7.99561e-06, 0.13] status=unconfirmed_sensor_location
- effective presets written: NC sigma16, NC client_ln_sigma, stirring
  - frozen mu @ 1273.15 K = 14.532 Pa·s
  - sigma16 @ ref = 16 ; client_ln @ ref = 0.113454
  - generator_input_power_W=30000 (not absorbed)
  - coil_current: **unknown** (no invented kA)

## Explicitly NOT claimed

- NC thermocouple 1000°C quantitative PASS
- EM absorbed power equal to 30 kW input
- Doc228 elastic paddle E0–F1
- same-current sigma16 vs client_ln EM Q compare (needs confirmed current + EM run)
