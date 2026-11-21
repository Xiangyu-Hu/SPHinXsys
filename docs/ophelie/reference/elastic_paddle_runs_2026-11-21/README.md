# Doc228 elastic paddle run pack (2026-11-21)

| Gate | Status | Notes |
|------|--------|-------|
| S0 solid reference | **use existing SPHINXSYS cases** | custom `doc228_s0_cantilever` abandoned; SYCL: `test_3d_passive_cantilever_LG` |
| T0 thermal | **use existing** | `test_3d_ophelie_thermal_1d_conduction` |
| E0 paddle hold | **PASS** 2026-11-21 | n=10119 n_clamp=649; stainless E; coarse dp WARN |
| client stir relax | **PASS** 2026-11-21 | Glass/Rotor/Wall Reload; n_rotor=4332 dp≈0.005 |
| F0 isothermal FSI | **PASS** 2026-11-21 | `ELASTIC_FSI_SHORT_PASS=1`; probe wired; e_off≈0.22 mm |
| F1 heated elastic | **PASS** 2026-11-21 | `HEATED_ELASTIC_INTEGRATION_PASS=1`; dT≈1.14 K; energy_rel≈1.4% |
| E1 paddle load | **PASS** 2026-11-21 | `LOAD_RESPONSE_PASS=1` `ELASTIC_PADDLE_ACTIVE=1`; e_off_ratio≈1.99 vs 2 |
| E2 paddle spin | **PASS** 2026-11-21 | `DRIVEN_ELASTIC_ROTATION_SHORT_PASS=1`; θ≈0.026 rad; full 2π not claimed |
| EM Q | later | real EM Q when current known |

Do not invent generic solid/fluid/heat harnesses for module-only gates.
