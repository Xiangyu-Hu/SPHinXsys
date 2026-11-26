# Doc248 client STIR IC seeds (ASSUMED)

Remapped from French stirring digitization onto client melt geometry
(`R=0.1415 m`, `H=0.230 m`, `z_bottom=0.001 m`).

| File | Use |
|------|-----|
| `SEED_CLIENT_STIR_Fig8a_T_rz.csv` | `--init-T-rz-csv=` (absolute Kelvin from Fig.8a) |
| `SEED_CLIENT_STIR_Fig10_Q_rz.csv` | `--init-Q-rz-csv=` (Q_norm; runtime scales to `--heat-power`) |
| `remap_french_stir_seeds_to_client.py` | regenerate from `docs/ophelie/搅拌caase起始/` |
| `figures/CLIENT_STIR_Fig8a_Fig10_remap_QA.png` | source vs remap QA |

**Not** literature PASS. Radial origin = crucible axis (not eccentric paddle).
