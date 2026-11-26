# Doc247 production seeds (tracked)

| File | Use |
|------|-----|
| `SEED_NC_Fig8_hotmelt_101_ASSUMED_T_rz.csv` | French NC Fig8 hot-melt IC (101², T∈[1400,1550] K) |
| `SEED_NC_Fig6_Qshape_ASSUMED_rz.csv` | French NC Fig6 Q **shape** (ASSUMED); use with `--init-Q-rz-csv=` + `--prescribed-heat-power=` (EM off) |
| `SEED_STIR_Fig8a_axisym_ASSUMED_T_rz.csv` | French stirring Fig8a axisymmetric IC (81²) |
| `SEED_STIR_Fig10_Qshape_ASSUMED_rz.csv` | Stirring Fig.10 **Q shape** for `--fig10-q-correct` — **true visual digitization** (PLACEHOLDER_NOT_FIG10=false) |
| `SEED_STIR_Fig10_Qshape_PLACEHOLDER_from_Fig6_*.csv` | Old NC-Fig6 remapped PLACEHOLDER (kept for audit only) |

`*_dense201*` grids stay **local-only** (gitignore).  
Do **not** use `standard_data/fig8_T/fig8_T_rz_Tmax1500_ASSUMED.csv` (750–1425 K cold-shell mapping).

Fig6 Q path is **not** live EM: shape from digitization INITIAL_GUESS (inherits ~71% extrapolation), scaled to prescribed watts.  
Fig10: regenerate via `docs/ophelie/reference/standard_data/doc247_french_digitization/scripts/digitize_stir_fig10.py`.  
Fig10 is **Q shape only**, not a temperature field — pair with Fig8a T IC.
