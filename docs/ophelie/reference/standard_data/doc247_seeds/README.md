# Doc247 production seeds (tracked)

| File | Use |
|------|-----|
| `SEED_NC_Fig8_hotmelt_101_ASSUMED_T_rz.csv` | French NC Fig8 hot-melt IC (101², T∈[1400,1550] K) |
| `SEED_NC_Fig6_Qshape_ASSUMED_rz.csv` | French NC Fig6 Q **shape** (ASSUMED); use with `--init-Q-rz-csv=` + `--prescribed-heat-power=` (EM off) |
| `SEED_STIR_Fig8a_axisym_ASSUMED_T_rz.csv` | French stirring Fig8a axisymmetric IC (81²) |

`*_dense201*` grids stay **local-only** (gitignore).  
Do **not** use `standard_data/fig8_T/fig8_T_rz_Tmax1500_ASSUMED.csv` (750–1425 K cold-shell mapping).

Fig6 Q path is **not** live EM: shape from digitization INITIAL_GUESS (inherits ~71% extrapolation), scaled to prescribed watts.
