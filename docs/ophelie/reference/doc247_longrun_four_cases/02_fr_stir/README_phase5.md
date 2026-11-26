# Doc247 法国搅拌 Phase5（合入说明）

法国桨是 **Simbody 刚体规定转动**（不是甲方 C2 的 NeoHookean 弹性桨）。

| 脚本 | 作用 |
|------|------|
| `run_fr_stir.sh` | 基线：Fig8a T + 真 Fig.10 Q-correct + 60 kW / 10 rpm / 常粘 |
| `run_fr_stir_muT_const10s.sh` | 常粘对照 μ≡4 |
| `run_fr_stir_muT_varT10s.sh` | 变粘 μ(T)/cp/k（CES Fig.2） |
| `run_fr_stir_height_relax_and_10s.sh` | 桨高 Δz=±20 mm relax + 10 s |
| `run_fr_stir_phase5_matrix_60s.sh` | rpm / 功率 / skull 矩阵 |

说明书：`MANUAL_phase5_FR_stir_CN.md`、`SUMMARY_muT_setup_CN.md`。  
种子：`standard_data/doc247_seeds/SEED_STIR_Fig{8a,10}_*`。
