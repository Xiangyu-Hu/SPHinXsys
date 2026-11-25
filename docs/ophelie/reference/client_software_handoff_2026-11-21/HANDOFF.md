# 甲方两案 → 软件机交接包（2026-11-21）

**软件机请先读：[`SOFTWARE_SYNC.md`](./SOFTWARE_SYNC.md)**（最短同步清单）。

分支：`feature/EM-Biot-Savart`（以 GitHub 最新提交为准）。  
用途：桌面 M1/M2 接入的**开发短跑预设**，不是合同定量终验。

## 门控一览

| 案 | 项 | 标签 | 证据 |
|----|----|------|------|
| NC | 审计 | exit 0 | `AUDIT_SUMMARY.md` + `effective_client_nc_sigma16.json` |
| NC | 等温静水 | passed=1 | `nc/client_nc_isothermal_hydro_summary.json` |
| NC | 受控 ΔT=20 | passed=1（A/B） | `nc/client_nc_dT20_*.json` |
| 搅拌 | E0 | `PADDLE_HOLD_PASS` | `stirring/E0_summary.json` |
| 搅拌 | E1 | `ELASTIC_PADDLE_ACTIVE` | `stirring/E1_summary.json` |
| 搅拌 | E2 短窗 | `DRIVEN_ELASTIC_ROTATION_SHORT_PASS` | `stirring/E2_summary.json`（**非整圈**） |
| 搅拌 | F0 | `ELASTIC_FSI_SHORT_PASS` | `stirring/F0_summary.json` |
| 搅拌 | F1 | `HEATED_ELASTIC_INTEGRATION_PASS` | `stirring/F1_summary.json`（规定体积热，**非 EM**） |

`effective_client_stir.json`：`mechanics_mode=elastic_driven`，`stirrer_rpm=10`，`generator_input_power_W=30000`，电流 unknown。

## 软件机约束

- 不要把短窗标成整圈 PASS；不要把规定 Q 标成已算 EM。
- 1000℃ / 1400℃ 定量窗口未齐，不可钉测点温度冒充 PASS。
- CAD：`client_nc_glass_relax/data/`、`client_stirring_glass_relax/data/`。

## 复跑命令

见同目录 `COMPILE_AND_RUN.sh` 与下方用户消息中的分步命令。完整物理说明：`docs/ophelie/archive/DOCS_201_250_MERGED.md`（Doc230）。
