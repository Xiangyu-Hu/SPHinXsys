# 甲方自然对流（Doc248 C0）说明书

## 几何与功率

| 项 | 值 |
|---|---:|
| R | 0.1415 m |
| H | 0.230 m（液底 0.001，液面 0.231） |
| 发电机账本 | 30 kW |
| EM 吸收目标 | 18 kW（工程假设；可用 17–20） |
| 线圈 | 4 匝；默认丝带 15%–92%H；短烟用 68%–97%H |

轴线测温（物理深度锁定）：20 mm → 817.54 °C；70 mm → 1068.82 °C；120 mm → 980.66 °C。

## 代码要点

1. `client_geometry.h`：H=230；`adaptClientNcCoilStackToMeltHeight(lo,hi)`。
2. `frozen_q.cpp`：CLI 液深与线圈 frac；吸收功率覆盖法国 50 kW。
3. `glass_relax.cpp`：改液深后适配线圈。
4. `02_run_C0_nc.sh`：`TARGET_ABSORBED`、Robin、`INIT_T_RZ`、`Q_UPPER_KAPPA`、`COIL_Z_FRAC_*`。

Doc246 `Q_UPPER_KAPPA` 原二进制已有；短烟 κ=3.5。

## 初场

首选 `configs/client_nc_gpt_H230_T_rz.csv`（Fig.10 类 GPT 分层场）。  
GPT 按 230 mm 构造：0–150 mm 原剖面；150–230 mm 平滑接 810 °C 底。峰值约 1474 °C（r≈91 mm，深≈44 mm）。

## 短烟口径

```bash
OUT=docs/ophelie/reference/doc248_client_cases/runs/C0_nc_gptIC_10 \
END_TIME=10 TARGET_ABSORBED=18000 H_SIDE=70 H_FREE=10 EMISSIVITY=0.4 \
INIT_T_RZ=docs/ophelie/reference/doc248_client_cases/configs/client_nc_gpt_H230_T_rz.csv \
Q_UPPER_KAPPA=3.5 COIL_Z_FRAC_LO=0.68 COIL_Z_FRAC_HI=0.97 \
bash docs/ophelie/reference/doc248_client_cases/02_run_C0_nc.sh
```

结果摘要见 `evidence/C0_nc_gptIC_10/`。完整 VTP 仅留在本机 runs（不入库）。  
作图用最后一帧 CSV，丢掉 `valid=0`。
