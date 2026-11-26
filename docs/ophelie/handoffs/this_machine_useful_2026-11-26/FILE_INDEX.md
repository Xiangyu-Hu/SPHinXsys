# 文件说明

路径均相对本目录 `docs/ophelie/handoffs/this_machine_useful_2026-11-26/`。

## 文档

| 文件 | 用途 |
|---|---|
| `README.md` | 包入口与优先级 |
| `FILE_INDEX.md` | 本表 |
| `APPLY_ON_OTHER_MACHINE.md` | 另一台 Cursor 合入步骤 |
| `MANUAL_nc_zh.md` | 甲方自然对流详细说明书（几何、功率、初场、短烟） |

## 补丁 `patches/`

相对基线 `e9a280ff8`。可用 `git apply`；失败则改用 `overlay/` 同路径文件。

### `patches/nc/` 甲方自然对流（建议先合）

| 文件 | 落到仓库的路径 | 干什么 |
|---|---|---|
| `0001-client-geometry-H230-coil-frac.patch` | `tests/.../electromagnetic_ophelie_client_geometry.h` | 液深默认 230 mm；`adaptClientNcCoilStackToMeltHeight`；线圈丝带分数 |
| `0002-client-case-params-H230.patch` | `tests/.../electromagnetic_ophelie_client_case_params.h` | 参数默认液深 230；含搅拌辅助函数（若对方已有搅拌改动需手合） |
| `0003-frozen-q-client-nc-cli.patch` | `tests/.../test_3d_ophelie_french_natural_convection_frozen_q.cpp` | `--client-nc-melt-height`、`--client-nc-coil-z-frac-lo/hi`；吸收功率防 50 kW 误用 |
| `0004-glass-relax-adapt-coil.patch` | `tests/.../test_3d_ophelie_client_nc_glass_relax.cpp` | 改液深后重映射线圈 |
| `0005-doc248-C0-scripts-configs.patch` | `docs/.../doc248_client_cases/` 下 C0 脚本与 JSON | `02_run_C0_nc.sh`、`C0_effective`、探针几何、COMPILE/README |

### `patches/stir/` 搅拌（选合）

| 文件 | 落到仓库的路径 | 干什么 |
|---|---|---|
| `0001-client-stirring-geometry.patch` | `electromagnetic_ophelie_client_stirring_geometry.h` | 甲方搅拌几何调整 |
| `0002-french-stirring-geometry.patch` | `electromagnetic_ophelie_french_stirring_geometry.h` | 法国搅拌几何共用调整 |
| `0003-french-stirring-em.patch` | `test_3d_ophelie_french_stirring_em.cpp` | EM 搅拌求解器改动（C5 等） |
| `0004-client-stirring-f1-heated-elastic.patch` | `test_3d_ophelie_client_stirring_f1_heated_elastic.cpp` | 弹性桨加热搅拌 |
| `0005-C2-script-and-effective.patch` | `03_run_C2_stir.sh`、`C2_effective.json` | C2 运行脚本与有效配置 |

## 完整快照 `overlay/`

与补丁同一批文件的**整文件**副本，目录结构与仓库根相对路径一致。  
补丁冲突时：打开 `overlay/<path>`，与对方当前文件做合并，再写回正式路径。

另含**新文件**（无补丁，直接拷贝即可）：

- `overlay/docs/.../04_run_C3_stir_stainless.sh`
- `overlay/docs/.../05_run_C4_stir_rigid.sh`
- `overlay/docs/.../06_run_C5_stir_em30.sh`
- `overlay/docs/.../configs/C3_effective.json`
- `overlay/docs/.../configs/C4_effective.json`
- `overlay/docs/.../configs/C5_effective.json`

## 初场与配置 `configs/`

| 路径 | 用途 |
|---|---|
| `gpt_ic/client_nc_gpt_H230_T_rz.csv` | **求解器初场**（`r_m,z_m,T_K`），C0 用这个 |
| `gpt_ic/client_nc_gpt_H230.meta.json` | 映射说明与轴线检查 |
| `gpt_ic/client_nc_gpt_H230_preview.png` | 初场预览 |
| `gpt_ic/OPHELIE_Client_NC_Initial_Temperature_Field-230mm.zip` | GPT 源包（图、公式头、参数） |
| `gpt_ic/source/client_nc_initial_temperature.h` | 可嵌入的温度公式（可选；当前走 CSV） |
| `gpt_ic/source/initial_temperature_parameters.json` | 节点与峰值诊断 |
| `gpt_ic/source/README.md` | GPT 包原说明 |
| `gpt_ic/source/OPHELIE_Client_NC_Initial_Temperature.png` | GPT 原图 |
| `ic_history/client_nc_hybrid_tongue_v6_*` | 上一版几何舌形初场备份（非首选） |

合入后建议拷到：`docs/ophelie/reference/doc248_client_cases/configs/`。

## 脚本 `scripts/`

| 路径 | 用途 |
|---|---|
| `apply_all.sh` | 一键 `git apply` 全部 nc+stir 补丁 |
| `revert_shared_sources.sh` | 把本包涉及的共享源码/脚本撤回到 HEAD |
| `build_hybrid_ic_and_fig6_q.py` | 历史 hybrid 初场生成器（可选；正式初场用 GPT CSV） |

运行脚本正文在 `overlay/docs/ophelie/reference/doc248_client_cases/`，不在此重复存放。

## 证据 `evidence/`（小，非完整 runs）

| 路径 | 用途 |
|---|---|
| `C0_nc_gptIC_10/summary.json` | 短烟通过与功率摘要 |
| `C0_nc_gptIC_10/french_T_rz_t0009.98.csv` | 最后一帧温度剖面（作图用；`valid=0` 为空箱） |
| `C0_nc_gptIC_10/T_rz_gptIC_vs_t10.png` | 初场 vs 10 s 对比图 |
| `C0_nc_gptIC_10/WHERE_FULL_RUN.txt` | 本机完整 runs/VTP 路径（不随包提交大结果） |

**不在包内**：所有 `doc248_client_cases/runs/**`、VTP、hybrid v1–v5、`local-deps/`。
