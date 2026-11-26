# 另一台电脑 / 另一个 Cursor：如何合入

目标：把本包内容写进**正式仓库路径**，并尽量少和本地已有改动冲突。

## 0. 先对齐分支

```bash
cd <SPHinXsys>
git status
git log -1 --oneline
# 理想情况：已包含或基于 e9a280ff8；若已超前，下面补丁可能需 3-way / 手合
```

## 1. 只合甲方自然对流（推荐先做）

```bash
PACK=docs/ophelie/handoffs/this_machine_useful_2026-11-26
git apply --check $PACK/patches/nc/*.patch && git apply $PACK/patches/nc/*.patch
```

某条失败时：

```bash
# 看失败补丁对应的 overlay 文件，与当前工作树同路径做合并
# 例如 geometry：
#   $PACK/overlay/tests/extra_source_and_tests/extra_src/shared/electromagnetic_ophelie/electromagnetic_ophelie_client_geometry.h
```

拷初场与 C0 配置（若补丁 0005 未完全应用）：

```bash
cp -a $PACK/configs/gpt_ic/client_nc_gpt_H230_* \
      docs/ophelie/reference/doc248_client_cases/configs/
# 可选源包
mkdir -p docs/ophelie/reference/doc248_client_cases/configs/gpt_nc_ic
cp -a $PACK/configs/gpt_ic/source/* \
      docs/ophelie/reference/doc248_client_cases/configs/gpt_nc_ic/
cp -a $PACK/configs/gpt_ic/OPHELIE_Client_NC_Initial_Temperature_Field-230mm.zip \
      docs/ophelie/  # 或 configs/ 下
```

重编并短烟（路径按对方机器改 `BUILD`）：

```bash
# 编译 frozen_q + client_nc_glass_relax
OUT=docs/ophelie/reference/doc248_client_cases/runs/C0_nc_gptIC_10 \
END_TIME=10 TARGET_ABSORBED=18000 H_SIDE=70 H_FREE=10 EMISSIVITY=0.4 \
INIT_T_RZ=docs/ophelie/reference/doc248_client_cases/configs/client_nc_gpt_H230_T_rz.csv \
Q_UPPER_KAPPA=3.5 COIL_Z_FRAC_LO=0.68 COIL_Z_FRAC_HI=0.97 \
bash docs/ophelie/reference/doc248_client_cases/02_run_C0_nc.sh
```

细节见 `MANUAL_nc_zh.md`。

## 2. 再合搅拌（仅当对方搅拌文件较干净）

```bash
git apply --check $PACK/patches/stir/*.patch && git apply $PACK/patches/stir/*.patch
cp -a $PACK/overlay/docs/ophelie/reference/doc248_client_cases/04_run_C3_stir_stainless.sh \
      $PACK/overlay/docs/ophelie/reference/doc248_client_cases/05_run_C4_stir_rigid.sh \
      $PACK/overlay/docs/ophelie/reference/doc248_client_cases/06_run_C5_stir_em30.sh \
      docs/ophelie/reference/doc248_client_cases/
cp -a $PACK/overlay/docs/ophelie/reference/doc248_client_cases/configs/C{3,4,5}_effective.json \
      docs/ophelie/reference/doc248_client_cases/configs/
```

若 `french_stirring_em.cpp` / `stirring_f1` 对方已大改：**不要强行 apply**，按 `overlay/` 里同名文件把「功率口径 / 几何 / CLI」相关片段手挪过去。

## 3. 合入后提交建议（对方机器）

- 正式代码改动与本 handoff 目录可以同 PR，也可以分两次。
- **不要**提交 `doc248_client_cases/runs/`、`*.vtp`、本机 `local-deps/`。
- `C0_effective.json` 的 `preferred_IC` 应指向 `client_nc_gpt_H230_T_rz.csv`。

## 4. 给 Cursor 的任务提示（可复制）

```
请阅读 docs/ophelie/handoffs/this_machine_useful_2026-11-26/ 下的
README.md、FILE_INDEX.md、APPLY_ON_OTHER_MACHINE.md、MANUAL_nc_zh.md。
优先把 patches/nc/ 合入正式源码与 doc248 C0 脚本；把 configs/gpt_ic/
拷到 doc248_client_cases/configs/。搅拌 patches/stir/ 仅在与本地搅拌改动
不冲突时合入，否则按 overlay 手合并。不要提交 runs/ 与 vtp。
```
