# 到另一台电脑后如何查看 / 合入

本机已把法国搅拌 Phase5 工作隔离到独立分支与独立目录，避免继续污染共享路径。你可在另一台电脑上**先审阅，再按同样方式挪进正常代码树**。

## 1. 两台机器分别用什么分支

| 机器 | 建议分支 | 说明 |
|------|----------|------|
| 另一台（你正在改的） | `feature/EM-Biot-Savart` | 继续自己的修改 |
| 本机（4090）搅拌成果 | `feature/ophelie-fr-stir-phase5-4090` | 只含 handoff 文档树 + 真 Fig.10 种子/脚本 |

**不要**把本机后续搅拌提交再 push 到 `feature/EM-Biot-Savart`。

注意：共享分支 tip `a502e16d8` **已经包含** Fig.10 Q-correct 的 C++（`--fig10-q-correct`，**默认关闭**）。另一台 `git pull` 后会带上该能力，但不打开 flag 则行为与原先一致。若另一台也改了同一个 `test_3d_ophelie_french_stirring_em.cpp`，合并时可能只在该文件冲突。

## 2. 在另一台查看本机 handoff（推荐）

```bash
cd <你的 sphinxsys 仓库>
git fetch origin
git checkout feature/ophelie-fr-stir-phase5-4090
# 或不动当前分支，只读：
git show origin/feature/ophelie-fr-stir-phase5-4090:docs/ophelie/reference/doc247_longrun_four_cases/handoff_4090_fr_stir_phase5/MANUAL_CN.md | less
```

关键目录：

```text
docs/ophelie/reference/doc247_longrun_four_cases/handoff_4090_fr_stir_phase5/
  MANUAL_CN.md              # 说明书
  SUMMARY_muT_setup_CN.md   # 变粘报告摘要
  seeds/                    # 真 Fig.10 + Fig.8a
  scripts/                  # 跑批 / 数字化（种子已指向本 handoff）
  figures/                  # Fig.10 QA
```

本机大数据包（gitignore，需单独拷贝或 Drive）：

```text
docs/ophelie/reference/doc247_longrun_four_cases/overnight_6h/
  phase5_fr_stir_6cases_pack_20261126.tgz
  phase5_fr_stir_muT10s_pack_20261126.tgz
  portable_fr_stir_phase5_4090_FULL.tgz   # 若已生成：说明书+代码快照一体包
```

## 3. 审阅满意后，挪进「正常代码」的建议步骤

按本机同样口径，**不要两边同时改共享种子路径**：

1. 在另一台建/切到专用分支（或 merge `feature/ophelie-fr-stir-phase5-4090`）。
2. 需要真 Fig.10 时，从 `handoff_.../seeds/` 复制或让跑批脚本指向该路径（handoff 内脚本已默认指向 handoff seeds）。
3. C++ Q-correct：若共享分支已有 `a502`，一般无需再打补丁；若你的另一台分支**没有**该提交，可从本机专用分支 cherry-pick `a502e16d8`，或拷贝：
   - `electromagnetic_ophelie_doc247_q_shape_correct.h`（新文件）
   - `test_3d_ophelie_french_stirring_em.cpp` 中与 `--fig10-q-correct*` 相关的接线
4. 跑通后再考虑是否把 handoff 内容合并进标准 `doc247_seeds/` / `02_fr_stir/`（合并前与本机约定，避免双写冲突）。

## 4. 冲突规避一句话

- 另一台日常开发：停在 `feature/EM-Biot-Savart`，避开搅拌 handoff 路径即可。  
- 要本机搅拌成果：看 / merge `feature/ophelie-fr-stir-phase5-4090`。  
- 大结果包走文件拷贝或 Drive，不进 git。
