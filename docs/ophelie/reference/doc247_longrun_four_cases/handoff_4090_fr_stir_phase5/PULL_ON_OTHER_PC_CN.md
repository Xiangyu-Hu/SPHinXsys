# 另一台电脑：怎样拽下本机 push，且不动你正在用的同名分支

本机搅拌相关后续只在：

```text
origin/feature/ophelie-fr-stir-phase5-4090
```

你另一台日常开发仍用：

```text
feature/EM-Biot-Savart
```

**目标**：把专用分支拽下来看说明书 / 种子 / 脚本，**不要 merge、不要 checkout 覆盖你当前的 `feature/EM-Biot-Savart` 工作区**。

---

## 推荐做法：单独 worktree（最安全）

在另一台电脑的仓库里执行（先把未提交改动 commit 或 stash，避免手滑切分支）：

```bash
cd /path/to/your/sphinxsys   # 你另一台的仓库根目录

# 1) 只下载远程分支信息，不会改你当前工作区文件
git fetch origin feature/ophelie-fr-stir-phase5-4090

# 2) 在旁边目录开一个只读式工作树，专门看 handoff
git worktree add ../sphinxsys-fr-stir-handoff origin/feature/ophelie-fr-stir-phase5-4090
```

然后去旁边目录看：

```bash
cd ../sphinxsys-fr-stir-handoff
ls docs/ophelie/reference/doc247_longrun_four_cases/handoff_4090_fr_stir_phase5/
less docs/ophelie/reference/doc247_longrun_four_cases/handoff_4090_fr_stir_phase5/MANUAL_CN.md
```

此时：

- 原仓库仍停在你的 `feature/EM-Biot-Savart`，文件不动  
- handoff 在 `../sphinxsys-fr-stir-handoff`，互不影响  

看完可删掉 worktree：

```bash
cd /path/to/your/sphinxsys
git worktree remove ../sphinxsys-fr-stir-handoff
```

---

## 备选：不建目录，只“偷看”文件

```bash
cd /path/to/your/sphinxsys
git fetch origin feature/ophelie-fr-stir-phase5-4090

# 看说明书（不切换分支）
git show origin/feature/ophelie-fr-stir-phase5-4090:docs/ophelie/reference/doc247_longrun_four_cases/handoff_4090_fr_stir_phase5/MANUAL_CN.md | less

# 导出整个 handoff 到临时目录（仍不改当前分支）
mkdir -p /tmp/fr_stir_handoff_view
git archive origin/feature/ophelie-fr-stir-phase5-4090 \
  docs/ophelie/reference/doc247_longrun_four_cases/handoff_4090_fr_stir_phase5 \
  | tar -x -C /tmp/fr_stir_handoff_view
ls /tmp/fr_stir_handoff_view/docs/ophelie/reference/doc247_longrun_four_cases/handoff_4090_fr_stir_phase5/
```

---

## 明确不要做的事（会碰到你正在用的分支）

```bash
# 不要（除非你真的要合进去）
git checkout feature/ophelie-fr-stir-phase5-4090   # 会离开 EM-Biot-Savart
git merge origin/feature/ophelie-fr-stir-phase5-4090  # 会改你当前分支历史/文件
git pull origin feature/ophelie-fr-stir-phase5-4090   # 在错误分支上 pull 可能混进提交
```

你日常继续：

```bash
git checkout feature/EM-Biot-Savart   # 保持你自己的工作
# 正常 commit / push 你自己的改动即可
```

---

## 结果大包怎么拿（不在 GitHub）

`.tgz` **没有** push 到 GitHub（在本机 `overnight_6h/`，gitignore）。  
要结果/画图包请用本机拷贝或 Drive，例如：

```text
.../overnight_6h/portable_fr_stir_phase5_4090_FULL.tgz
.../overnight_6h/phase5_fr_stir_muT10s_pack_20261126.tgz
.../overnight_6h/phase5_fr_stir_6cases_pack_20261126.tgz
```

---

## 以后真要合进“正常代码”时

审阅满意后，再在**你选定的分支**上按 `APPLY_ON_OTHER_MACHINE_CN.md` / `MANUAL_CN.md` 合入；  
合入前先在另一台把 `feature/EM-Biot-Savart` 的工作 commit 干净，再单独开分支 merge `origin/feature/ophelie-fr-stir-phase5-4090`，不要直接在脏工作区硬合。
