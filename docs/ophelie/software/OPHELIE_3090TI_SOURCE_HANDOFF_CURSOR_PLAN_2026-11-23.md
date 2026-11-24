# OPHELIE 3090 Ti 软件源码整理与上传任务书

日期 2026-11-23  
执行对象：目前连接 3090 Ti 软件机的 Cursor。  
本轮交付：一条可获取的 Linux 软件交接分支、一份可复建说明和可选的小型回放样本。接收端是用户身边运行 OPHELIE 物理工作的 RTX 3090。

## 0. 开工依据

3090 Ti 上一次交接记录的工作树是 `/home/ubuntu/worktree/SPHinXsys-Ophelie-software`，当时的基线为 `feature/EM-Biot-Savart` 上的 `832dbc782`，`software/` 仍有未跟踪内容。以上是历史定位线索，先核对实际 checkout、当前 HEAD、已提交和未提交内容、活跃任务及最新验收结果。

上次证据覆盖 Linux GPU 的自然对流等温短窗新算、搅拌粒子准备、Qt5/VTK 桌面编译和打包骨架。SSH 下的实际 GUI 交互、3090 Ti 本机搅拌 F0 新算和干净环境发行包验收当时待补。历史回放、编译和本机新算各自记账。若本次已经补做，用对应的命令、日志和结果更新结论。

参考：[总任务书](https://chatgpt.com/library?file_id=libfile_9d8f3991c3c481919d89cf1a0c593be6&account_id=personal) · [230 Linux 桌面与第三方审查任务](https://chatgpt.com/library?file_id=libfile_35edbed0b9c48191a04d6f5d81e6d8f5&account_id=personal) · [3090 Ti 历史精简包](https://chatgpt.com/library?file_id=libfile_4fff34e55c3881918d4f602b392097e7&account_id=personal)

## 1. T0｜查清实际源码和基线

在 3090 Ti 的远程终端进入软件工作树。以下为只读核查，先记录输出并复核实际路径：

```bash
git rev-parse --show-toplevel
git branch --show-current
git rev-parse HEAD
git status --short --untracked-files=all
git diff --stat
git diff --cached --stat
git ls-files --others --exclude-standard -- software
git worktree list
git log -8 --oneline --decorate
```

- 确认当前软件相对于物理集成分支的**共同祖先、已提交的软件改动、未暂存改动、已暂存改动和未跟踪源码**。若某个必须的源码被 ignore，先定位规则，再精确处理该文件。
- 核对 GUI、示例配置、适配层、安装脚本、CMake 与 native solver 之间的调用关系。每个实际运行入口记录目标名、完整命令、工作目录、需要的 STL/Reload、输出目录和所用提交。核查 Qt 主版本、VTK 的 Qt 集成模块、oneAPI/Codeplay/CUDA、生成器、驱动、可见设备和关键系统包的实际版本。记录需要重建的路径与运行时依赖来源。
- 原工作区若有正在运行的任务或混有物理 WIP，先在旁边建专门的交接 worktree，从确认的已提交基线起步，选择性复制经检查的未提交文件或应用相关补丁。保留原工作区及其暂存区。交接 worktree 完成后核对它覆盖了仍需保留的新文件。

## 2. T1｜整理提交内容

| 放入 Git 的内容 | 核查要求 |
|---|---|
| Linux 桌面源码和构建入口 | `software/desktop/` 的源码、资源、CMake 和实际启动所需配置。检查 Qt5/VTK 方案能在目标环境复建 |
| Linux CLI 适配层和示例 | `software/coldcrucible_solver/`、必要的 `software/bin/` 入口、两案的 `case.json`、已有 schema/测试。确认关键文件已跟踪 |
| 依赖与发行脚本 | `software/packaging/` 下的脚本、依赖检查说明、可复现的构建及启动步骤。将原机绝对路径改为项目相对路径、配置覆盖或明确的开发机参数 |
| 共享求解器或 CMake 补丁 | 仅提交 Linux 软件实际依赖且经核查的改动。按功能独立提交，列出依赖顺序、涉及目标和测试。重复实现及物理实验 WIP 保留在相应工作区 |
| 交接记录 | 建议在 `docs/ophelie/software/` 增加短说明，列出基线、提交序列、输入清单、重建命令、证据矩阵和已知缺口 |

版本化文件以可复建源码为主。构建树、CMake 缓存、可执行文件、系统库、完整发行包、Reload、批量 VTP 和大日志保留在原机或既有成果位置。甲方资产优先引用已获准的版本和相对路径，必要时按项目允许的渠道单独补传。记录所需 STL、物性文件、配置、校验和及 Reload 生成方式。可另给一份很小的自然对流和 F0 回放样本，含真实 VTP、帧时间和 `result-index.json`，并注明来源和计算平台。只有 summary 时标为摘要证据。

**暂存时用明确路径或逐块选择。** 提交前检查 `git diff --cached --name-status`、`git diff --cached --stat`、`git diff --cached --check` 和实际 diff。确认样例没有隐藏依赖原机绝对路径，已排除二进制、密钥、个人配置及未经许可的甲方原件。附一份未纳入交接但仍应保留的改动摘要，防止把已有物理结果误当作遗漏。

## 3. T2｜复建与最小验收

1. 以交接候选提交为基准，记录实际 CMake 配置、编译目标和 `info --json`。能在既有环境独立重跑时，验证两份配置及准备路径，做一次自然对流等温 GPU 短跑。查实际设备与物理时间推进，保留最小日志、`summary.json` 和 `result-index.json`。搅拌记录其实际达到的阶段，长窗或新算有现成证据时附上。
2. 给已有历史和新算标注原机器、版本、配置与物理时长。若运行依赖 WIP，记录对应补丁或提交。先确认可复建性，SSH 会话缺图形窗口的项目标为未测，并交给 3090 本地 Ubuntu 验收。
3. 记录 3090 继续工作所需的接续信息。包括 GPU 路径、Qt/VTK 选型、CLI 入口、预处理与 case 定位、进程取消方式、安装目录和运行目录区分、当前未完成项。对仍存在的机器绝对路径列出具体文件和字段。

## 4. T3｜推送并给接收端定位

在项目已批准的**私有远程**上新建交接分支。建议分支名 `ophelie-linux-3090ti-handoff`。若同名已存在，先检查它是否是同一交接历史，也可用带日期的新名字。由 Cursor 按实际远程名替换命令中的占位符，并从确认的交接分支运行：

```bash
git rev-parse HEAD
git status --short
git log --oneline --decorate -8
git push -u <private-remote> HEAD:refs/heads/<handoff-branch>
git ls-remote <private-remote> refs/heads/<handoff-branch>
```

核对远端返回 SHA 等于本地交接 HEAD。使用新交接分支，集成分支的更新由接收端验证后安排。若私有远程不可用，保留本地交接提交，使用获准的方式交付该分支的 Git bundle，并记录创建与校验结果。先确认接收端是否已有共同祖先。需要传输历史时，核对传输的历史范围。

## 5. 最终回传格式

- 已核对的私有远程名或 bundle 名、交接分支、**完整交接 SHA、共同基线 SHA、可合入提交列表及顺序**。
- 实际交接文件清单、尚未提交或单独传输的项目、甲方资产定位与校验和。
- 最新可复现的环境、构建、启动、测试命令，已通过、失败、未测记录和最小证据。
- 当前限制及给 3090 的第一步。把本任务书的配套文件 `OPHELIE_3090_LINUX_SOFTWARE_CONTINUATION_CURSOR_PLAN_2026-11-23.md` 与交接信息一起交给接收端。

**完成判据**：接收端能以交接 SHA 找到全部 Linux 软件源码与依赖说明，已提交和未提交改动的状态明确。远端可见性或 bundle 可用性已经验证。
