# Linux 软件交接说明（3090 Ti → 身边 3090）

日期：2026-11-24  
交接分支建议：`ophelie-linux-3090ti-handoff`  
共同基线：`feature/EM-Biot-Savart` @ `832dbc782`

## 1. 放在哪里

全部 Linux 软件源码与说明在仓库内：

```text
docs/ophelie/software/
  LINUX_SOFTWARE_HANDOFF_3090TI.md          # 本文件
  230_OPHELIE_LINUX_DESKTOP_...md
  ENVIRONMENT_AND_BASELINE.md
  OPHELIE_3090TI_SOURCE_HANDOFF_...md
  HEADLESS_STATUS.txt
  linux_src/                               # 可复建源码（无大二进制）
    coldcrucible_solver/                   # Python CLI 适配层
    bin/coldcrucible-solver
    bin/run_detached.sh
    desktop/                               # Qt5 + VTK GUI 源码 + CMake
    packaging/                             # 发行打包脚本与依赖说明
    examples/{natural-convection,mechanical-stirring}/
      case.json
      assets/*.STL  -> 相对软链到 tests/.../data/
```

**不上 Git：** `build/`、`dist/`、Reload、大批 VTP、SYCL/Qt 系统库。对面机器自行编译与 `prepare`。

## 2. 另一台机器第一步

```bash
git fetch origin
git checkout ophelie-linux-3090ti-handoff   # 或 merge 进你的物理开发分支（只增 docs/ 路径，冲突概率低）
# 可选：合并进 feature/EM-Biot-Savart 时用
# git merge origin/ophelie-linux-3090ti-handoff

# 配置（与物理机同款）后编译求解器目标，并：
cmake -G Ninja -B docs/ophelie/software/linux_src/desktop/build \
  -S docs/ophelie/software/linux_src/desktop -DCMAKE_BUILD_TYPE=Release
cmake --build docs/ophelie/software/linux_src/desktop/build -j"$(nproc)"

source /opt/intel/oneapi/setvars.sh --force   # 若需要
export ONEAPI_DEVICE_SELECTOR=cuda:gpu
./docs/ophelie/software/linux_src/bin/coldcrucible-solver info --json
./docs/ophelie/software/linux_src/bin/coldcrucible-solver validate \
  --config docs/ophelie/software/linux_src/examples/natural-convection/case.json
```

环境变量（可选，无绝对硬编码）：

- `OPHELIE_REPO_ROOT`：仓库根（一般自动探测）
- `OPHELIE_BUILD`：求解器 build 目录（默认 `$REPO/build`）
- `OPHELIE_LINUX_SRC`：本目录 `linux_src` 的绝对路径（GUI 示例按钮用）
- `OPHELIE_SOLVER`：CLI 入口路径
- `OPHELIE_EXAMPLES_ROOT`：示例根目录

## 3. 证据与缺口

| 项 | 状态 |
|----|------|
| NC CLI prepare+run（3090 Ti） | 曾 PASS（summary 未随源码入库；对面请重跑） |
| 搅拌 prepare | 曾 PASS；Reload 未入库，请重跑 prepare |
| 搅拌 F0 本机新算 | 断线未完成；可用交接包历史 F0 演示，勿夸大 |
| GUI | 已编译过；SSH 无 DISPLAY 未测 → 请在本地 Ubuntu 图形会话测 |
| 发行包 dist | 未入库；可用 `linux_src/packaging/make_linux_bundle.sh` 重生 |

## 4. 禁止夸大

- E2 短窗 ≠ 整圈；F1 规定热 ≠ 真实 EM；30 kW = 电源输入  
- 搅拌必须 `elastic_driven`，禁止刚体动画冒充  

## 5. 与物理开发的关系

本交接只新增 `docs/ophelie/software/**` 下文件，**不改求解器物理源码路径**。推送使用独立分支；合入 `feature/EM-Biot-Savart` 时通常只是文档/软件树增量，不覆盖物理 WIP。
