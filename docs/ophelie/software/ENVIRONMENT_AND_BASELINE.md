# ENVIRONMENT_AND_BASELINE（L0）

日期：2026-11-21（路径说明已于 2026-11-24 改为相对/可移植）  
工作树：仓库根（任意 checkout；勿写死本机绝对路径）  
对照：`230_OPHELIE_LINUX_DESKTOP_THIRD_PARTY_REVIEW_CURSOR_PLAN_2026-11-21.md` §2

## 1. 环境（已查）

| 项 | 值 |
|----|-----|
| OS | Ubuntu 24.04.4 LTS (noble) |
| CPU | 16 threads |
| GPU | NVIDIA GeForce RTX 3090 Ti, 24564 MiB |
| Driver | 580.178.04 |
| SYCL | `sycl-ls` → `[cuda:gpu] RTX 3090 Ti 8.6 [CUDA 13.0]` |
| `ONEAPI_DEVICE_SELECTOR` | `cuda:gpu` |
| oneAPI | `/opt/intel/oneapi/setvars.sh` 存在 |
| CMake | 3.28.3 |
| g++ | 13.3.0 |
| Python | 3.12.3 |
| Qt6 | **未安装**（无 qmake6 / Qt6Widgets pkg） |
| VTK | **未安装** |
| 图形会话 | `DISPLAY`/`WAYLAND` **unset**（当前像纯 SSH；桌面验收需图形会话） |
| 构建目录 | **尚无** `build/` |

## 2. Git 基线

| 项 | 值 |
|----|-----|
| 分支 | `feature/EM-Biot-Savart` |
| SHA | `832dbc782` |
| 跟踪 | `origin/feature/EM-Biot-Savart` |
| 本地未跟踪 | `docs/ophelie/software/230_OPHELIE_LINUX_DESKTOP_THIRD_PARTY_REVIEW_CURSOR_PLAN_2026-11-21.md`（及本文件写入后） |

审查包锁定此 SHA；后续软件提交另记，不自动追物理机 tip。

## 3. 求解器目标（源码已在树）

| 角色 | CMake target |
|------|----------------|
| NC relax | `test_3d_ophelie_client_nc_glass_relax` |
| NC 推进 | `test_3d_ophelie_french_natural_convection_frozen_q` + `--client-nc-*` |
| 搅拌 relax | `test_3d_ophelie_client_stirring_glass_relax` |
| F0 / F1 | `test_3d_ophelie_client_stirring_f0_elastic_fsi` / `..._f1_heated_elastic` |
| 门控 E0–E2 | `test_3d_ophelie_doc228_e0/e1/e2_*` |

几何：`tests/extra_source_and_tests/3d_examples/test_3d_ophelie_client_{nc,stirring}_glass_relax/data/`  
交接预设：`docs/ophelie/reference/client_software_handoff_2026-11-21/effective_client_{nc_sigma16,stir}.json`  
交接脚本（路径默认旧机）：`.../COMPILE_AND_RUN.sh` → 本机需改 `SRC`/`BUILD`。

## 4. 本机短跑

| 案 | 状态 | 证据 |
|----|------|------|
| NC relax + iso / ΔT20 | **本机复跑 passed=1**（2026-11-21 ~22:50） | `docs/.../handoff/.../nc/*.json`；reload 已生成 |
| 搅拌 F0 | 本机曾启动；满窗可后补。交接 summary `ELASTIC_FSI_SHORT_PASS=1` | 不锈钢亚步慢；软件演示可用已落盘结果 |

CMake（已验证）：Ninja + icx/icpx + SYCL `nvptx64-nvidia-cuda` + `CMAKE_PREFIX_PATH=$HOME/deps`（或本机第三方前缀）+ `-DFMT_HEADER_ONLY=1`；`SPHINXSYS_BUILD_3D_EXAMPLES=ON`。

## 5. 阻塞 / 缺项

1. 早期记录：Qt6/VTK 未装 → 后改用 **Qt5 + VTK9.1**（见 `linux_src/desktop/`）。  
2. **无 DISPLAY** → GUI 试跑需图形会话；CLI 可继续。  
3. L1 适配层现位于 `docs/ophelie/software/linux_src/`（info/validate/prepare/run）。

## 6. 下一步

1. 安装 Qt5 Widgets + VTK（与 Qt5 同 ABI），确认图形会话。  
2. L1：用 `docs/ophelie/software/linux_src/bin/coldcrucible-solver` 对两示例做 validate；需要时再 prepare/run。  
3. L2：NC 桌面闭环（见 `LINUX_SOFTWARE_HANDOFF_3090TI.md`）。
