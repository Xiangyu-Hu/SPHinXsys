# 执行进度

最后更新：2026-11-09

状态约定：`not_started` | `in_progress` | `blocked` | `deferred` | `done` | `na_other_branch`

## 1. 阶段总览

| 阶段 | 内容 | 状态 | 备注 |
|---|---|---|---|
| Phase 0 | 盘点与隔离 | done | 本机 2×4090 + 5995WX 已确认 |
| Phase 0b | 开跑前采集门禁（282 字段骨架+observer+missing 政策） | done | recorder + missing 政策 + A1/A2 副本 |
| Phase 1 | 公共构建与短测验收（非正式点） | done | A1 管线已验；A2 待用户编译短测 |
| Phase 2 | 二维 Poiseuille→梁→窗框→两相溃坝 正式精度 | done | A1–A4 正式矩阵完成；A5 deferred |
| Phase 3 | B 组规模 / group / CPU | in_progress | B2 主档已有；B3 启动；B1 几何后置 |
| Phase 4 | 本分支冻结包 | not_started | 不含工业 C |
| Phase 5 | 外卡 COMMON | not_started | 可选 |
| Phase 6 | A5 FSI2 移植 | done | L1–L3 share validate 三端齐 |
| Phase 7 | 工业 C1–C3 | in_progress | C1 GPU×3 完成（median 26477s）；CPU/外卡待续；见 `results/C1_gearbox/C1_DEVICE_HANDOFF.md` |

## 2. 案例进度

| ID | 状态 | 阻塞 / 依赖 | 下次动作 |
|---|---|---|---|
| A1 二维 Poiseuille | done | L1–L3 × GPU/CPU32/CPU64 validate 均 EXIT:0 | 详见 `results/A1_poiseuille/NOTES.md` |
| A2 二维梁振动 | done | L1–L3 × GPU/CPU32/CPU64 validate 均 EXIT:0 | 详见 `results/A2_oscillating_beam/NOTES.md` |
| A3 窗框传热 D4 | done | L1–L3 × GPU/CPU32/CPU64 validate 均 EXIT:0 | 详见 `results/A3_windows_frame_D4/NOTES.md` |
| A4 二维两相溃坝 | done | L3 CPU64 r01 segfault，正式用 r02 | 详见 `results/A4_two_phase_dambreak/NOTES.md` |
| A5 圆柱+柔性尾板 FSI2 | done | **L1–L3 × 三端 share validate 全 ok** | 本机 P0 收束；可进 B×3 |
| B1 | in_progress | GPU L1–L3 bench OK；L4 GPU=`ILLEGAL_ADDRESS`@inner_relation（同 uint32） | L4/L5 GPU blocked；可×3 或收束写笔记 |
| B2 | in_progress | L4/L5 CPU 不稳定 | GPU L1–L5 + CPU L1–L3 可用；L5 CPU 已杀；转 B3 |
| B3 | in_progress | **15/15 r01 矩阵齐**；可用 L1–L2 三端+L3-64；L4/L5 全 collapse | ×3 或查阻尼另议 |
| B4 | done | **L1–L5 × GPU/CPU32/CPU64 bench 全 ok** | ×3 或写笔记收束 |
| C1 | in_progress | T=**1.05 s**；CF关；**GPU FP32 ×3 完成** median **26477.2 s** | CPU×3 / 补 VTP；交接：`results/C1_gearbox/C1_DEVICE_HANDOFF.md` |
| C2 | na_other_branch | 工业分支 | — |
| C3 | na_other_branch | 工业分支 | — |

## 3. 本机硬件记录

| 项 | 值 |
|---|---|
| GPU | 2× NVIDIA GeForce RTX 4090，各 24564 MiB |
| 正式用卡 | **GPU 1**，UUID `GPU-6d9732cb-5aed-bfad-6d66-736689de9900`（`CUDA_VISIBLE_DEVICES=1`） |
| 勿用 | GPU 0（常被占用） |
| CPU | AMD Ryzen Threadripper PRO 5995WX，64 核 / 128 线程 |
| 主机角色 | 主开发 + 本机 CPU/GPU 配对；用户执行命令 |

## 4. 讨论结论日志

| 日期 | 结论 |
|---|---|
| 2026-11-03 | 主执行改到本机 4090+64 核；3090 不稳定降为可选 |
| 2026-11-03 | 工业 C 组不在本分支跑，另开分支 |
| 2026-11-03 | A5 等到其余内容跑完再做 GPU 移植 |
| 2026-11-03 | B1 几何后补；正确性 / 跨硬件 / 工业三块对论文同等重要 |
| 2026-11-03 | 产出 V3 测试计划：`01_TEST_PLAN_V3.md` |
| 2026-11-03 | 强制 GPU 1；助手不自跑；case 进 `extra_source_and_tests` 的 `gpu_paper` 副本；结果进 `results/`；tmux 断线可续看 → `03_RUN_CONVENTIONS.md` |
| 2026-11-03 | 明确采样频率政策：性能汇总 + A 组密探针；写入 `04_SAMPLING_POLICY.md` |
| 2026-11-03 | 收到并通读 V4 任务包；笔记见 `05_V4_REVIEW_NOTES.md`；采样以 V4 的 11/12/case 卡为准待确认 |
| 2026-11-03 | 用户确认：频率跟 V4；正式跑前采集/observer 一次就绪以免重跑；工具类统一 missing；FSI2 后置 → `06_PRE_RUN_INSTRUMENTATION_GATE.md` |
| 2026-11-03 | 采纳时长/物理输出 V2：物理时长以该文件为准；A1 正式终点改为 4 s；见 `08_DURATION_OUTPUTS_V2_ALIGNMENT.md` |
| 2026-11-04 | A1 L1–L3 × GPU FP32 / CPU FP32 / CPU FP64 正式 validate 全部完成；矩阵见 `results/A1_poiseuille/NOTES.md` |
| 2026-11-04 | 运行脚本改为 `tee`：用户终端/`tmux attach` 直接看时间步；助手不代跑仿真 |
| 2026-11-04 | 启动 A2：`test_2d_gpu_paper_a2_oscillating_beam`，T=2，PH/dp L1/L2/L3=10/20/40 |
| 2026-11-04 | A2 L1–L3 × GPU/CPU FP32/CPU FP64 正式 validate 全部完成；见 `results/A2_oscillating_beam/NOTES.md` |
| 2026-11-04 | 按顺序启动 A3：**D4**（D7 后置）；伪时间 0.02；`test_2d_gpu_paper_a3_windows_frame_D4` |
| 2026-11-04 | A3 D4 L1–L3 × GPU/CPU FP32/CPU FP64 正式 validate 全部完成 |

| 2026-11-04 | A4 L1–L3 正式矩阵完成（L3 CPU64 以 r02 为准） |
| 2026-11-04 | 启动 B2：`test_3d_gpu_paper_b2_dambreak`，T=10，LH/dp L1–L5=10/20/40/80/160 |

| 2026-11-05 | 中止 B2 L5 CPU；启动 B3 `test_3d_gpu_paper_b3_twisting_column` T=0.5 |

| 2026-11-05 | 启动 A5：`test_2d_gpu_paper_a5_fsi2`（经典 SPH；T=200；GPU 卡1） |

| 2026-11-06 | B1：接入 StanfordBunny_fixed.stl（sha256 a1b1badb…），创建 `test_3d_gpu_paper_b1_particle_relaxation`
| 2026-11-06 | B1：按 bunny AABB 重算域；去掉 Extrude，直接 STL 实体贴体；分记 levelset / particle_gen / relax 墙钟 |
| 2026-11-05 | A3/A5 正式迁 SYCL/CK：库补 Robin+LocalIsotropic InterParticleCoeff；A5 用 emitter/disposer |
| 2026-11-07 | B4 L1–L5 三端 bench 全齐；本机收束策略：先 A5 L2/L3，再 B×3/group，外卡仅 COMMON×3 |
| 2026-11-07 | A5 L1–L3 share 三端齐；B GPU 可用档 r02：B1×3+B2×5+B3×2+B4×5 = 15/15 ok |
| 2026-11-07 | B GPU 可用档 r03 亦 15/15 ok；除 B1-L3/B2-L3 外 spread 多 <4%；进 CPU FP32 ×3 |
| 2026-11-08 | CPU32 ×3：B1/B3/B4(L1–L4) 全齐；B2-L3 r03/r04 non_finite，以 r01/r02 报×2 |
| 2026-11-08 | 迁出 A3/A5 对 `src/` 的库补 → `extra_src/shared/gpu_paper/`；接线 B2+B4 `GPU_PAPER_GROUP` |

## 5. 本机（4090）收束清单 — 外卡工作量最小化

目标：凡能在本机跑完的都跑完；外卡只留 **B 组 COMMON×3**（每卡每例约 3 次），不在外卡做分辨率扫描 / A 组 / A5。

| 优先级 | 项 | 状态 | 说明 |
|---|---|---|---|
| P0 | **A5 L2/L3** 三端 validate（shared Reload） | **done** | L1–L3 全齐 |
| P1 | B 可用档 **×3**（r02/r03） | GPU 齐；CPU32 可用档齐（**B2-L3 仅×2**） | 可选：B4-L5 CPU×3 / group / COMMON 包 |
| P1 | B **L3 group** 扫描（B2+B4） | **24/24 ok** | G↑略加速；主表仍 G=64；见 `GROUP_SCAN_NOTES.md` |
| P2 | COMMON 冻结 + handoff 包 | 待 P0/P1 后 | 外卡只跑该点 |
| — | B1 L4+/B2 L4+CPU/B3 细档 collapse | **blocked 收束** | 不重跑；笔记已记 |
| — | C1–C3 | 另分支 | 不在本机本分支 |

## 6. 待你确认 / 提供的事项

- [ ] A3 最终用 D4 还是 D7（取决于本地 GPU 成熟度与 COMSOL 是否齐全）
- [x] B1 复杂几何文件（StanfordBunny_fixed.stl 已冻结）
- [ ] 外卡清单与优先级（2080S / 3090 / H20 哪些要配对）

更新本文件时：改对应行状态，并在「讨论结论日志」追加一行，避免口头约定丢失。
