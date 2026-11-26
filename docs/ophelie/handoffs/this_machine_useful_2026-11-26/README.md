# 本机有用工作归档（2026-11-26）

给**另一台电脑上的 Cursor**：本目录是这台机器上值得保留的改动集合。  
**正式仓库路径里不要直接改的共享文件**，已撤回 `HEAD`；这里用补丁 + 完整文件快照保存。

| 优先级 | 内容 |
|---|---|
| 必合 | 甲方自然对流 C0：H=230、线圈丝带 CLI、GPT 初场、C0 脚本 |
| 选合 | 搅拌 C2–C5 脚本与搅拌源码改动（另一台若已改搅拌，请人工合并） |
| 不合 | `runs/`、VTP、hybrid v1–v5、`local-deps/` |

下一步请读：

1. [`FILE_INDEX.md`](FILE_INDEX.md) — 每个文件干什么  
2. [`APPLY_ON_OTHER_MACHINE.md`](APPLY_ON_OTHER_MACHINE.md) — 如何挪进正式代码  

基线提交：`e9a280ff8`（`feature/EM-Biot-Savart`）。若对方已超前，优先用 `overlay/` 做三路合并。
