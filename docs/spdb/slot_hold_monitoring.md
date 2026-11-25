# 有槽保压长算：断线后如何检测

长算用 `nohup` 启动，不依赖 Cursor 会话。进程掉线后仍可用下面文件判断进度。

## 启动

```bash
# 默认：dp=1 mm，Δt=60 s，120 步 → 目标 7200 s
docs/spdb/scripts/run_slot_hold_background.sh

# 自定义
HOLD_STEPS=10 PHYSICAL_DT=60 docs/spdb/scripts/run_slot_hold_background.sh
```

默认输出目录：`docs/spdb/results/slot_hold_dt60_n120/`。

## 检测（任选）

```bash
# 一键摘要：是否存活、保压时钟、心跳年龄、日志尾部
docs/spdb/scripts/check_slot_hold_status.sh

# 或指定目录
docs/spdb/scripts/check_slot_hold_status.sh docs/spdb/results/slot_hold_dt60_n120

# 最小检查
cat docs/spdb/results/slot_hold_dt60_n120/hold_status.txt
cat docs/spdb/results/slot_hold_dt60_n120/hold_heartbeat
kill -0 $(cat docs/spdb/results/slot_hold_dt60_n120/job.pid) && echo alive
tail -f docs/spdb/results/slot_hold_dt60_n120/run.log
```

## 文件含义

| 文件 | 作用 |
|------|------|
| `job.pid` | 主进程 PID |
| `job_meta.txt` | 启动参数与目标保压时间 |
| `hold_status.txt` | 程序每 ~500 声学步 / 每接受一步覆盖写入：`hold_physical_s`、`phase`、`progress_frac`、`finished` |
| `hold_heartbeat` | 同上时刻的一行心跳；若进程仍在但心跳 >15 min，脚本会警告 |
| `run.log` | 完整 stdout/stderr |
| `spdb_slot_hold_*_metrics.txt` | **只有正常结束或中止收尾时才有** |

`phase` 常见值：`elastic`、`elastic_balance`、`post`、`accepted`、`finished`、`aborted`。

## 上次未完成长算

完整 7200 s 已于 2026-10-23 跑完：`results/slot_hold_dt60_n120/`，`finished=1`，`hold_is_7200s=1`。数字见 [`slot_hold_results.md`](slot_hold_results.md)。

此前一次在约 **3300 s** 被中断（无最终 metrics）的目录已挪到：

`docs/spdb/results/slot_hold_dt60_n120_incomplete_3300s_*`

不要把那次标成 7200 s 完成。