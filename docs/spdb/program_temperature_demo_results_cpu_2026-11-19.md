# 程序温度双槽 Demo — CPU 结果（供先看）

来源：`build-cpu-icpx/.../test_3d_spdb_programmed_temperature/bin/`
材料：`demo_assumption`（E=120 GPa，α=12e-6）。**不是** 20Cr13 全温实测预测。
GPU 通路暂搁置；以本文件为当前可交付结果。

## 已完成矩阵

| 工况 | physical_s | T 终 | gate | aborted | 孔净高终值 | 保压段孔 Δ | 穿透终值 | mean 非弹性 | wall |
|------|------------|------|------|---------|------------|------------|----------|-------------|------|
| geometry | — | — | — | 0 | 初值 ~21 mm | — | — | — | 快 |
| F1 | 0（预平衡） | 20℃ | 0 | 0 | ~21.00 mm | 弹性闭合约 −1.91 µm | ~0.340 µm | 0 | 短 |
| F0 | 60 | ~76℃ | 0 | 0 | 指标失效记 0 | 勿当孔塌 | 0 | 0 | ~2951 s |
| F2 | 60 | ~76℃ | 0 | 0 | 0.0209981 | ~−6 nm | 0.340 µm | 0 | ~827 s |
| F3 | 60 | ~76℃ | 0 | 0 | 0.0209981 | ~−8 nm | 0.340 µm | 0 | ~813 s |
| F2 | 300 | 300℃ | 0 | 0 | 0.0209981 | ~−16 nm | 0.340 µm | 0 | ~4375 s |
| F3 | 300 | 300℃ | 0 | 0 | 0.0209981 | ~−18 nm | 0.340 µm | 0 | ~4366 s |
| **F3** | **1200** | **950℃** | **1** | **0** | **0.0209980** | **~−48 nm** | **0.340 µm** | **~2.7e-9** | **~16955 s (~4.7 h)** |

`F3 1800`：曾因 `output_backup` 秒挂，**尚未成功跑完**。收尾跑法见 [`program_temperature_demo_cpu_finish.md`](program_temperature_demo_cpu_finish.md)。
F3 1200 已归档：`docs/spdb/results/cpu_program_temp_F3_end1200/`。

## F3 1200（主看点）

- 完整跑通：`physical_s=1200`，`aborted=0`，`finished=1`
- 程序升温到 **950℃**；约 **1080 s** 过 850℃ 门，`gate` 升至 1
- 门后出现极小非弹性应变（mean ~1e-9 量级）——`demo_assumption` 下量级仅作通路验证，**不能**当工艺塌陷预测
- 保压段孔变化远小于弹性闭合约 1.9 µm；穿透 ~0.34 µm，须与孔变分开报
- 日志：`run_F3_1200.log`；metrics：`spdb_program_temp_F3_thermal_arrhenius_dp0.001_end1200_metrics.txt`

## 如何本地打开

```bash
BIN=/home/ubuntu/software/SPHinXsys-SPDB/build-cpu-icpx/tests/extra_source_and_tests/3d_examples/test_3d_spdb_programmed_temperature/bin
cd "$BIN"
cat hold_status.txt
cat spdb_program_temp_F3_thermal_arrhenius_dp0.001_end1200_metrics.txt
# VTK（若仍保留）
ls output/*.vtp | tail
```

## 若还要跑满 1800 s（950℃ 保温段）

单独、串行（勿与其它算例并行抢同一 `bin/`）：

```bash
cd /home/ubuntu/software/SPHinXsys-SPDB
source /opt/intel/oneapi/setvars.sh
BIN=build-cpu-icpx/tests/extra_source_and_tests/3d_examples/test_3d_spdb_programmed_temperature/bin
cd "$BIN"

tmux new -s spdb_cpu_f3_1800
./test_3d_spdb_programmed_temperature --F3 --dp 0.001 --end_time 1800 --macro_dt 10 \
  2>&1 | tee run_F3_1800.log
# Ctrl+b d 脱离
```

粗估墙钟：1200 s 约 4.7 h → 1800 s 可能再多数小时（宏步更多）。
