# CPU Demo 收尾跑法（主缺口）

对照任务书：主缺口是 **F3 1800 s**（950℃ 保温段）和 **F2 长时对照**。GPU 暂搁置。

已归档的 F3 1200 结果：`docs/spdb/results/cpu_program_temp_F3_end1200/`。

## 启动前清理（避免 output_backup 秒挂）

库在存在 `./output` 时会 `rename → output_backup`；若 `output_backup` 半残会抛异常。每次开跑前在 `bin/` 执行：

```bash
BIN=/home/ubuntu/software/SPHinXsys-SPDB/build-cpu-icpx/tests/extra_source_and_tests/3d_examples/test_3d_spdb_programmed_temperature/bin
cd "$BIN"
# 可选：先把当前 output 挪走备份
if [ -d output ]; then
  mv output "output_saved_$(date +%Y%m%d_%H%M%S)"
fi
rm -rf output_backup
mkdir -p output
```

## 顺序（串行，同一 bin，勿并行）

### 1) 主 demo：F3 1800（优先）

```bash
cd /home/ubuntu/software/SPHinXsys-SPDB
source /opt/intel/oneapi/setvars.sh
BIN=build-cpu-icpx/tests/extra_source_and_tests/3d_examples/test_3d_spdb_programmed_temperature/bin
cd "$BIN"
# 先做上面的清理

tmux new -s spdb_cpu_f3_1800
./test_3d_spdb_programmed_temperature --F3 --dp 0.001 --end_time 1800 --macro_dt 10 \
  2>&1 | tee run_F3_1800.log
# Ctrl+b d 脱离；回来：tmux attach -t spdb_cpu_f3_1800
```

成功标志：`hold_status.txt` 中 `finished=1`、`physical_s=1800`、`aborted=0`、`temperature_C=950`。
产物：`spdb_program_temp_F3_thermal_arrhenius_dp0.001_end1800_metrics.txt`、`program_temperature_observables.csv`。

跑完后立刻归档：

```bash
ARCH=/home/ubuntu/software/SPHinXsys-SPDB/docs/spdb/results/cpu_program_temp_F3_end1800
mkdir -p "$ARCH/vtk"
cp -a spdb_program_temp_F3_thermal_arrhenius_dp0.001_end1800_metrics.txt hold_status.txt \
  program_temperature_observables.csv run_F3_1800.log "$ARCH/"
cp -a output/*.vtp "$ARCH/vtk/" 2>/dev/null || true
```

### 2) F2 1200 对照（F3−F2 非弹性差）

```bash
# 仍在同一 BIN；先清理 output / output_backup（同上）
tmux new -s spdb_cpu_f2_1200
./test_3d_spdb_programmed_temperature --F2 --dp 0.001 --end_time 1200 --macro_dt 10 \
  2>&1 | tee run_F2_1200.log
```

归档到 `docs/spdb/results/cpu_program_temp_F2_end1200/`。

### 3)（次要）敏感性 — 主跑通后再做

```bash
# 门函数硬开关
./test_3d_spdb_programmed_temperature --F3 --hard_gate_900 --dp 0.001 --end_time 1200 --macro_dt 10 \
  2>&1 | tee run_F3_1200_hardgate.log

# 宏步减半
./test_3d_spdb_programmed_temperature --F3 --dp 0.001 --end_time 1200 --macro_dt 5 --max_macro_dt 15 \
  2>&1 | tee run_F3_1200_macro5.log
```

## 监控

```bash
tail -f hold_status.txt
# 或
tail -f run_F3_1800.log
```

## 墙钟粗估（dp=1 mm，CPU）

| 工况 | 参考 |
|------|------|
| F3 1200 | ~4.7 h（已完成） |
| F3 1800 | 可能 ~6–8 h+（多 600 s 保温宏步） |
| F2 1200 | 与 F3 同量级、略少（无 Arrhenius） |
