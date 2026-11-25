# 程序温度双槽 Demo — 编译与运行（由你在本机执行）

工作区：`/home/ubuntu/software/SPHinXsys-SPDB`
算例：`tests/extra_source_and_tests/3d_examples/test_3d_spdb_programmed_temperature/`
后端：CPU，`build-cpu-icpx`，`SPHINXSYS_USE_SYCL=OFF`
材料：`demo_assumption`（E=120 GPa，α=12e-6）；**不是** 950℃ 实测预测。

**GPU/SYCL 版**（独立算例与构建目录）见 [`program_temperature_demo_sycl_run.md`](program_temperature_demo_sycl_run.md)。
## 1. 编译

```bash
cd /home/ubuntu/software/SPHinXsys-SPDB
source /opt/intel/oneapi/setvars.sh
# 若尚未把新算例编进 build（首次或 CMake 未识别新目录时）：
cmake -S . -B build-cpu-icpx -DSPHINXSYS_USE_SYCL=OFF
cmake --build build-cpu-icpx --target test_3d_spdb_programmed_temperature -j"$(nproc)"
```

二进制：

```text
build-cpu-icpx/tests/extra_source_and_tests/3d_examples/test_3d_spdb_programmed_temperature/bin/test_3d_spdb_programmed_temperature
```

工作目录必须是该 `bin/`（以便读到 `input/` 下 STL/CSV）。

```bash
BIN=build-cpu-icpx/tests/extra_source_and_tests/3d_examples/test_3d_spdb_programmed_temperature/bin
cd "$BIN"
```

## 2. 建议运行顺序

```bash
# 温度插值 + 门函数单元
./test_3d_spdb_programmed_temperature --temperature_unit

# STL 粒子 / 体积 / 初始空腔
./test_3d_spdb_programmed_temperature --geometry_only

# F1：20℃、5 MPa 弹性平衡
./test_3d_spdb_programmed_temperature --F1 --dp 0.001

# F0：自由热膨胀（先短时）
./test_3d_spdb_programmed_temperature --F0 --dp 0.001 --end_time 60 --macro_dt 10
./test_3d_spdb_programmed_temperature --F0 --dp 0.001 --end_time 300 --macro_dt 10

# F2：升温 + 5 MPa，关非弹性
./test_3d_spdb_programmed_temperature --F2 --dp 0.001 --end_time 60 --macro_dt 10
./test_3d_spdb_programmed_temperature --F2 --dp 0.001 --end_time 300 --macro_dt 10

# F3：主 demo（门控 Arrhenius）
./test_3d_spdb_programmed_temperature --F3 --dp 0.001 --end_time 60 --macro_dt 10
./test_3d_spdb_programmed_temperature --F3 --dp 0.001 --end_time 300 --macro_dt 10
```

通过后再拉长（成本可能数小时级，建议 `nohup`）：

```bash
nohup ./test_3d_spdb_programmed_temperature --F3 --dp 0.001 --end_time 1200 --macro_dt 10 \
  > run_F3_1200.log 2>&1 &
echo $! > job.pid

nohup ./test_3d_spdb_programmed_temperature --F3 --dp 0.001 --end_time 1800 --macro_dt 10 \
  > run_F3_1800.log 2>&1 &
echo $! > job.pid
```

门函数敏感性（900℃ 硬开关）：

```bash
./test_3d_spdb_programmed_temperature --F3 --hard_gate_900 --dp 0.001 --end_time 300 --macro_dt 10
```

宏步减半：

```bash
./test_3d_spdb_programmed_temperature --F3 --dp 0.001 --end_time 300 --macro_dt 5 --max_macro_dt 15
```

## 3. 监控与产物

```bash
cat hold_status.txt
cat hold_heartbeat
kill -0 "$(cat job.pid)" && echo alive   # 若用了后台 pid
tail -f run_F3_1800.log
```

主要输出（均在 `bin/`）：

| 文件 | 含义 |
|------|------|
| `program_temperature_observables.csv` | 逐步观测 |
| `spdb_program_temp_*_metrics.txt` | 终态摘要 |
| `spdb_program_temp_geometry_metrics.txt` | 几何校核 |
| `hold_status.txt` / `hold_heartbeat` | 断线可查进度 |
| VTK（库默认写出） | 关键帧；**不要**大批量提交进 git |

## 4. 口径提醒

- 拆开报告：自由热胀（F0）/ 加载弹性（F1）/ F3−F2 非弹性差
- 穿透用 μm 绝对值对比孔变化；勿把热胀或弹性闭合标成蠕变塌陷
- 未冷却 ≠ 室温永久孔尺寸；末态是 950℃ 有载状态
- STL **不要再 ×0.001**（已是米）

任务书与输入包：`docs/spdb/program_temperature_demo_2026-11-16/`
复用审计：`docs/spdb/reuse_audit.md`
