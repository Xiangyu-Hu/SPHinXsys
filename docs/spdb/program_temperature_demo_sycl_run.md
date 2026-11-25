# 程序温度双槽 Demo — SYCL/GPU 编译与运行

CPU 版（`build-cpu-icpx`）保留作参考。GPU 是**独立算例**：

`tests/tests_sycl/3d_examples/test_3d_spdb_programmed_temperature_sycl/`

接触用 **`ContactPenaltyForceCK`（沿 e_ij，对齐 CPU `ContactForce`）**，**不用** 库内 `RepulsionForceCK`（齐平面已验证会炸）。

## 为何不能只给 CPU case 开 SYCL 开关

- `SPHINXSYS_USE_SYCL=ON` 时默认不编 `extra_source_and_tests` 里的 CPU 算例
- 旧 `ContactForce` 是 host/TBB，不会上 GPU
- 已改 CMake：SYCL 构建仍会编 `extra_src`（材料 + CK 接触），供本 GPU 算例链接

## 配置与编译（独立目录，勿覆盖 CPU 树）

```bash
cd /home/ubuntu/software/SPHinXsys-SPDB
source /opt/intel/oneapi/setvars.sh

cmake -G Ninja -B build-sycl-gpu \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER=icx \
  -DCMAKE_CXX_COMPILER=icpx \
  -DSPHINXSYS_USE_SYCL=ON \
  -DSPHINXSYS_SYCL_TARGETS=nvptx64-nvidia-cuda \
  -DSPHINXSYS_BUILD_2D_EXAMPLES=OFF \
  -DSPHINXSYS_BUILD_OPTIMIZATION_EXAMPLES=OFF \
  -DSPHINXSYS_BUILD_PYTHON_INTERFACE=OFF \
  -DSPHINXSYS_BUILD_UNIT_TESTS=OFF \
  -DSPHINXSYS_BUILD_MODULES=OFF \
  -DSPHINXSYS_3D=ON \
  -DSPHINXSYS_BUILD_3D_EXAMPLES=ON \
  -DSPHINXSYS_BUILD_EXTRA_SOURCE_AND_TESTS=ON \
  -DSPHINXSYS_BUILD_TESTS=ON

cmake --build build-sycl-gpu --target test_3d_spdb_programmed_temperature_sycl -j"$(nproc)"
```

注意：

- `SPHINXSYS_BUILD_3D_EXAMPLES` 必须为 **ON**，否则会出现 `unknown target 'test_3d_spdb_programmed_temperature_sycl'`。
- 终端里 `CL/sycl.hpp is deprecated` 只是警告；真正会导致失败的是 **libfmt 9 + icpx SYCL** 对 `long double` 的特化缺失。仓库已用 `cmake/fmt_sycl_fix/` 阴影头修复；需重新 `cmake -B build-sycl-gpu ...` 后再编译。
- 若本机 CUDA 目标名不同，用 `sycl-ls` 核对后改 `SPHINXSYS_SYCL_TARGETS`。## 运行

```bash
BIN=build-sycl-gpu/tests/tests_sycl/3d_examples/test_3d_spdb_programmed_temperature_sycl/bin
cd "$BIN"

./test_3d_spdb_programmed_temperature_sycl --geometry_only
./test_3d_spdb_programmed_temperature_sycl --F1 --dp 0.001
./test_3d_spdb_programmed_temperature_sycl --F2 --dp 0.001 --end_time 60 --macro_dt 10
./test_3d_spdb_programmed_temperature_sycl --F3 --dp 0.001 --end_time 60 --macro_dt 10
```

### 加长阶梯（同一 `bin/`，**串行**，勿并行抢 `output/`）

温度曲线：0→20℃ … 900 s→700℃ … **~1080 s 过 850℃ 门** … 1200 s→950℃ … 1800 s 保温。

**阻尼说明**（相对 `test_3d_sliding`）：
- 现已接入 **host** `DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner>>`（与 sliding / CPU SPDB 同族；**不是**人工体积粘性）
- 流程：device 声学步 + `VelocityDecayCK`；每 `--host_damp_stride` 步（默认 10）同步 pos/vel → host CLL → split damp → 回写 vel
- 热膨胀后先做若干次 host damp；`vmax` 升高时额外加 damp（不清零速度）
- **T≥700℃** 仍自动缩小宏步（5 / 2.5 / 1 s）

```bash
# 编译（在仓库根）
source /opt/intel/oneapi/setvars.sh
cmake --build build-sycl-gpu --target test_3d_spdb_programmed_temperature_sycl -j"$(nproc)"

BIN=build-sycl-gpu/tests/tests_sycl/3d_examples/test_3d_spdb_programmed_temperature_sycl/bin
cd "$BIN"

# 1200 s（建议 tmux/nohup）
nohup ./test_3d_spdb_programmed_temperature_sycl --F3 --dp 0.001 --end_time 1200 --macro_dt 10 \
  --relax_cap 1.5e-3 --velocity_decay 5e4 --vmax_calm 2e-3 \
  --physical_viscosity 5e6 --damp_random 0.5 --host_damp_stride 10 \
  > f3_end1200_hostdamp.log 2>&1 &
```

监控：`tail -f hold_status.txt` 或对应 `f3_end*.log`。成功标志：`finished=1` 且 `aborted=0`（不要只看日志里的 `done mode=`）。

与 CPU 算例**不要共用同一 `bin/`**。

## VTP 场（ParaView）

写出：`Temperature`、`EquivalentInelasticStrain`、`VonMisesStress`、`VonMisesStrain`、`Displacement`、`Ren2014EquivalentStress`（含塑性历史的 J2 等效应力）。

- Representation 用 **Points / Point Gaussian**（不是 Surface）
- 应力/应变色标请手动缩到实际量级（应变常为 1e-6～1e-3 量级）
- `--vtp_every 100`：每 100 s 物理时间写一帧（默认）；末态总会再写一次

```bash
# 开跑前清 output，避免 output_backup 秒挂
cd "$BIN"
rm -rf output_backup
if [ -d output ]; then mv output "output_saved_$(date +%Y%m%d_%H%M%S)"; fi
mkdir -p output

./test_3d_spdb_programmed_temperature_sycl --F3 --dp 0.001 --end_time 60 --macro_dt 10 --vtp_every 30
```

## 注意

- SYCL 强制 **float**；与 CPU double 参考会对不齐到最后一位
- 准静态：`VelocityDecayCK` + 稀疏 host `InteractionSplit` 成对阻尼
- 材料仍是 `demo_assumption`；未做传热/电场
- 空腔指标：绝对 z 带失败时回退到槽内 z 极值（GPU 漂移后仍应给出净高）
- 既往长时 F3 曾在升温段 `bad_dt`/`vmax` 中止；先短时（60/300）确认应力场与稳定性，再拉长

CPU 跑法仍见 [`program_temperature_demo_run.md`](program_temperature_demo_run.md)。
