# 需要你执行的命令（按顺序）

物理时长以 Duration V2 为准。A1 L1 GPU FP32 T=4 已通过。

---

## A) 重编译现有 GPU build（分辨率 env + VTP 收紧）

```bash
cd /home/yongchuan/SPHinXsys-gpu-paper/build
ninja test_2d_gpu_paper_a1_poiseuille 2>&1 | tee /home/yongchuan/SPHinXsys-gpu-paper/docs/gpu_paper_discussion/logs/build_a1_gpu_res_env.txt
```

---

## B) 新建 CPU FP32 构建（必须带 vcpkg toolchain）

上次失败原因：未传 `CMAKE_TOOLCHAIN_FILE`，找不到 Simbody。

```bash
ROOT=/home/yongchuan/SPHinXsys-gpu-paper
VCPKG=$ROOT/../vcpkg/scripts/buildsystems/vcpkg.cmake
mkdir -p $ROOT/docs/gpu_paper_discussion/logs
rm -rf $ROOT/build_cpu_fp32
mkdir -p $ROOT/build_cpu_fp32
tmux new-session -d -s gp_build_cpu32 "
  set -o pipefail
  source /opt/intel/oneapi/setvars.sh --force
  cd $ROOT/build_cpu_fp32
  cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=icpx \
    -DCMAKE_TOOLCHAIN_FILE=$VCPKG \
    -DSPHINXSYS_USE_SYCL=OFF -DSPHINXSYS_USE_FLOAT=ON \
    -DSPHINXSYS_BUILD_EXTRA_SOURCE_AND_TESTS=ON \
    -DSPHINXSYS_BUILD_2D_EXAMPLES=ON -DSPHINXSYS_BUILD_3D_EXAMPLES=ON \
    2>&1 | tee $ROOT/docs/gpu_paper_discussion/logs/build_a1_cpu_fp32.txt
  ninja test_2d_gpu_paper_a1_poiseuille 2>&1 | tee -a $ROOT/docs/gpu_paper_discussion/logs/build_a1_cpu_fp32.txt
  echo EXIT:\${PIPESTATUS[0]} | tee -a $ROOT/docs/gpu_paper_discussion/logs/build_a1_cpu_fp32.txt
"
```

查看：`tail -n 40 $ROOT/docs/gpu_paper_discussion/logs/build_a1_cpu_fp32.txt`

---

## C) 新建 CPU FP64 构建

```bash
ROOT=/home/yongchuan/SPHinXsys-gpu-paper
VCPKG=$ROOT/../vcpkg/scripts/buildsystems/vcpkg.cmake
rm -rf $ROOT/build_cpu_fp64
mkdir -p $ROOT/build_cpu_fp64
tmux new-session -d -s gp_build_cpu64 "
  set -o pipefail
  source /opt/intel/oneapi/setvars.sh --force
  cd $ROOT/build_cpu_fp64
  cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=icpx \
    -DCMAKE_TOOLCHAIN_FILE=$VCPKG \
    -DSPHINXSYS_USE_SYCL=OFF -DSPHINXSYS_USE_FLOAT=OFF \
    -DSPHINXSYS_BUILD_EXTRA_SOURCE_AND_TESTS=ON \
    -DSPHINXSYS_BUILD_2D_EXAMPLES=ON -DSPHINXSYS_BUILD_3D_EXAMPLES=ON \
    2>&1 | tee $ROOT/docs/gpu_paper_discussion/logs/build_a1_cpu_fp64.txt
  ninja test_2d_gpu_paper_a1_poiseuille 2>&1 | tee -a $ROOT/docs/gpu_paper_discussion/logs/build_a1_cpu_fp64.txt
  echo EXIT:\${PIPESTATUS[0]} | tee -a $ROOT/docs/gpu_paper_discussion/logs/build_a1_cpu_fp64.txt
"
```

---

## D) A1 L1 CPU FP32 正式验证（B 编过后）

```bash
ROOT=/home/yongchuan/SPHinXsys-gpu-paper
BIN=$ROOT/build_cpu_fp32/tests/extra_source_and_tests/2d_examples/test_2d_gpu_paper_a1_poiseuille/bin/test_2d_gpu_paper_a1_poiseuille
RUN=$ROOT/docs/gpu_paper_discussion/results/A1_poiseuille/runs/A1_L1_cpu_fp32_validate_r01
mkdir -p "$RUN" "$ROOT/docs/gpu_paper_discussion/logs"
SESSION=gp_A1_L1_cpu32

tmux new-session -d -s "$SESSION" "
  set -o pipefail
  export GPU_PAPER_RUN_ID=A1_L1_cpu_fp32_validate_r01
  export GPU_PAPER_MODE=validate
  export GPU_PAPER_RUN_DIR=$RUN
  export GPU_PAPER_END_TIME=4.0
  export GPU_PAPER_DH_OVER_DP=20
  export GPU_PAPER_OUTPUT_INTERVAL=4.0
  export GPU_PAPER_ENABLE_VTP=0
  export GPU_PAPER_SCHEDULE_JSON=$ROOT/docs/gpu_paper_discussion/results/A1_poiseuille/schedules/A1_validate_T4_bound.json
  cd \$(dirname $BIN)
  $BIN 2>&1 | tee $ROOT/docs/gpu_paper_discussion/logs/${SESSION}.log
  echo EXIT:\${PIPESTATUS[0]} | tee -a $ROOT/docs/gpu_paper_discussion/logs/${SESSION}.log
"
```

---

## E) A1 L1 CPU FP64 正式验证（C 编过后）

同 D，把路径/`RUN_ID`/`SESSION` 换成 `build_cpu_fp64` 与 `A1_L1_cpu_fp64_validate_r01` / `gp_A1_L1_cpu64`。

---

## 说明（A1）

- `GPU_PAPER_DH_OVER_DP`：L1=20，L2=40，L3=80  
- A 组 CPU 轨迹默认关 VTP（仅 L2 GPU 代表场）  
- 进度查看：`bash docs/gpu_paper_discussion/scripts/check_runs.sh`  
- **看时间步**：脚本已 `tee`；前台跑，或 `tmux attach -t <session>`  

---

## F) 编译 A2（GPU 现有 build，增量）

```bash
cd /home/yongchuan/SPHinXsys-gpu-paper/build
source /opt/intel/oneapi/setvars.sh --force
cmake ..   # 若 ninja 尚无 A2 target
ninja test_2d_gpu_paper_a2_oscillating_beam 2>&1 | tee /home/yongchuan/SPHinXsys-gpu-paper/docs/gpu_paper_discussion/logs/build_a2_gpu.txt
```

CPU（已有 `build_cpu_fp32` / `build_cpu_fp64`）：

```bash
cd /home/yongchuan/SPHinXsys-gpu-paper/build_cpu_fp32 && ninja test_2d_gpu_paper_a2_oscillating_beam
cd /home/yongchuan/SPHinXsys-gpu-paper/build_cpu_fp64 && ninja test_2d_gpu_paper_a2_oscillating_beam
```

（若提示缺 target：先在该目录 `cmake ..` 再 ninja。）

---

## G) A2 短测（可选，T=0.05）

```bash
ROOT=/home/yongchuan/SPHinXsys-gpu-paper
BIN=$ROOT/build/tests/extra_source_and_tests/2d_examples/test_2d_gpu_paper_a2_oscillating_beam/bin/test_2d_gpu_paper_a2_oscillating_beam
source /opt/intel/oneapi/setvars.sh --force
export CUDA_VISIBLE_DEVICES=1
export GPU_PAPER_RUN_ID=A2_smoke_r01 GPU_PAPER_MODE=validate
export GPU_PAPER_RUN_DIR=$ROOT/docs/gpu_paper_discussion/results/A2_oscillating_beam/runs/A2_smoke_r01
export GPU_PAPER_END_TIME=0.05 GPU_PAPER_DH_OVER_DP=10 GPU_PAPER_ENABLE_VTP=0
export GPU_PAPER_SCHEDULE_JSON=$ROOT/docs/gpu_paper_discussion/results/A2_oscillating_beam/schedules/A2_validate_T2_bound.json
mkdir -p "$GPU_PAPER_RUN_DIR"
cd $(dirname $BIN) && $BIN
```

---

## H) A2 L1 GPU 正式（T=2，终端可见时间步）

```bash
bash /home/yongchuan/SPHinXsys-gpu-paper/docs/gpu_paper_discussion/scripts/run_a2_validate.sh L1 gpu_fp32
```

或 tmux：

```bash
tmux new-session -s gp_A2_L1_gpu "bash /home/yongchuan/SPHinXsys-gpu-paper/docs/gpu_paper_discussion/scripts/run_a2_validate.sh L1 gpu_fp32"
```

通用：`run_a2_validate.sh <L1|L2|L3> <gpu_fp32|cpu_fp32|cpu_fp64>`  
（A2 的 DH 含义为 **PH/dp**：L1=10，L2=20，L3=40；仅 L2 GPU 开 VTP。）

---

## I) 编译 A4（GPU / CPU FP32 / CPU FP64，增量）

GPU（已有 `build/`，缺 target 时先 `cmake ..`）：

```bash
cd /home/yongchuan/SPHinXsys-gpu-paper/build
cmake .. >/dev/null
ninja test_2d_gpu_paper_a4_two_phase_dambreak 2>&1 | tee /home/yongchuan/SPHinXsys-gpu-paper/docs/gpu_paper_discussion/logs/build_a4_gpu.txt
echo EXIT:${PIPESTATUS[0]} | tee -a /home/yongchuan/SPHinXsys-gpu-paper/docs/gpu_paper_discussion/logs/build_a4_gpu.txt
```

CPU FP32 / FP64：

```bash
ROOT=/home/yongchuan/SPHinXsys-gpu-paper
cd $ROOT/build_cpu_fp32 && cmake .. >/dev/null && ninja test_2d_gpu_paper_a4_two_phase_dambreak 2>&1 | tee $ROOT/docs/gpu_paper_discussion/logs/build_a4_cpu_fp32.txt
echo EXIT:${PIPESTATUS[0]} | tee -a $ROOT/docs/gpu_paper_discussion/logs/build_a4_cpu_fp32.txt
cd $ROOT/build_cpu_fp64 && cmake .. >/dev/null && ninja test_2d_gpu_paper_a4_two_phase_dambreak 2>&1 | tee $ROOT/docs/gpu_paper_discussion/logs/build_a4_cpu_fp64.txt
echo EXIT:${PIPESTATUS[0]} | tee -a $ROOT/docs/gpu_paper_discussion/logs/build_a4_cpu_fp64.txt
```

## J) A4 L1 validate（T=10，终端可见 N=/Time=）

```bash
bash /home/yongchuan/SPHinXsys-gpu-paper/docs/gpu_paper_discussion/scripts/run_a4_validate.sh L1 gpu_fp32
```

tmux：

```bash
tmux new-session -s gp_A4_L1_gpu "bash /home/yongchuan/SPHinXsys-gpu-paper/docs/gpu_paper_discussion/scripts/run_a4_validate.sh L1 gpu_fp32"
```

通用：`run_a4_validate.sh <L1|L2|L3> <gpu_fp32|cpu_fp32|cpu_fp64>`  
（`GPU_PAPER_DH_OVER_DP` = **LH/dp**：L1=20，L2=40，L3=80；仅 L2 GPU 开 VTP。几何为 CK 源码非论文匹配。）

---

## K) 编译 B2（各端一条；缺 target 时先 cmake）

```bash
cd /home/yongchuan/SPHinXsys-gpu-paper/build && cmake .. >/dev/null && ninja test_3d_gpu_paper_b2_dambreak 2>&1 | tee /home/yongchuan/SPHinXsys-gpu-paper/docs/gpu_paper_discussion/logs/build_b2_gpu.txt; echo EXIT:${PIPESTATUS[0]} | tee -a /home/yongchuan/SPHinXsys-gpu-paper/docs/gpu_paper_discussion/logs/build_b2_gpu.txt
```

```bash
cd /home/yongchuan/SPHinXsys-gpu-paper/build_cpu_fp32 && cmake .. >/dev/null && ninja test_3d_gpu_paper_b2_dambreak 2>&1 | tee /home/yongchuan/SPHinXsys-gpu-paper/docs/gpu_paper_discussion/logs/build_b2_cpu_fp32.txt; echo EXIT:${PIPESTATUS[0]} | tee -a /home/yongchuan/SPHinXsys-gpu-paper/docs/gpu_paper_discussion/logs/build_b2_cpu_fp32.txt
```

```bash
cd /home/yongchuan/SPHinXsys-gpu-paper/build_cpu_fp64 && cmake .. >/dev/null && ninja test_3d_gpu_paper_b2_dambreak 2>&1 | tee /home/yongchuan/SPHinXsys-gpu-paper/docs/gpu_paper_discussion/logs/build_b2_cpu_fp64.txt; echo EXIT:${PIPESTATUS[0]} | tee -a /home/yongchuan/SPHinXsys-gpu-paper/docs/gpu_paper_discussion/logs/build_b2_cpu_fp64.txt
```

## L) B2 短测 / 正式（脚本）

```bash
bash /home/yongchuan/SPHinXsys-gpu-paper/docs/gpu_paper_discussion/scripts/run_b2.sh L2 gpu_fp32 smoke 0.2
```

```bash
bash /home/yongchuan/SPHinXsys-gpu-paper/docs/gpu_paper_discussion/scripts/run_b2.sh L3 gpu_fp32 benchmark
```
