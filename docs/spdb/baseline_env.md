# SPDB baseline environment

**Date:** 2026-10-15
**Worktree:** `/home/ubuntu/software/SPHinXsys-SPDB`
**Branch:** `feature/spdb-thermo-viscoplastic` (tracks `origin/master`)
**Baseline SHA:** `4e686ee357de1cb798ac47eb760dd79e2fb407ff`
**Upstream remote:** `git@github.com:Xiangyu-Hu/SPHinXsys.git`
**Primary worktree (untouched):** `/home/ubuntu/software/SPHinXsys` on `master`

## Build used for first-round verification

| Item | Value |
|------|--------|
| Build dir | `build-cpu-icpx` |
| Compiler | Intel oneAPI `icpx` 2025.2.2 |
| SYCL | **OFF** (CPU reference; SYCL path not verified this round) |
| CMake generator | Ninja |
| Target | `test_3d_spdb_geometry` |
| GPU present | NVIDIA GeForce RTX 3090 Ti (available, **not used** this round) |

### Configure command

```bash
source /opt/intel/oneapi/setvars.sh
cmake -G Ninja -B build-cpu-icpx \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER=icx \
  -DCMAKE_CXX_COMPILER=icpx \
  -DSPHINXSYS_USE_SYCL=OFF \
  -DSPHINXSYS_BUILD_2D_EXAMPLES=OFF \
  -DSPHINXSYS_BUILD_OPTIMIZATION_EXAMPLES=OFF \
  -DSPHINXSYS_BUILD_PYTHON_INTERFACE=OFF \
  -DSPHINXSYS_BUILD_UNIT_TESTS=OFF \
  -DSPHINXSYS_BUILD_EXTRA_SOURCE_AND_TESTS=ON \
  -DSPHINXSYS_BUILD_MODULES=OFF \
  -DSPHINXSYS_3D=ON \
  -DSPHINXSYS_BUILD_3D_EXAMPLES=ON \
  -DSPHINXSYS_BUILD_TESTS=ON
cmake --build build-cpu-icpx --target test_3d_spdb_geometry -j$(nproc)
```

2026-10-18 起，SPDB 源码在 `tests/extra_source_and_tests/extra_src`，算例在 `tests/extra_source_and_tests/3d_examples/test_3d_spdb_geometry`。配置必须打开 `SPHINXSYS_BUILD_EXTRA_SOURCE_AND_TESTS`。可执行文件路径变为 `build-cpu-icpx/tests/extra_source_and_tests/3d_examples/test_3d_spdb_geometry/bin/test_3d_spdb_geometry`。首轮当时把算例放在 `tests/3d_examples` 且该开关为 OFF，那个目录已经移除。

### Notes

- GNU `g++` configure failed on this host due to Eigen `maybe-uninitialized` promoted by `-Werror`.
- Submodule `3rd_party/wasmtbb` initialized (empty practical impact for this case).
- System deps: Boost, Eigen3, Simbody, TBB, GTest, fmt from distro packages (no vcpkg on this server).
- **Not run this round:** SYCL/GPU solid twin, creep, electro-thermal, LLM interface.
