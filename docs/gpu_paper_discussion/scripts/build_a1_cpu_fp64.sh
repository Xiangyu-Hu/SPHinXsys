#!/usr/bin/env bash
# Clean rebuild A1 CPU FP64 with Ninja (avoids leftover Unix Makefiles cache).
ROOT=/home/yongchuan/SPHinXsys-gpu-paper
VCPKG=$ROOT/../vcpkg/scripts/buildsystems/vcpkg.cmake
LOG=$ROOT/docs/gpu_paper_discussion/logs/build_a1_cpu_fp64.txt
BUILD=$ROOT/build_cpu_fp64
mkdir -p "$(dirname "$LOG")"
{
  echo "START:$(date -Is)"
  echo "cleaning $BUILD"
  rm -rf "$BUILD"
  mkdir -p "$BUILD"
} >"$LOG" 2>&1

# setvars may alter shell options; source after clean
# shellcheck disable=SC1091
source /opt/intel/oneapi/setvars.sh --force >>"$LOG" 2>&1 || true
set -uo pipefail

cd "$BUILD" || { echo "EXIT:2" | tee -a "$LOG"; exit 2; }
cmake .. -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=icpx \
  -DCMAKE_TOOLCHAIN_FILE="$VCPKG" \
  -DSPHINXSYS_USE_SYCL=OFF \
  -DSPHINXSYS_USE_FLOAT=OFF \
  -DSPHINXSYS_BUILD_EXTRA_SOURCE_AND_TESTS=ON \
  -DSPHINXSYS_BUILD_2D_EXAMPLES=ON \
  -DSPHINXSYS_BUILD_3D_EXAMPLES=ON \
  >>"$LOG" 2>&1
cmake_rc=$?
if [[ "$cmake_rc" -ne 0 ]]; then
  echo "EXIT:${cmake_rc}" | tee -a "$LOG"
  exit "$cmake_rc"
fi
if [[ ! -f "$BUILD/build.ninja" ]]; then
  echo "ERROR: build.ninja missing after cmake" | tee -a "$LOG"
  echo "EXIT:3" | tee -a "$LOG"
  exit 3
fi
ninja test_2d_gpu_paper_a1_poiseuille >>"$LOG" 2>&1
echo "EXIT:$?" | tee -a "$LOG"
