source /opt/intel/oneapi/setvars.sh --include-intel-llvm && cmake   -G "Unix Makefiles"                                                     \
    -D CMAKE_BUILD_TYPE=Debug                                                   \
    -D CMAKE_C_COMPILER=icx -D CMAKE_CXX_COMPILER=icpx                          \
    -D CMAKE_TOOLCHAIN_FILE="$HOME/vcpkg/scripts/buildsystems/vcpkg.cmake"      \
    -D CMAKE_C_COMPILER_LAUNCHER=ccache -D CMAKE_CXX_COMPILER_LAUNCHER=ccache   \
    -D SPHINXSYS_USE_SYCL=ON                                                    \
    -D SPHINXSYS_SYCL_TARGETS=spir64_x86_64                                            \
    -S .                                                                        \
    -B ./build
cmake --build build/ --target test_2d_dambreak_sycl -j$(nproc)
