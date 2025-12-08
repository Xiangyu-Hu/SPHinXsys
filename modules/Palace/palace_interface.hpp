// modules/Palace/palace_interface.hpp
#pragma once

#include <string>

// 前向声明 MPI_Comm 即可，不一定要在这里包含 <mpi.h>
// 也可以 #include <mpi.h>，随你
#include <mpi.h>

namespace sphinxsys_palace
{
    // 单个 Palace 磁静态算例：
    // - config_file: Palace 的 json 配置路径
    // - verbose: 是否打印 Palace 自己的一些信息
    //
    // 注意：
    //   - 这里内部使用 MPI_COMM_WORLD，不做 MPI_Init / MPI_Finalize
    //   - MPI 的生命周期由外层 main 负责（比如 test_3d_palace_rings.cpp）
    int RunMagnetostaticCase(const std::string &config_file,
                             bool verbose = true);
}
