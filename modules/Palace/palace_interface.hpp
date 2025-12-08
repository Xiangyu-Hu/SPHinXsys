// modules/Palace/palace_interface.hpp
#pragma once

#include <string>

// 前向声明 MPI_Comm 即可，不一定要在这里包含 <mpi.h>
// 也可以 #include <mpi.h>，随你
#include <mpi.h>

namespace sphinxsys_palace
{
   
    int RunMagnetostaticCase(const std::string &config_file,
                             bool verbose = true);
}
