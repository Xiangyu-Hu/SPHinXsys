// modules/Palace/palace_interface.hpp
#pragma once

#include <string>

// We only need MPI_Comm here, so including <mpi.h> is sufficient.
// MPI is expected to be initialized and finalized by the caller (e.g. SPHinXsys).
#include <mpi.h>

namespace sphinxsys_palace
{
   
    // Run the Palace magnetostatic solver (e.g. rings case) using the given
    // JSON configuration file. The function assumes:
    //   - MPI_Init has already been called by the caller,
    //   - config_file is a valid Palace input JSON,
    //   - all required Palace libraries and data files are accessible.
    //
    // Parameters:
    //   config_file : path to Palace JSON configuration (e.g. "data/rings.json")
    //   verbose     : if true, print information messages on the root rank
    //
    // Returns:
    //   0 on success (currently always returns 0 if no fatal error occurs).
    
    int RunMagnetostaticCase(const std::string &config_file,
                             bool verbose = true);
}
