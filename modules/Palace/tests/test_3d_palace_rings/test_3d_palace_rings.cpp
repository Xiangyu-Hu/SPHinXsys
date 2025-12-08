// modules/Palace/tests/test_3d_palace_rings/test_3d_palace_rings.cpp
#include <iostream>
#include <string>

#include "palace_interface.hpp"
#include "palace/utils/communication.hpp"
#include "palace/linalg/hypre.hpp"
#include "palace/linalg/slepc.hpp"
#include "palace/fem/libceed/ceed.hpp"

#ifndef PALACE_JSON_PATH
#define PALACE_JSON_PATH "data/rings.json"
#endif

int main(int argc, char* argv[])
{
 
    // ---------------- MPI initialization ----------------
    // Use Palace's MPI wrapper, which internally calls MPI_Init.
    // SPHinXsys will follow the same pattern when integrating.
    palace::Mpi::Init(argc, argv);


    // Select JSON config file: default (from CMake macro) or command-line override.
    std::string json_path = PALACE_JSON_PATH;
    if (argc > 1)
    {
        json_path = argv[1];
    }

    std::cout << "[SphinxSys-Palace] Magnetostatic rings test\n";
    std::cout << "  Using config file: " << json_path << std::endl;


    // ---------------- Run Palace solver ----------------
    // Calls the interface function defined in palace_interface.cpp.
    int status = sphinxsys_palace::RunMagnetostaticCase(json_path, /*verbose=*/true);

    // ---------------- Report result ----------------
    if (status != 0)
    {
        std::cerr << "[SphinxSys-Palace] Rings test FAILED with status "
                  << status << std::endl;
    }
    else
    {
        std::cout << "[SphinxSys-Palace] Rings test finished successfully."
                  << std::endl;
    }


    // NOTE:
    // We intentionally do NOT call Mpi::Finalize() here.
    // Palace may perform further internal finalization when destructors run.
    // Uncomment only if running this in a standalone environment.
    // palace::Mpi::Finalize();

    return status;
}
