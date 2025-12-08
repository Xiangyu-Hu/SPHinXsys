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
   
    palace::Mpi::Init(argc, argv);

    std::string json_path = PALACE_JSON_PATH;
    if (argc > 1)
    {
        json_path = argv[1];
    }

    std::cout << "[SphinxSys-Palace] Magnetostatic rings test\n";
    std::cout << "  Using config file: " << json_path << std::endl;


    int status = sphinxsys_palace::RunMagnetostaticCase(json_path, /*verbose=*/true);

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


    // palace::Mpi::Finalize();

    return status;
}
