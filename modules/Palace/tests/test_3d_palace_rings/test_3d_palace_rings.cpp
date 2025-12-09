#include <iostream>
#include "sphinxsys_palace_interface.hpp"

int main(int argc, char* argv[])
{
    std::string json = "data/rings.json";
    if (argc > 1)
        json = argv[1];

    std::cout << "[SPHinXsys-Palace] Running rings test …\n";
    std::cout << "  Using config: " << json << "\n";

    int status = sphinxsys_palace::RunPalaceFromJson(json);

    if (status == 0)
        std::cout << "[OK] Palace simulation completed.\n";
    else
        std::cout << "[ERROR] Palace simulation failed, status = " << status << "\n";

    return status;
}
