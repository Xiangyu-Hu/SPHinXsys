#include "palace_magnetostatic_interface.hpp"
#include "palace/utils/communication.hpp"

int main(int argc, char *argv[])
{
    palace::Mpi::Init(argc, argv);

    std::string json_path = (argc > 1) ? argv[1] : "data/rings.json";

    sphinxsys_palace::MagnetostaticCase problem(
        json_path, palace::Mpi::World(), /*verbose=*/true);

    problem.Run();

    // Example 1: export B field for the first source to CSV
    problem.ExportBFieldToCsv("rings_B_source0.csv", /*source_index=*/0);

    // Example 2: access B as a ParGridFunction and sample at custom coordinates
    const mfem::ParGridFunction &B = problem.GetBField(0);

    // Suppose you have particle positions from SPHinXsys:
    // std::vector<std::array<double,3>> particle_pos = ...;
    // for each position, use B.GetVectorValue(...) to interpolate.

    return 0;
}
