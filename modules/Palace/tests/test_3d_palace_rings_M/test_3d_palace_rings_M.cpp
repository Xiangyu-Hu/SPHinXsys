#include <iostream>

#include "palace_magnetostatic_interface.hpp"
#include "palace/utils/communication.hpp"

int main(int argc, char *argv[])
{
    palace::Mpi::Init(argc, argv);

    std::string json_path = (argc > 1) ? argv[1] : "data/rings.json";

    sphinxsys_palace::MagnetostaticCase problem(
        json_path, palace::Mpi::World(), /*verbose=*/true);

    problem.Run();

    int n_src = problem.GetNumSources();
    std::cout << "[SphinxSys-Palace] number of sources = " << n_src << std::endl;

    if (n_src > 0)
    {
        // 1) True DOFs of B (for own processing)
        const palace::Vector &B0 = problem.GetBTrueDofs(0);
        std::cout << "  B(true dofs) size = " << B0.Size() << std::endl;

        // 2) GridFunction for B, ready for evaluation / VTK output
        const mfem::ParGridFunction &Bg = problem.GetBField(0);

        // 3) Save as VTU/PVTU
        problem.SaveBFieldParaView("rings_B_source0", "postpro/sphinxsys", 0);
    }

    return 0;
}
