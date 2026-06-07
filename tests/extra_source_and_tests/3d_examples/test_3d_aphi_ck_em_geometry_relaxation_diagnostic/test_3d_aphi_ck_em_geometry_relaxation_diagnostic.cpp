/** Stage 10.14-B P7: three-body TEAM7 lattice/contact geometry scaffold (no EM solve). */
#include "electromagnetic_dynamics/diagnostics/aphi_em_geometry_relaxation_diagnostic_helpers.h"

#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics::test;

int main(int ac, char *av[])
{
    const AphiEmGeometryRelaxationSpec spec;
    const AphiEmGeometryRelaxationMetrics metrics = runEmGeometryRelaxationDiagnostic(ac, av, spec);
    const bool passed = emGeometryRelaxationDiagnosticPassed(metrics);
    std::cout << "test_3d_aphi_ck_em_geometry_relaxation_diagnostic"
              << " passed=" << (passed ? 1 : 0) << " particles=" << metrics.particles
              << " min_body_particles=" << metrics.min_body_particle_count
              << " max_body_particles=" << metrics.max_body_particle_count
              << " geometry_scaffold_only=1" << std::endl;
    return passed ? 0 : 1;
}
