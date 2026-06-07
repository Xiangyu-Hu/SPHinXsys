/**
 * Stage 10.7: impressed-current source divergence diagnostic (div J_s).
 */
#include "electromagnetic_dynamics/diagnostics/aphi_a_gauge_diagnostic_helpers.h"
#include "electromagnetic_dynamics/benchmark/aphi_benchmark_case_ck.hpp"
#include "electromagnetic_dynamics/test_helpers/aphi_lhs_test_helpers.h"

#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics;
using namespace SPH::electromagnetics::benchmark;
using namespace SPH::electromagnetics::test;

int main(int ac, char *av[])
{
    const Real dp_0 = 0.1;
    const Real body_length = 1.0;
    const Real body_height = 1.0;
    const Real body_width = 1.0;
    const Real boundary_width = 3.0 * dp_0;
    const Real sigma = 2.0;
    const Real nu = 1.5;
    const Real impressed_current_amplitude = 5.0;
    const AphiBoxRegion source_region{0.35, 0.65, 0.35, 0.65, 0.35, 0.65};
    const Vecd current_real(0.0, 0.0, 1.0);
    const Vecd current_imag(0.0, 0.0, 0.1);

    AphiLhsTestBody test_body(dp_0, body_length, body_height, body_width, boundary_width, ac, av);
    IOEnvironment io_environment(test_body.sph_system);

    AphiVariableNames names;
    const std::string div_j_real_name = "DivJReal";
    const std::string div_j_imag_name = "DivJImag";

    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_aphi(test_body.body, sigma, nu, names);
    StateDynamics<MainExecutionPolicy, SetAphiMaterialPropertiesCK> set_material(test_body.body, sigma, nu, names.material);
    StateDynamics<MainExecutionPolicy, AssignImpressedCurrentRhsCK> assign_impressed_current(
        test_body.body, names.rhs, source_region, current_real, current_imag, impressed_current_amplitude);

    InteractionDynamicsCK<MainExecutionPolicy, LinearCorrectionMatrix<Inner<WithUpdate>>> linear_correction_matrix(
        test_body.inner());
    InteractionDynamicsCK<MainExecutionPolicy, LinearGradient<Inner<Vecd>>> j_real_gradient(test_body.inner(),
                                                                                            names.rhs.a_real);
    InteractionDynamicsCK<MainExecutionPolicy, LinearGradient<Inner<Vecd>>> j_imag_gradient(test_body.inner(),
                                                                                            names.rhs.a_imag);
    StateDynamics<MainExecutionPolicy, AphiVectorGradientDivergenceCK> div_j_real(
        test_body.body, names.rhs.a_real + "Gradient", div_j_real_name);
    StateDynamics<MainExecutionPolicy, AphiVectorGradientDivergenceCK> div_j_imag(
        test_body.body, names.rhs.a_imag + "Gradient", div_j_imag_name);

    initialize_aphi.exec();
    set_material.exec();
    assign_impressed_current.exec();
    test_body.updateRelations();
    linear_correction_matrix.exec();
    j_real_gradient.exec();
    j_imag_gradient.exec();
    div_j_real.exec();
    div_j_imag.exec();

    BaseParticles &particles = test_body.body.getBaseParticles();
    const size_t total_real_particles = particles.TotalRealParticles();
    const AphiSourceCurrentDivergenceMetrics metrics = hostSourceCurrentDivergenceMetrics(
        particles, names.rhs, div_j_real_name, div_j_imag_name, source_region, total_real_particles);

    const bool passed = metrics.source_j_l2 > TinyReal && metrics.source_region_count > 0;

    std::cout << "test_3d_aphi_ck_impressed_current_divergence_diagnostic"
              << " source_J_L2=" << metrics.source_j_l2 << " source_div_J_L2=" << metrics.source_div_j_l2
              << " source_div_J_relative=" << metrics.source_div_j_relative
              << " source_div_J_Linf=" << metrics.source_div_j_linf
              << " source_region_count=" << metrics.source_region_count
              << " total_particles=" << metrics.total_particle_count << " passed=" << (passed ? 1 : 0) << std::endl;

    return passed ? 0 : 1;
}
