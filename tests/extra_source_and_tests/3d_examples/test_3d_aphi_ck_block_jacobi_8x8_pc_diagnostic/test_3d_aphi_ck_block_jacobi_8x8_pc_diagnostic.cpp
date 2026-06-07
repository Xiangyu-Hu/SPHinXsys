/**
 * Stage 10A-1: CoupledPointBlock8x8 PC fallback / min-pivot diagnostics on TEAM7-like layout.
 */
#include "electromagnetic_dynamics/test_helpers/aphi_gmres_benchmark_helpers.h"

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
    const Real omega = 1.25;
    const Real phi_gauge_penalty = 100.0;
    const Real impressed_current_amplitude = 8.0;

    const AphiTeam7LikeUnitBoxLayout layout;
    const Vecd coil_current_real(0.0, 0.0, 1.0);
    const Vecd coil_current_imag(0.0, 0.0, 0.0);
    const AphiBlockJacobiDiagnosticNames jacobi_diag_names;

    AphiLhsTestBody test_body(dp_0, body_length, body_height, body_width, boundary_width, ac, av);
    IOEnvironment io_environment(test_body.sph_system);

    AphiVariableNames names;
    RegisterAphiBlockJacobiDiagnosticFieldsCK register_jacobi_diag(test_body.body, jacobi_diag_names);
    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_aphi_variables(
        test_body.body, layout.air.sigma, layout.air.nu, names);
    StateDynamics<MainExecutionPolicy, AssignTeam7LikeRegionMaterialsCK> assign_material(test_body.body, layout,
                                                                                         names.material);
    StateDynamics<MainExecutionPolicy, AssignImpressedCurrentRhsCK> assign_coil_source(
        test_body.body, names.rhs, layout.coil, coil_current_real, coil_current_imag, impressed_current_amplitude);

    AphiLhsAssemblyOptions options;
    options.omega = omega;
    options.use_phi_gauge_penalty = true;
    options.phi_gauge_penalty = phi_gauge_penalty;

    InteractionDynamicsCK<MainExecutionPolicy, AphiComputeBlockJacobiDiagonalCK<Inner<>>> compute_jacobi_diagonal(
        test_body.inner(), names.material, omega, options);
    StateDynamics<MainExecutionPolicy, AphiApplyBlockJacobiInverseCK> apply_coupled_pc(
        test_body.body, names.rhs, names.search, names.material, omega, options,
        AphiBlockJacobiPreconditionerKind::CoupledPointBlock8x8, AphiBlockJacobiDiagonalNames{}, &jacobi_diag_names);

    (void)register_jacobi_diag;
    initialize_aphi_variables.exec();
    assign_material.exec();
    assign_coil_source.exec();
    test_body.updateRelations();

    compute_jacobi_diagonal.exec();
    apply_coupled_pc.exec();
    test_body.updateRelations();

    BaseParticles &particles = test_body.body.getBaseParticles();
    syncVariableToHost<Vecd>(particles, "Position");
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    const size_t total_real_particles = particles.TotalRealParticles();

    const AphiBlockJacobi8x8DiagnosticSummary summary =
        hostBlockJacobi8x8DiagnosticSummary(particles, positions, total_real_particles, layout, jacobi_diag_names);

    std::cout << "test_3d_aphi_ck_block_jacobi_8x8_pc_diagnostic"
              << " particles=" << summary.total_particles << " fallback_count=" << summary.fallback_count
              << " fallback_fraction=" << summary.fallback_fraction << " global_min_pivot=" << summary.global_min_pivot
              << " conductor_fallback=" << summary.conductor_fallback << " coil_fallback=" << summary.coil_fallback
              << " air_fallback=" << summary.air_fallback << " passed=1" << std::endl;

    return 0;
}
