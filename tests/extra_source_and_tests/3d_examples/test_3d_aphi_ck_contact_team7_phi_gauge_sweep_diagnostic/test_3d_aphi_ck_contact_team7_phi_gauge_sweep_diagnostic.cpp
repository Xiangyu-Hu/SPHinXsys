/**
 * Stage 10.5/10.6: three-body TEAM7 phi gauge penalty sweep + divA diagnostic.
 */
#include "electromagnetic_dynamics/diagnostics/aphi_contact_phi_gauge_sweep_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_gmres_benchmark_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_team7_contact_test_helpers.h"

#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics;
using namespace SPH::electromagnetics::benchmark;
using namespace SPH::electromagnetics::test;

namespace
{

inline AphiTeam7PlateObservables runMonolithicTeam7PlateObservables(int ac, char *av[], Real dp_0, Real body_length,
                                                                    Real body_height, Real body_width,
                                                                    Real boundary_width, Real core_shell, Real omega,
                                                                    Real phi_gauge_penalty, Real tolerance,
                                                                    UnsignedInt restart_dimension,
                                                                    UnsignedInt max_outer_iterations,
                                                                    const Vecd &coil_current_real,
                                                                    const Vecd &coil_current_imag,
                                                                    Real impressed_current_amplitude)
{
    const AphiTeam7LikeUnitBoxLayout layout = buildTeam7LayoutForBox(body_length, body_height, body_width);
    AphiLhsTestBody test_body(dp_0, body_length, body_height, body_width, boundary_width, ac, av);
    AphiVariableNames names;
    AphiJouleHeatingFieldNames joule_names;

    AphiLhsAssemblyOptions options;
    options.omega = omega;
    options.use_phi_gauge_penalty = true;
    options.phi_gauge_penalty = phi_gauge_penalty;
    const AphiMatrixFreeSolverOptions solver_options =
        defaultMatrixFreeGMRESOptions(tolerance, restart_dimension, max_outer_iterations);

    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_aphi(test_body.body, layout.air.sigma,
                                                                                  layout.air.nu, names);
    StateDynamics<MainExecutionPolicy, AssignTeam7LikeRegionMaterialsCK> assign_material(test_body.body, layout,
                                                                                         names.material);
    RegisterAphiJouleHeatingFieldsCK register_joule(test_body.body, joule_names);
    StateDynamics<MainExecutionPolicy, AssignImpressedCurrentRhsCK> assign_coil_source(
        test_body.body, names.rhs, layout.coil, coil_current_real, coil_current_imag, impressed_current_amplitude);
    StateDynamics<MainExecutionPolicy, AphiZeroBlockCK> zero_solution(test_body.body, names.solution);

    AphiMatrixFreeSolveCK<MainExecutionPolicy> solver(test_body.body, test_body.inner(), names, options,
                                                        solver_options);
    InteractionDynamicsCK<MainExecutionPolicy, AphiComputeScalarPhiGradientCK<Inner<>>> compute_grad_phi(
        test_body.inner(), names.solution, joule_names);
    StateDynamics<MainExecutionPolicy, AphiComputeFrequencyElectricFieldCK> compute_electric_field(
        test_body.body, omega, names.solution, joule_names);
    StateDynamics<MainExecutionPolicy, AphiComputeJouleHeatSourceCK> compute_joule(test_body.body, names.material,
                                                                                     joule_names);

    initialize_aphi.exec();
    assign_material.exec();
    assign_coil_source.exec();
    test_body.updateRelations();
    zero_solution.exec();
    test_body.updateRelations();
    (void)solver.solve();
    compute_grad_phi.exec();
    compute_electric_field.exec();
    compute_joule.exec();

    return hostTeam7PlateObservablesFromBody(test_body.body, names, joule_names, layout, body_length, body_height,
                                             body_width, core_shell);
}

inline void printTeam7Row(const char *test_name, const AphiTeam7PhiGaugeSweepRow &row)
{
    std::cout << test_name << " penalty=" << row.phi_gauge_penalty << " converged=" << (row.converged ? 1 : 0)
              << " outer=" << row.outer_iterations << " global_true_rel=" << row.global_true_rel
              << " max_bodywise_true_rel=" << row.max_bodywise_true_rel
              << " plate_joule_power_gap=" << row.plate_joule_power_gap << " plate_j_L2_gap=" << row.plate_j_L2_gap
              << " div_A_relative=" << row.global_div_a.div_a_relative
              << " plate_div_A_relative=" << row.plate_div_a.div_a_relative
              << " div_A_level=" << divAGaugeDiagnosticLevel(row.global_div_a.div_a_relative) << std::endl;
}

} // namespace

int main(int ac, char *av[])
{
    const Real dp_0 = 0.1;
    const Real body_length = 1.0;
    const Real body_height = 1.0;
    const Real body_width = 1.0;
    const Real boundary_width = 3.0 * dp_0;
    const Real core_shell = 2.0 * dp_0;
    const Real omega = 1.25;
    const Real tolerance = 5.0e-4;
    const Real impressed_current_amplitude = 8.0;
    const UnsignedInt restart_dimension = 50;
    const UnsignedInt max_outer_iterations = 100;
    const Vecd coil_current_real(0.0, 0.0, 1.0);
    const Vecd coil_current_imag(0.0, 0.0, 0.0);

    const Real penalty_values[] = {10.0, 30.0, 100.0, 300.0};
    std::cout << "test_3d_aphi_ck_contact_team7_phi_gauge_sweep_diagnostic begin" << std::endl;
    for (const Real penalty : penalty_values)
    {
        const AphiTeam7PlateObservables mono_plate_at_penalty = runMonolithicTeam7PlateObservables(
            ac, av, dp_0, body_length, body_height, body_width, boundary_width, core_shell, omega, penalty,
            tolerance, restart_dimension, max_outer_iterations, coil_current_real, coil_current_imag,
            impressed_current_amplitude);
        const AphiTeam7PhiGaugeSweepRow row = runTeam7PhiGaugeSweepRow(
            ac, av, dp_0, body_length, body_height, body_width, boundary_width, core_shell, omega, penalty, tolerance,
            coil_current_real, coil_current_imag, impressed_current_amplitude, mono_plate_at_penalty);
        printTeam7Row("test_3d_aphi_ck_contact_team7_phi_gauge_sweep_diagnostic", row);
    }

    std::cout << "test_3d_aphi_ck_contact_team7_phi_gauge_sweep_diagnostic completed_rows=4 passed=1" << std::endl;
    return 0;
}
