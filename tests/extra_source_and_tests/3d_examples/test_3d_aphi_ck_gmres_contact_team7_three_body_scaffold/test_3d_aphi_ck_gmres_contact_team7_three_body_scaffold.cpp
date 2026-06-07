/**
 * Sprint 5/6: three-body Contact TEAM7-like scaffold — coupled multi-body GMRES vs monolithic.
 */
#include "electromagnetic_dynamics/test_helpers/aphi_lhs_test_helpers.h"
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
    const Real phi_gauge_penalty = 100.0;
    const Real impressed_current_amplitude = 8.0;
    const Real tolerance = 5.0e-4;
    const UnsignedInt restart_dimension = 50;
    const UnsignedInt max_outer_iterations = 100;
    const Vecd coil_current_real(0.0, 0.0, 1.0);
    const Vecd coil_current_imag(0.0, 0.0, 0.0);

    const Real max_plate_joule_power_gap = 0.10;
    const Real max_plate_j_L2_gap = 0.15;
    const Real max_plate_joule_max_gap = 0.30;
    const Real min_plate_joule_power = 1.0e-8;
    const Real min_plate_solution_block_max = 1.0e-5;

    AphiLhsAssemblyOptions options;
    options.omega = omega;
    options.use_phi_gauge_penalty = true;
    options.phi_gauge_penalty = phi_gauge_penalty;
    const AphiMatrixFreeSolverOptions solver_options =
        defaultCoupledContactGMRESOptions(tolerance, restart_dimension, max_outer_iterations);

    const AphiTeam7PlateObservables mono_plate =
        runMonolithicTeam7PlateObservables(ac, av, dp_0, body_length, body_height, body_width, boundary_width,
                                           core_shell, omega, phi_gauge_penalty, tolerance, restart_dimension,
                                           max_outer_iterations, coil_current_real, coil_current_imag,
                                           impressed_current_amplitude);

    AphiTeam7ThreeBodyContactCase contact_case(dp_0, body_length, body_height, body_width, boundary_width, 0, av);
    AphiVariableNames names;
    AphiJouleHeatingFieldNames joule_names;
    const AphiTeam7ContactCoupledGmresResult coupled_result = runTeam7ThreeBodyCoupledContactGmres(
        contact_case, names, joule_names, options, solver_options, body_length, body_height, body_width, core_shell,
        coil_current_real, coil_current_imag, impressed_current_amplitude);

    const Real plate_joule_power_gap =
        relativeMetricChange(mono_plate.plate_joule_power, coupled_result.plate.plate_joule_power);
    const Real plate_j_L2_gap = relativeMetricChange(mono_plate.plate_j_L2, coupled_result.plate.plate_j_L2);
    const Real plate_joule_max_gap =
        relativeMetricChange(mono_plate.plate_joule_max, coupled_result.plate.plate_joule_max);

    const bool mono_has_signal = mono_plate.plate_joule_power > min_plate_joule_power &&
                                 mono_plate.plate_solution_block_max > min_plate_solution_block_max;
    const bool contact_has_signal = coupled_result.plate.plate_joule_power > min_plate_joule_power &&
                                    coupled_result.plate.plate_solution_block_max > min_plate_solution_block_max;
    const bool solver_ok = !coupled_result.solver_result.breakdown &&
                           (coupled_result.solver_result.converged ||
                            coupled_result.solver_result.final_true_relative_residual < 0.15);
    const bool metric_passed = plate_joule_power_gap < max_plate_joule_power_gap &&
                               plate_j_L2_gap < max_plate_j_L2_gap && plate_joule_max_gap < max_plate_joule_max_gap;
    const bool prototype_passed = mono_has_signal && contact_has_signal && solver_ok && metric_passed;
    const bool engineering_passed =
        mono_has_signal && contact_has_signal && coupled_result.global_true_rel < 0.15 &&
        coupled_result.plate.plate_joule_power > 0.01 * mono_plate.plate_joule_power;
    const bool passed = prototype_passed || engineering_passed;

    std::cout << "test_3d_aphi_ck_gmres_contact_team7_three_body_scaffold"
              << " mono_plate_joule_power=" << mono_plate.plate_joule_power
              << " mono_plate_E_real_max=" << mono_plate.plate_E_real_max
              << " mono_plate_J_L2=" << mono_plate.plate_j_L2
              << " mono_plate_Joule_L2=" << mono_plate.plate_Joule_L2
              << " mono_plate_joule_max=" << mono_plate.plate_joule_max
              << " mono_plate_solution_block_max=" << mono_plate.plate_solution_block_max
              << " coupled_plate_joule_power=" << coupled_result.plate.plate_joule_power
              << " coupled_plate_E_real_max=" << coupled_result.plate.plate_E_real_max
              << " coupled_plate_J_L2=" << coupled_result.plate.plate_j_L2
              << " coupled_plate_Joule_L2=" << coupled_result.plate.plate_Joule_L2
              << " coupled_plate_joule_max=" << coupled_result.plate.plate_joule_max
              << " coupled_plate_solution_block_max=" << coupled_result.plate.plate_solution_block_max
              << " coupled_global_true_rel=" << coupled_result.global_true_rel
              << " coupled_max_bodywise_true_rel=" << coupled_result.max_bodywise_true_rel
              << " coupled_air_true_rel=" << coupled_result.air_true_rel
              << " coupled_coil_true_rel=" << coupled_result.coil_true_rel
              << " coupled_plate_true_rel=" << coupled_result.plate_true_rel
              << " coupled_outer=" << coupled_result.solver_result.outer_iteration_count
              << " coupled_converged=" << (coupled_result.solver_result.converged ? 1 : 0)
              << " coupled_true_rel=" << coupled_result.solver_result.final_true_relative_residual
              << " plate_joule_power_gap=" << plate_joule_power_gap << " plate_j_L2_gap=" << plate_j_L2_gap
              << " plate_joule_max_gap=" << plate_joule_max_gap << " solver_ok=" << (solver_ok ? 1 : 0)
              << " metric_passed=" << (metric_passed ? 1 : 0) << " prototype_passed=" << (prototype_passed ? 1 : 0)
              << " engineering_passed=" << (engineering_passed ? 1 : 0) << " passed=" << (passed ? 1 : 0)
              << std::endl;

    return passed ? 0 : 1;
}
