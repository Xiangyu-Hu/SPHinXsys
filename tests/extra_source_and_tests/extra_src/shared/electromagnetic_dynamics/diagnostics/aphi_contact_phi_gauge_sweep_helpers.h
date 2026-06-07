#ifndef APHI_CONTACT_PHI_GAUGE_SWEEP_HELPERS_H
#define APHI_CONTACT_PHI_GAUGE_SWEEP_HELPERS_H

#include "electromagnetic_dynamics/test_helpers/aphi_contact_gmres_test_helpers.h"
#include "electromagnetic_dynamics/diagnostics/aphi_contact_left_field_error_helpers.h"
#include "electromagnetic_dynamics/diagnostics/aphi_div_a_diagnostic_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_gmres_benchmark_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_gmres_test_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_team7_contact_test_helpers.h"

namespace SPH
{
namespace electromagnetics
{
namespace test
{

inline AphiLhsAssemblyOptions fullInterfaceMmsLhsOptions(Real omega, Real phi_gauge_penalty,
                                                         bool use_phi_gauge_penalty = true)
{
    AphiLhsAssemblyOptions options;
    options.omega = omega;
    options.use_phi_gauge_penalty = use_phi_gauge_penalty;
    options.phi_gauge_penalty = phi_gauge_penalty;
    return options;
}

inline AphiLhsAssemblyOptions phiOnlyLaplacePenaltyLhsOptions(Real omega, Real phi_gauge_penalty,
                                                              bool use_phi_gauge_penalty = true)
{
    AphiLhsAssemblyOptions options;
    options.omega = omega;
    options.terms.laplace_a = false;
    options.terms.laplace_phi = true;
    options.terms.reaction = false;
    options.terms.grad_phi_coupling = false;
    options.terms.div_sigma_a_coupling = false;
    options.use_phi_gauge_penalty = use_phi_gauge_penalty;
    options.phi_gauge_penalty = phi_gauge_penalty;
    return options;
}

struct AphiContactPhiGaugeSweepRow
{
    Real phi_gauge_penalty = 0.0;
    bool use_phi_gauge_penalty = true;
    bool phi_only_operator = false;
    bool converged = false;
    UnsignedInt outer_iterations = 0;
    Real global_true_rel = 0.0;
    Real left_true_rel = 0.0;
    Real left_continuous = 0.0;
    AphiLeftBodyGaugeErrorMetrics full_left{};
    AphiLeftBodyPhiComponentErrorMetrics phi_component{};
    AphiLeftBodyGaugeErrorMetrics core_interior{};
    AphiLeftBodyPhiComponentErrorMetrics core_interior_phi{};
    AphiDivAReductionMetrics global_div_a{};
    AphiDivAReductionMetrics left_div_a{};
    AphiDivAReductionMetrics right_div_a{};
};

inline AphiContactPhiGaugeSweepRow runContactPhiGaugeSweepRow(
    int ac, char *av[], Real dp_0, Real body_length, Real body_height, Real body_width, Real boundary_width,
    Real core_shell, Real sigma_conductor, Real sigma_air, Real nu, Real omega, Real tolerance,
    Real x_interface, Real interface_band_half_width, const AphiLhsAssemblyOptions &options,
    const AphiMatrixFreeSolverOptions &solver_options, bool phi_only_operator)
{
    AphiContactPhiGaugeSweepRow row;
    row.phi_gauge_penalty = options.phi_gauge_penalty;
    row.use_phi_gauge_penalty = options.use_phi_gauge_penalty;
    row.phi_only_operator = phi_only_operator;

    AphiVariableNames names;
    AphiTwoBodyInterfaceCase contact_case(dp_0, body_length, body_height, body_width, boundary_width, 0, av);
    setupTwoBodyInterfaceMmsFields(contact_case, names, sigma_conductor, sigma_air, nu, x_interface);
    contact_case.updateRelations();

    const AphiInterfaceMmsRunMetrics metrics = runTwoBodyCoupledContactInterfaceMms(
        contact_case, names, options, solver_options, body_length, body_height, body_width, core_shell, x_interface,
        interface_band_half_width, 0);

    BaseParticles &left_particles = contact_case.left_body.getBaseParticles();
    syncVariableToHost<Vecd>(left_particles, "Position");
    const size_t left_total_real_particles = left_particles.TotalRealParticles();
    const Vecd *left_positions = left_particles.getVariableDataByName<Vecd>("Position");

    row.converged = metrics.solver_result.converged;
    row.outer_iterations = metrics.solver_result.outer_iteration_count;
    row.global_true_rel = metrics.global_true_rel;
    row.left_true_rel = metrics.left_true_rel;
    row.left_continuous = metrics.left_continuous_error;
    row.full_left = hostLeftBodyGaugeErrorMetrics(left_particles, names.solution, names.r_hat, left_positions,
                                                  left_total_real_particles, body_length, body_height, body_width,
                                                  core_shell, x_interface);
    row.phi_component = hostLeftBodyPhiComponentErrorMetrics(left_particles, names.solution, names.r_hat, left_positions,
                                                             left_total_real_particles, body_length, body_height,
                                                             body_width, core_shell, x_interface);
    row.core_interior =
        hostLeftBodyGaugeErrorMetrics(left_particles, names.solution, names.r_hat, left_positions,
                                      left_total_real_particles, body_length, body_height, body_width, core_shell,
                                      x_interface, AphiLeftBodyErrorZone::CoreInterior, interface_band_half_width);
    row.core_interior_phi = hostLeftBodyPhiComponentErrorMetrics(
        left_particles, names.solution, names.r_hat, left_positions, left_total_real_particles, body_length,
        body_height, body_width, core_shell, x_interface, AphiLeftBodyErrorZone::CoreInterior,
        interface_band_half_width);

    row.global_div_a = hostTwoBodyContactDivAReductionMetrics(contact_case, names);
    row.left_div_a = hostTwoBodyContactDivAReductionMetrics(
        contact_case, names, [&](const Vecd &position) {
            return isLeftHalfCoreParticle(position, body_length, body_height, body_width, core_shell, x_interface);
        });
    row.right_div_a = hostTwoBodyContactDivAReductionMetrics(
        contact_case, names, [&](const Vecd &position) {
            return position[0] > x_interface + TinyReal &&
                   isCoreParticle(position, body_length, body_height, body_width, core_shell);
        });
    return row;
}

struct AphiTeam7PhiGaugeSweepRow
{
    Real phi_gauge_penalty = 0.0;
    bool converged = false;
    UnsignedInt outer_iterations = 0;
    Real global_true_rel = 0.0;
    Real max_bodywise_true_rel = 0.0;
    Real plate_joule_power_gap = 0.0;
    Real plate_j_L2_gap = 0.0;
    AphiDivAReductionMetrics global_div_a{};
    AphiDivAReductionMetrics plate_div_a{};
};

inline AphiTeam7PhiGaugeSweepRow runTeam7PhiGaugeSweepRow(
    int ac, char *av[], Real dp_0, Real body_length, Real body_height, Real body_width, Real boundary_width,
    Real core_shell, Real omega, Real phi_gauge_penalty, Real tolerance, const Vecd &coil_current_real,
    const Vecd &coil_current_imag, Real impressed_current_amplitude, const AphiTeam7PlateObservables &mono_plate)
{
    AphiTeam7PhiGaugeSweepRow row;
    row.phi_gauge_penalty = phi_gauge_penalty;

    AphiVariableNames names;
    AphiJouleHeatingFieldNames joule_names;
    AphiLhsAssemblyOptions options;
    options.omega = omega;
    options.use_phi_gauge_penalty = true;
    options.phi_gauge_penalty = phi_gauge_penalty;
    const AphiMatrixFreeSolverOptions solver_options = defaultCoupledContactGMRESOptions(tolerance, 50, 100);

    AphiTeam7ThreeBodyContactCase contact_case(dp_0, body_length, body_height, body_width, boundary_width, 0, av);
    const AphiTeam7ContactCoupledGmresResult coupled_result = runTeam7ThreeBodyCoupledContactGmres(
        contact_case, names, joule_names, options, solver_options, body_length, body_height, body_width, core_shell,
        coil_current_real, coil_current_imag, impressed_current_amplitude);

    row.converged = coupled_result.solver_result.converged;
    row.outer_iterations = coupled_result.solver_result.outer_iteration_count;
    row.global_true_rel = coupled_result.global_true_rel;
    row.max_bodywise_true_rel = coupled_result.max_bodywise_true_rel;
    row.plate_joule_power_gap =
        relativeMetricChange(mono_plate.plate_joule_power, coupled_result.plate.plate_joule_power);
    row.plate_j_L2_gap = relativeMetricChange(mono_plate.plate_j_L2, coupled_result.plate.plate_j_L2);
    row.global_div_a = hostTeam7ContactDivAReductionMetrics(contact_case, names);
    row.plate_div_a = hostTeam7ContactDivAReductionMetrics(
        contact_case, names, [&](const Vecd &position) {
            const benchmark::AphiTeam7LikeUnitBoxLayout layout =
                benchmark::buildTeam7LayoutForBox(body_length, body_height, body_width);
            return team7ParticleInRegion(position, layout, AphiBenchmarkMaterialRegion::Conductor);
        });
    return row;
}

} // namespace test
} // namespace electromagnetics
} // namespace SPH

#endif // APHI_CONTACT_PHI_GAUGE_SWEEP_HELPERS_H
