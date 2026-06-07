/**
 * Sprint 5 Step 5.1: monolithic TEAM7-like Inner vs three-body Contact apply + Joule equivalence.
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

inline AphiTeam7RegionalApplyMetrics runMonolithicTeam7Apply(int ac, char *av[], Real dp_0, Real body_length,
                                                               Real body_height, Real body_width, Real boundary_width,
                                                               Real core_shell, const AphiTeam7LikeUnitBoxLayout &layout,
                                                               Real omega)
{
    AphiLhsTestBody test_body(dp_0, body_length, body_height, body_width, boundary_width, ac, av);
    AphiVariableNames names;
    AphiJouleHeatingFieldNames joule_names;
    AphiLhsAssemblyOptions options;
    options.omega = omega;

    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_aphi(test_body.body, layout.air.sigma,
                                                                                  layout.air.nu, names);
    StateDynamics<MainExecutionPolicy, AssignTeam7LikeRegionMaterialsCK> assign_material(test_body.body, layout,
                                                                                         names.material);
    RegisterAphiJouleHeatingFieldsCK register_joule(test_body.body, joule_names);
    StateDynamics<MainExecutionPolicy, AssignSeparableAphiFieldsCK> assign_fields(test_body.body, names);

    AphiAssembleLhsDebugDynamicsBundle<MainExecutionPolicy> assemble_debug(test_body.body, test_body.inner(), names,
                                                                           options);
    InteractionDynamicsCK<MainExecutionPolicy, AphiComputeScalarPhiGradientCK<Inner<>>> compute_grad_phi(
        test_body.inner(), names.solution, joule_names);
    StateDynamics<MainExecutionPolicy, AphiComputeFrequencyElectricFieldCK> compute_electric_field(
        test_body.body, omega, names.solution, joule_names);
    StateDynamics<MainExecutionPolicy, AphiComputeJouleHeatSourceCK> compute_joule(test_body.body, names.material,
                                                                                     joule_names);

    initialize_aphi.exec();
    assign_material.exec();
    assign_fields.exec();
    test_body.updateRelations();
    assemble_debug.exec();
    compute_grad_phi.exec();
    compute_electric_field.exec();
    compute_joule.exec();

    return hostTeam7RegionalApplyMetricsFromMonolithic(test_body.body, names, joule_names, layout, body_length,
                                                       body_height, body_width, core_shell, dp_0);
}

inline AphiTeam7RegionalApplyMetrics runSplitTeam7ContactApply(int ac, char *av[], Real dp_0, Real body_length,
                                                                 Real body_height, Real body_width, Real boundary_width,
                                                                 Real core_shell, Real omega)
{
    AphiTeam7ThreeBodyContactCase case_setup(dp_0, body_length, body_height, body_width, boundary_width, ac, av);
    AphiVariableNames names;
    AphiJouleHeatingFieldNames joule_names;
    AphiLhsAssemblyOptions options;
    options.omega = omega;

    initializeTeam7ThreeBodyFields(case_setup, names);
    RegisterAphiJouleHeatingFieldsCK register_air_joule(case_setup.air_body, joule_names);
    RegisterAphiJouleHeatingFieldsCK register_coil_joule(case_setup.coil_body, joule_names);
    RegisterAphiJouleHeatingFieldsCK register_plate_joule(case_setup.plate_body, joule_names);
    StateDynamics<MainExecutionPolicy, AssignSeparableAphiFieldsCK> assign_air_fields(case_setup.air_body, names);
    StateDynamics<MainExecutionPolicy, AssignSeparableAphiFieldsCK> assign_coil_fields(case_setup.coil_body, names);
    StateDynamics<MainExecutionPolicy, AssignSeparableAphiFieldsCK> assign_plate_fields(case_setup.plate_body, names);
    (void)register_air_joule;
    (void)register_coil_joule;
    (void)register_plate_joule;

    AphiAssembleLhsDebugContactDynamicsBundle<MainExecutionPolicy> assemble_air(
        case_setup.air_body, case_setup.air_inner(), case_setup.air_contact(), names, options);
    AphiAssembleLhsDebugContactDynamicsBundle<MainExecutionPolicy> assemble_coil(
        case_setup.coil_body, case_setup.coil_inner(), case_setup.coil_to_air(), names, options);
    AphiAssembleLhsDebugContactDynamicsBundle<MainExecutionPolicy> assemble_plate(
        case_setup.plate_body, case_setup.plate_inner(), case_setup.plate_to_air(), names, options);

    assign_air_fields.exec();
    assign_coil_fields.exec();
    assign_plate_fields.exec();
    case_setup.updateRelations();
    assemble_air.exec();
    assemble_coil.exec();
    assemble_plate.exec();
    runTeam7ContactJoulePipeline(case_setup, names, joule_names, omega);

    return hostTeam7RegionalApplyMetricsFromContactCase(case_setup, names, joule_names, body_length, body_height,
                                                        body_width, core_shell, dp_0);
}

} // namespace

int main(int ac, char *av[])
{
    const Real dp_0 = 0.1;
    const Real body_length = 1.0;
    const Real body_height = 1.0;
    const Real body_width = 1.0;
    const Real boundary_width = 3.0 * dp_0;
    const Real core_shell = 2.5 * dp_0;
    const Real omega = 1.25;
    const AphiTeam7LikeUnitBoxLayout layout = buildTeam7LayoutForBox(body_length, body_height, body_width);
    const Real max_regional_lhs_gap = 0.05;
    const Real max_regional_joule_gap = 0.10;

    const AphiTeam7RegionalApplyMetrics mono = runMonolithicTeam7Apply(ac, av, dp_0, body_length, body_height, body_width,
                                                                       boundary_width, core_shell, layout, omega);
    const AphiTeam7RegionalApplyMetrics contact =
        runSplitTeam7ContactApply(ac, av, dp_0, body_length, body_height, body_width, boundary_width, core_shell, omega);

    const Real conductor_lhs_gap = relativeMetricChange(mono.mono_conductor_lhs_max, contact.contact_conductor_lhs_max);
    const Real coil_lhs_gap = relativeMetricChange(mono.mono_coil_lhs_max, contact.contact_coil_lhs_max);
    const Real air_lhs_gap = relativeMetricChange(mono.mono_air_lhs_max, contact.contact_air_lhs_max);
    const Real conductor_joule_gap =
        relativeMetricChange(mono.mono_conductor_joule_max, contact.contact_conductor_joule_max);

    const bool passed = mono.mono_conductor_lhs_max > TinyReal && mono.mono_coil_lhs_max > TinyReal &&
                        mono.mono_air_lhs_max > TinyReal && conductor_lhs_gap < max_regional_lhs_gap &&
                        coil_lhs_gap < max_regional_lhs_gap && air_lhs_gap < max_regional_lhs_gap &&
                        conductor_joule_gap < max_regional_joule_gap;

    std::cout << "test_3d_aphi_ck_contact_team7_three_body_apply_vs_monolithic"
              << " mono_conductor_lhs_max=" << mono.mono_conductor_lhs_max
              << " contact_conductor_lhs_max=" << contact.contact_conductor_lhs_max
              << " mono_coil_lhs_max=" << mono.mono_coil_lhs_max << " contact_coil_lhs_max=" << contact.contact_coil_lhs_max
              << " mono_air_lhs_max=" << mono.mono_air_lhs_max << " contact_air_lhs_max=" << contact.contact_air_lhs_max
              << " mono_conductor_joule_max=" << mono.mono_conductor_joule_max
              << " contact_conductor_joule_max=" << contact.contact_conductor_joule_max
              << " conductor_lhs_gap=" << conductor_lhs_gap << " coil_lhs_gap=" << coil_lhs_gap
              << " air_lhs_gap=" << air_lhs_gap << " conductor_joule_gap=" << conductor_joule_gap
              << " passed=" << (passed ? 1 : 0) << std::endl;

    return passed ? 0 : 1;
}
