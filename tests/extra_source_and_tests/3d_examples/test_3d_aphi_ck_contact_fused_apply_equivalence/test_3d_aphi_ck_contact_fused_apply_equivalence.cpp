/**
 * Stage 10-contact-3: fused Inner+Contact apply vs debug assemble; monolithic vs split equivalence.
 */
#include "sphinxsys.h"
#include "electromagnetic_dynamics/all_electromagnetic_dynamics_ck.h"
#include "electromagnetic_dynamics/diagnostics/aphi_assemble_lhs_debug_ck.hpp"
#include "electromagnetic_dynamics/test_helpers/aphi_contact_test_helpers.h"

#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics;
using namespace SPH::electromagnetics::test;
using namespace SPH::electromagnetics::benchmark;

namespace
{

using MainExecutionPolicy = execution::MainExecutionPolicy;

struct Summary
{
    size_t mono_vs_split_matched = 0;
    Real mono_vs_split_max_diff = 0.0;
    Real split_debug_vs_fused_max_diff = 0.0;
};

inline void runMonolithicDebug(AphiBlockMapByPosition &lhs_map, int ac, char *av[], Real dp_0, Real body_length,
                               Real body_height, Real body_width, Real boundary_width, Real core_shell,
                               Real x_interface, Real sigma_left, Real sigma_right, Real nu, Real omega)
{
    BoundingBoxd system_bounds(Vecd(-boundary_width, -boundary_width, -boundary_width),
                               Vecd(body_length + boundary_width, body_height + boundary_width,
                                    body_width + boundary_width));

    SPHSystem sph_system(system_bounds, dp_0);
    sph_system.handleCommandlineOptions(ac, av);

    const Vecd center(0.5 * body_length, 0.5 * body_height, 0.5 * body_width);
    const Vecd halfsize(0.5 * body_length, 0.5 * body_height, 0.5 * body_width);
    SolidBody body(sph_system, makeShared<AphiHalfSpaceBoxShape>("MonolithicBody", center, halfsize));
    body.defineAdaptation<SPHAdaptation>(1.15, 1.0);
    body.defineMaterial<Solid>();
    body.defineBodyLevelSetShape();
    body.generateParticles<BaseParticles, Lattice>();

    Inner<> inner_ck(body);
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();

    AphiVariableNames names;
    AphiLhsAssemblyOptions options;
    options.omega = omega;

    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_aphi(body, sigma_left, nu, names);
    StateDynamics<MainExecutionPolicy, SetAphiMaterialPropertiesCK> set_material(body, sigma_left, nu, names.material);
    StateDynamics<MainExecutionPolicy, AssignPiecewiseSigmaHalfSpaceCK> assign_sigma(
        body, x_interface, sigma_left, sigma_right, nu, names.material);
    StateDynamics<MainExecutionPolicy, AssignSeparableAphiFieldsCK> assign_fields(body, names);

    UpdateCellLinkedList<MainExecutionPolicy, RealBody> update_cell_linked_list(body);
    UpdateRelation<MainExecutionPolicy, Inner<>> update_inner_relation(inner_ck);

    AphiAssembleLhsDebugDynamicsBundle<MainExecutionPolicy> assemble_debug(body, inner_ck, names, options);

    initialize_aphi.exec();
    set_material.exec();
    assign_sigma.exec();
    assign_fields.exec();
    update_cell_linked_list.exec();
    update_inner_relation.exec();
    assemble_debug.exec();

    collectCoreBlockByPosition(lhs_map, body.getBaseParticles(), names.lhs, body_length, body_height, body_width,
                               core_shell, x_interface, dp_0);
}

inline void runSplitContactFusedAndDebug(AphiBlockMapByPosition &fused_map, AphiBlockMapByPosition &debug_map, int ac,
                                         char *av[], Real dp_0, Real body_length, Real body_height, Real body_width,
                                         Real boundary_width, Real core_shell, Real x_interface, Real sigma_left,
                                         Real sigma_right, Real nu, Real omega)
{
    BoundingBoxd system_bounds(Vecd(-boundary_width, -boundary_width, -boundary_width),
                               Vecd(body_length + boundary_width, body_height + boundary_width,
                                    body_width + boundary_width));

    SPHSystem sph_system(system_bounds, dp_0);
    sph_system.handleCommandlineOptions(ac, av);

    const Vecd left_center(0.25 * body_length, 0.5 * body_height, 0.5 * body_width);
    const Vecd right_center(0.75 * body_length, 0.5 * body_height, 0.5 * body_width);
    const Vecd halfsize(0.25 * body_length, 0.5 * body_height, 0.5 * body_width);

    SolidBody left_body(sph_system, makeShared<AphiHalfSpaceBoxShape>("LeftBody", left_center, halfsize));
    SolidBody right_body(sph_system, makeShared<AphiHalfSpaceBoxShape>("RightBody", right_center, halfsize));
    for (auto *body_ptr : {&left_body, &right_body})
    {
        body_ptr->defineAdaptation<SPHAdaptation>(1.15, 1.0);
        body_ptr->defineMaterial<Solid>();
        body_ptr->defineBodyLevelSetShape();
        body_ptr->generateParticles<BaseParticles, Lattice>();
    }

    Inner<> left_inner(left_body);
    Inner<> right_inner(right_body);
    Contact<> left_to_right(left_body, {&right_body});
    Contact<> right_to_left(right_body, {&left_body});

    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();

    AphiVariableNames names;
    AphiLhsAssemblyOptions options;
    options.omega = omega;

    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_left(left_body, sigma_left, nu, names);
    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_right(right_body, sigma_right, nu, names);
    StateDynamics<MainExecutionPolicy, SetAphiMaterialPropertiesCK> set_left_material(left_body, sigma_left, nu, names.material);
    StateDynamics<MainExecutionPolicy, SetAphiMaterialPropertiesCK> set_right_material(right_body, sigma_right, nu, names.material);
    StateDynamics<MainExecutionPolicy, AssignConstantMaterialSigmaCK> assign_left_sigma(left_body, sigma_left, names.material);
    StateDynamics<MainExecutionPolicy, AssignConstantMaterialSigmaCK> assign_right_sigma(right_body, sigma_right, names.material);
    StateDynamics<MainExecutionPolicy, AssignSeparableAphiFieldsCK> assign_left_fields(left_body, names);
    StateDynamics<MainExecutionPolicy, AssignSeparableAphiFieldsCK> assign_right_fields(right_body, names);

    UpdateCellLinkedList<MainExecutionPolicy, RealBody> update_left_cell_linked_list(left_body);
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> update_right_cell_linked_list(right_body);
    UpdateRelation<MainExecutionPolicy, Inner<>> update_left_inner(left_inner);
    UpdateRelation<MainExecutionPolicy, Inner<>> update_right_inner(right_inner);
    UpdateRelation<MainExecutionPolicy, Contact<>> update_left_contact(left_to_right);
    UpdateRelation<MainExecutionPolicy, Contact<>> update_right_contact(right_to_left);

    AphiAssembleLhsDebugContactDynamicsBundle<MainExecutionPolicy> assemble_left_debug(
        left_body, left_inner, left_to_right, names, options);
    AphiAssembleLhsDebugContactDynamicsBundle<MainExecutionPolicy> assemble_right_debug(
        right_body, right_inner, right_to_left, names, options);
    AphiApplyContactDynamicsBundle<MainExecutionPolicy> fused_left(
        left_body, left_inner, left_to_right, names.solution, names.v, names.material, omega, options);
    AphiApplyContactDynamicsBundle<MainExecutionPolicy> fused_right(
        right_body, right_inner, right_to_left, names.solution, names.v, names.material, omega, options);

    initialize_left.exec();
    initialize_right.exec();
    set_left_material.exec();
    set_right_material.exec();
    assign_left_sigma.exec();
    assign_right_sigma.exec();
    assign_left_fields.exec();
    assign_right_fields.exec();
    update_left_cell_linked_list.exec();
    update_right_cell_linked_list.exec();
    update_left_inner.exec();
    update_right_inner.exec();
    update_left_contact.exec();
    update_right_contact.exec();

    assemble_left_debug.exec();
    assemble_right_debug.exec();
    fused_left.exec();
    fused_right.exec();

    for (auto *body_ptr : {&left_body, &right_body})
    {
        BaseParticles &particles = body_ptr->getBaseParticles();
        collectCoreBlockByPosition(fused_map, particles, names.v, body_length, body_height, body_width, core_shell,
                                   x_interface, dp_0);
        collectCoreBlockByPosition(debug_map, particles, names.lhs, body_length, body_height, body_width, core_shell,
                                   x_interface, dp_0);
    }
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
    const Real x_interface = 0.5;
    const Real sigma_left = 10.0;
    const Real sigma_right = 1.0;
    const Real nu = 1.5;
    const Real omega = 1.25;
    const Real mono_vs_split_threshold = 5.0e-4;
    const Real debug_vs_fused_threshold = 1.0e-3;

    AphiBlockMapByPosition monolithic_lhs;
    AphiBlockMapByPosition split_fused;
    AphiBlockMapByPosition split_debug;

    runMonolithicDebug(monolithic_lhs, ac, av, dp_0, body_length, body_height, body_width, boundary_width, core_shell,
                       x_interface, sigma_left, sigma_right, nu, omega);
    runSplitContactFusedAndDebug(split_fused, split_debug, ac, av, dp_0, body_length, body_height, body_width,
                                   boundary_width, core_shell, x_interface, sigma_left, sigma_right, nu, omega);

    Summary summary;
    size_t missing = 0;
    summary.mono_vs_split_max_diff =
        maxAbsBlockDifference(monolithic_lhs, split_fused, summary.mono_vs_split_matched, missing);
    size_t split_matched = 0;
    size_t split_missing = 0;
    summary.split_debug_vs_fused_max_diff =
        maxAbsBlockDifference(split_debug, split_fused, split_matched, split_missing);

    const bool passed = summary.mono_vs_split_matched > 0 && missing == 0 && split_matched > 0 && split_missing == 0 &&
                        summary.mono_vs_split_max_diff < mono_vs_split_threshold &&
                        summary.split_debug_vs_fused_max_diff < debug_vs_fused_threshold;

    std::cout << "test_3d_aphi_ck_contact_fused_apply_equivalence"
              << " mono_vs_split_matched=" << summary.mono_vs_split_matched
              << " mono_vs_split_max_diff=" << summary.mono_vs_split_max_diff
              << " split_debug_vs_fused_max_diff=" << summary.split_debug_vs_fused_max_diff
              << " passed=" << (passed ? 1 : 0) << std::endl;

    return passed ? 0 : 1;
}
