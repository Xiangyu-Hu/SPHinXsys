/**
 * Stage 10-contact-3: monolithic Inner vs split Inner+Contact block-Jacobi diagonal equivalence.
 */
#include "sphinxsys.h"
#include "electromagnetic_dynamics/all_electromagnetic_dynamics_ck.h"
#include "electromagnetic_dynamics/aphi_block_jacobi_preconditioner_ck.hpp"
#include "electromagnetic_dynamics/test_helpers/aphi_contact_test_helpers.h"

#include <iostream>
#include <map>

using namespace SPH;
using namespace SPH::electromagnetics;
using namespace SPH::electromagnetics::test;
using namespace SPH::electromagnetics::benchmark;

namespace
{

using MainExecutionPolicy = execution::MainExecutionPolicy;
using RealMap = std::map<AphiContactPositionKey, Real>;
using VecdMap = std::map<AphiContactPositionKey, Vecd>;

struct Summary
{
    size_t matched = 0;
    size_t missing = 0;
    Real max_laplace_a = 0.0;
    Real max_laplace_phi = 0.0;
    Real max_grad_phi = 0.0;
    Real max_div_a = 0.0;
};

inline void runMonolithicJacobi(RealMap &laplace_a, RealMap &laplace_phi, VecdMap &grad_phi, VecdMap &div_a, int ac,
                                char *av[], Real dp_0, Real body_length, Real body_height, Real body_width,
                                Real boundary_width, Real core_shell, Real x_interface, Real sigma_left, Real sigma_right,
                                Real nu, Real omega)
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
    const AphiBlockJacobiDiagonalNames diag_names;

    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_aphi(body, sigma_left, nu, names);
    StateDynamics<MainExecutionPolicy, SetAphiMaterialPropertiesCK> set_material(body, sigma_left, nu, names.material);
    StateDynamics<MainExecutionPolicy, AssignPiecewiseSigmaHalfSpaceCK> assign_sigma(
        body, x_interface, sigma_left, sigma_right, nu, names.material);

    UpdateCellLinkedList<MainExecutionPolicy, RealBody> update_cell_linked_list(body);
    UpdateRelation<MainExecutionPolicy, Inner<>> update_inner_relation(inner_ck);

    InteractionDynamicsCK<MainExecutionPolicy, AphiComputeBlockJacobiDiagonalCK<Inner<>>> compute_jacobi(
        DynamicsArgs(inner_ck, names.material, omega, options));

    initialize_aphi.exec();
    set_material.exec();
    assign_sigma.exec();
    update_cell_linked_list.exec();
    update_inner_relation.exec();
    compute_jacobi.exec();

    BaseParticles &particles = body.getBaseParticles();
    collectCoreRealByPosition(laplace_a, particles, diag_names.laplace_a_diag, body_length, body_height, body_width,
                              core_shell, x_interface, dp_0);
    collectCoreRealByPosition(laplace_phi, particles, diag_names.laplace_phi_diag, body_length, body_height, body_width,
                              core_shell, x_interface, dp_0);
    collectCoreVecdByPosition(grad_phi, particles, diag_names.grad_phi_coupling, body_length, body_height, body_width,
                              core_shell, x_interface, dp_0);
    collectCoreVecdByPosition(div_a, particles, diag_names.div_a_coupling, body_length, body_height, body_width,
                              core_shell, x_interface, dp_0);
}

inline void runSplitContactJacobi(RealMap &laplace_a, RealMap &laplace_phi, VecdMap &grad_phi, VecdMap &div_a, int ac,
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
    const AphiBlockJacobiDiagonalNames diag_names;

    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_left(left_body, sigma_left, nu, names);
    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_right(right_body, sigma_right, nu, names);
    StateDynamics<MainExecutionPolicy, SetAphiMaterialPropertiesCK> set_left_material(left_body, sigma_left, nu, names.material);
    StateDynamics<MainExecutionPolicy, SetAphiMaterialPropertiesCK> set_right_material(right_body, sigma_right, nu, names.material);
    StateDynamics<MainExecutionPolicy, AssignConstantMaterialSigmaCK> assign_left_sigma(left_body, sigma_left, names.material);
    StateDynamics<MainExecutionPolicy, AssignConstantMaterialSigmaCK> assign_right_sigma(right_body, sigma_right, names.material);

    UpdateCellLinkedList<MainExecutionPolicy, RealBody> update_left_cell_linked_list(left_body);
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> update_right_cell_linked_list(right_body);
    UpdateRelation<MainExecutionPolicy, Inner<>> update_left_inner(left_inner);
    UpdateRelation<MainExecutionPolicy, Inner<>> update_right_inner(right_inner);
    UpdateRelation<MainExecutionPolicy, Contact<>> update_left_contact(left_to_right);
    UpdateRelation<MainExecutionPolicy, Contact<>> update_right_contact(right_to_left);

    AphiComputeBlockJacobiContactDynamicsBundle<MainExecutionPolicy> compute_left_jacobi(
        left_body, left_inner, left_to_right, names.material, omega, options);
    AphiComputeBlockJacobiContactDynamicsBundle<MainExecutionPolicy> compute_right_jacobi(
        right_body, right_inner, right_to_left, names.material, omega, options);

    initialize_left.exec();
    initialize_right.exec();
    set_left_material.exec();
    set_right_material.exec();
    assign_left_sigma.exec();
    assign_right_sigma.exec();
    update_left_cell_linked_list.exec();
    update_right_cell_linked_list.exec();
    update_left_inner.exec();
    update_right_inner.exec();
    update_left_contact.exec();
    update_right_contact.exec();
    compute_left_jacobi.exec();
    compute_right_jacobi.exec();

    for (auto *body_ptr : {&left_body, &right_body})
    {
        BaseParticles &particles = body_ptr->getBaseParticles();
        collectCoreRealByPosition(laplace_a, particles, diag_names.laplace_a_diag, body_length, body_height, body_width,
                                  core_shell, x_interface, dp_0);
        collectCoreRealByPosition(laplace_phi, particles, diag_names.laplace_phi_diag, body_length, body_height,
                                  body_width, core_shell, x_interface, dp_0);
        collectCoreVecdByPosition(grad_phi, particles, diag_names.grad_phi_coupling, body_length, body_height, body_width,
                                  core_shell, x_interface, dp_0);
        collectCoreVecdByPosition(div_a, particles, diag_names.div_a_coupling, body_length, body_height, body_width,
                                  core_shell, x_interface, dp_0);
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
    const Real validation_threshold_laplace = 8.0e-4;
    const Real validation_threshold_coupling = 5.0e-4;

    RealMap mono_laplace_a, contact_laplace_a;
    RealMap mono_laplace_phi, contact_laplace_phi;
    VecdMap mono_grad_phi, contact_grad_phi;
    VecdMap mono_div_a, contact_div_a;

    runMonolithicJacobi(mono_laplace_a, mono_laplace_phi, mono_grad_phi, mono_div_a, ac, av, dp_0, body_length,
                        body_height, body_width, boundary_width, core_shell, x_interface, sigma_left, sigma_right, nu,
                        omega);
    runSplitContactJacobi(contact_laplace_a, contact_laplace_phi, contact_grad_phi, contact_div_a, ac, av, dp_0,
                          body_length, body_height, body_width, boundary_width, core_shell, x_interface, sigma_left,
                          sigma_right, nu, omega);

    Summary summary;
    summary.max_laplace_a =
        maxAbsScalarDifference(mono_laplace_a, contact_laplace_a, summary.matched, summary.missing);
    size_t matched_phi = 0;
    size_t missing_phi = 0;
    summary.max_laplace_phi =
        maxAbsScalarDifference(mono_laplace_phi, contact_laplace_phi, matched_phi, missing_phi);
    size_t matched_grad = 0;
    size_t missing_grad = 0;
    summary.max_grad_phi = maxAbsVecdMapDifference(mono_grad_phi, contact_grad_phi, matched_grad, missing_grad);
    size_t matched_div = 0;
    size_t missing_div = 0;
    summary.max_div_a = maxAbsVecdMapDifference(mono_div_a, contact_div_a, matched_div, missing_div);

    const bool passed = summary.matched > 0 && summary.missing == 0 && matched_phi == summary.matched &&
                        missing_phi == 0 && matched_grad == summary.matched && missing_grad == 0 &&
                        matched_div == summary.matched && missing_div == 0 &&
                        summary.max_laplace_a < validation_threshold_laplace &&
                        summary.max_laplace_phi < validation_threshold_laplace &&
                        summary.max_grad_phi < validation_threshold_coupling &&
                        summary.max_div_a < validation_threshold_coupling;

    std::cout << "test_3d_aphi_ck_contact_block_jacobi_diagonal_equivalence"
              << " matched=" << summary.matched
              << " max_laplace_a=" << summary.max_laplace_a
              << " max_laplace_phi=" << summary.max_laplace_phi
              << " max_grad_phi=" << summary.max_grad_phi
              << " max_div_a=" << summary.max_div_a
              << " passed=" << (passed ? 1 : 0) << std::endl;

    return passed ? 0 : 1;
}
