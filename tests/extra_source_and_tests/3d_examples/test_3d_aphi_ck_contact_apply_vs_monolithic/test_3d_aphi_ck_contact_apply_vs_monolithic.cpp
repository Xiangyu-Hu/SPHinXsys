/**
 * Stage 10-contact-2: monolithic Inner debug K(X) vs split two-body Inner+Contact.
 * Also compares Joule heat source from grad(phi) Inner+Contact path.
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

class ZeroJouleGradFieldsCK : public LocalDynamics
{
  public:
    ZeroJouleGradFieldsCK(SPHBody &sph_body, const AphiJouleHeatingFieldNames &field_names)
        : LocalDynamics(sph_body),
          dv_grad_phi_real_(particles_->template getVariableByName<Vecd>(field_names.grad_phi_real)),
          dv_grad_phi_imag_(particles_->template getVariableByName<Vecd>(field_names.grad_phi_imag))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : grad_phi_real_(encloser.dv_grad_phi_real_->DelegatedData(ex_policy)),
              grad_phi_imag_(encloser.dv_grad_phi_imag_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            grad_phi_real_[index_i] = Vecd::Zero();
            grad_phi_imag_[index_i] = Vecd::Zero();
        }

      protected:
        Vecd *grad_phi_real_;
        Vecd *grad_phi_imag_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_grad_phi_real_;
    DiscreteVariable<Vecd> *dv_grad_phi_imag_;
};

struct Summary
{
    size_t lhs_matched = 0;
    size_t lhs_missing = 0;
    Real lhs_max_abs_diff = 0.0;
    size_t joule_matched = 0;
    size_t joule_missing = 0;
    Real joule_max_abs_diff = 0.0;
};

inline void runMonolithicApply(AphiBlockMapByPosition &lhs_map, std::map<AphiContactPositionKey, Real> &joule_map,
                               int ac, char *av[], Real dp_0, Real body_length, Real body_height, Real body_width,
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
    AphiJouleHeatingFieldNames joule_names;
    AphiLhsAssemblyOptions options;
    options.omega = omega;

    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_aphi(body, sigma_left, nu, names);
    StateDynamics<MainExecutionPolicy, SetAphiMaterialPropertiesCK> set_material(body, sigma_left, nu, names.material);
    RegisterAphiJouleHeatingFieldsCK register_joule(body, joule_names);
    StateDynamics<MainExecutionPolicy, AssignPiecewiseSigmaHalfSpaceCK> assign_sigma(
        body, x_interface, sigma_left, sigma_right, nu, names.material);
    StateDynamics<MainExecutionPolicy, AssignSeparableAphiFieldsCK> assign_fields(body, names);

    UpdateCellLinkedList<MainExecutionPolicy, RealBody> update_cell_linked_list(body);
    UpdateRelation<MainExecutionPolicy, Inner<>> update_inner_relation(inner_ck);

    AphiAssembleLhsDebugDynamicsBundle<MainExecutionPolicy> assemble_debug(body, inner_ck, names, options);

    InteractionDynamicsCK<MainExecutionPolicy, AphiComputeScalarPhiGradientCK<Inner<>>> compute_grad_phi(
        inner_ck, names.solution, joule_names);
    StateDynamics<MainExecutionPolicy, AphiComputeFrequencyElectricFieldCK> compute_electric_field(
        body, omega, names.solution, joule_names);
    StateDynamics<MainExecutionPolicy, AphiComputeJouleHeatSourceCK> compute_joule(body, names.material, joule_names);

    initialize_aphi.exec();
    set_material.exec();
    assign_sigma.exec();
    assign_fields.exec();
    update_cell_linked_list.exec();
    update_inner_relation.exec();

    assemble_debug.exec();
    compute_grad_phi.exec();
    compute_electric_field.exec();
    compute_joule.exec();

    BaseParticles &particles = body.getBaseParticles();
    collectCoreBlockByPosition(lhs_map, particles, names.lhs, body_length, body_height, body_width, core_shell,
                               x_interface, dp_0);
    collectCoreScalarByPosition(joule_map, particles, joule_names.joule_heat_source, body_length, body_height,
                                body_width, core_shell, x_interface, dp_0);
}

inline void runSplitContactApply(AphiBlockMapByPosition &lhs_map, std::map<AphiContactPositionKey, Real> &joule_map,
                                 int ac, char *av[], Real dp_0, Real body_length, Real body_height, Real body_width,
                                 Real boundary_width, Real core_shell, Real x_interface, Real sigma_left, Real sigma_right,
                                 Real nu, Real omega)
{
    (void)sigma_right;
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
    AphiJouleHeatingFieldNames joule_names;
    AphiLhsAssemblyOptions options;
    options.omega = omega;

    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_left(left_body, sigma_left, nu, names);
    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_right(right_body, sigma_right, nu, names);
    StateDynamics<MainExecutionPolicy, SetAphiMaterialPropertiesCK> set_left_material(left_body, sigma_left, nu, names.material);
    StateDynamics<MainExecutionPolicy, SetAphiMaterialPropertiesCK> set_right_material(right_body, sigma_right, nu, names.material);
    RegisterAphiJouleHeatingFieldsCK register_left_joule(left_body, joule_names);
    RegisterAphiJouleHeatingFieldsCK register_right_joule(right_body, joule_names);
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

    AphiAssembleLhsDebugContactDynamicsBundle<MainExecutionPolicy> assemble_left(left_body, left_inner, left_to_right,
                                                                                 names, options);
    AphiAssembleLhsDebugContactDynamicsBundle<MainExecutionPolicy> assemble_right(right_body, right_inner, right_to_left,
                                                                                  names, options);

    InteractionDynamicsCK<MainExecutionPolicy, AphiComputeScalarPhiGradientCK<Inner<>, Contact<>>> left_grad_phi(
        DynamicsArgs(left_inner, names.solution, joule_names),
        DynamicsArgs(left_to_right, names.solution, joule_names));
    InteractionDynamicsCK<MainExecutionPolicy, AphiComputeScalarPhiGradientCK<Inner<>, Contact<>>> right_grad_phi(
        DynamicsArgs(right_inner, names.solution, joule_names),
        DynamicsArgs(right_to_left, names.solution, joule_names));

    StateDynamics<MainExecutionPolicy, ZeroJouleGradFieldsCK> zero_left_grad(left_body, joule_names);
    StateDynamics<MainExecutionPolicy, ZeroJouleGradFieldsCK> zero_right_grad(right_body, joule_names);
    StateDynamics<MainExecutionPolicy, AphiComputeFrequencyElectricFieldCK> left_electric_field(
        left_body, omega, names.solution, joule_names);
    StateDynamics<MainExecutionPolicy, AphiComputeFrequencyElectricFieldCK> right_electric_field(
        right_body, omega, names.solution, joule_names);
    StateDynamics<MainExecutionPolicy, AphiComputeJouleHeatSourceCK> left_joule(left_body, names.material, joule_names);
    StateDynamics<MainExecutionPolicy, AphiComputeJouleHeatSourceCK> right_joule(right_body, names.material, joule_names);

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

    assemble_left.exec();
    assemble_right.exec();
    zero_left_grad.exec();
    zero_right_grad.exec();
    left_grad_phi.exec();
    right_grad_phi.exec();
    left_electric_field.exec();
    right_electric_field.exec();
    left_joule.exec();
    right_joule.exec();

    for (auto *body_ptr : {&left_body, &right_body})
    {
        BaseParticles &particles = body_ptr->getBaseParticles();
        collectCoreBlockByPosition(lhs_map, particles, names.lhs, body_length, body_height, body_width, core_shell,
                                   x_interface, dp_0);
        collectCoreScalarByPosition(joule_map, particles, joule_names.joule_heat_source, body_length, body_height,
                                    body_width, core_shell, x_interface, dp_0);
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
    const Real lhs_threshold = 5.0e-4;
    const Real joule_threshold = 5.0e-4;

    AphiBlockMapByPosition monolithic_lhs;
    AphiBlockMapByPosition contact_lhs;
    std::map<AphiContactPositionKey, Real> monolithic_joule;
    std::map<AphiContactPositionKey, Real> contact_joule;

    runMonolithicApply(monolithic_lhs, monolithic_joule, ac, av, dp_0, body_length, body_height, body_width,
                       boundary_width, core_shell, x_interface, sigma_left, sigma_right, nu, omega);
    runSplitContactApply(contact_lhs, contact_joule, ac, av, dp_0, body_length, body_height, body_width, boundary_width,
                         core_shell, x_interface, sigma_left, sigma_right, nu, omega);

    Summary summary;
    summary.lhs_max_abs_diff = maxAbsBlockDifference(monolithic_lhs, contact_lhs, summary.lhs_matched, summary.lhs_missing);
    summary.joule_max_abs_diff =
        maxAbsScalarDifference(monolithic_joule, contact_joule, summary.joule_matched, summary.joule_missing);

    const bool passed = summary.lhs_matched > 0 && summary.lhs_missing == 0 && summary.joule_matched > 0 &&
                        summary.joule_missing == 0 && summary.lhs_max_abs_diff < lhs_threshold &&
                        summary.joule_max_abs_diff < joule_threshold;

    std::cout << "test_3d_aphi_ck_contact_apply_vs_monolithic"
              << " lhs_matched=" << summary.lhs_matched
              << " lhs_missing=" << summary.lhs_missing
              << " lhs_max_abs_diff=" << summary.lhs_max_abs_diff
              << " joule_matched=" << summary.joule_matched
              << " joule_missing=" << summary.joule_missing
              << " joule_max_abs_diff=" << summary.joule_max_abs_diff
              << " passed=" << (passed ? 1 : 0) << std::endl;

    return passed ? 0 : 1;
}
