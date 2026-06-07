#ifndef APHI_CONTACT_GRADDIV_PC_DIAGNOSTIC_HELPERS_H
#define APHI_CONTACT_GRADDIV_PC_DIAGNOSTIC_HELPERS_H

#include "electromagnetic_dynamics/aphi_block_jacobi_preconditioner_ck.hpp"
#include "electromagnetic_dynamics/aphi_block_zero_ck.hpp"
#include "electromagnetic_dynamics/aphi_contact_pairwise_a_divergence_penalty_pipeline.h"
#include "electromagnetic_dynamics/aphi_pairwise_a_divergence_penalty_pipeline.h"
#include "electromagnetic_dynamics/diagnostics/aphi_graddiv_block_diagnostic_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_contact_test_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_test_device_sync.h"

namespace SPH
{
namespace electromagnetics
{
namespace test
{

struct AphiContactGradDivPcConsistencyMetrics
{
    size_t core_pc_probed = 0;
    Real core_pc_max_abs_diff = 0.0;
    size_t interface_pc_probed = 0;
    Real interface_pc_max_abs_diff = 0.0;
};

class AssignProbeUnitBasisAFieldCK : public LocalDynamics
{
  public:
    AssignProbeUnitBasisAFieldCK(SPHBody &sph_body, const std::string &a_real_name, size_t probe_index,
                                 UnsignedInt basis_index)
        : LocalDynamics(sph_body),
          dv_a_real_(particles_->template getVariableByName<Vecd>(a_real_name)), probe_index_(probe_index),
          basis_index_(basis_index)
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : a_real_(encloser.dv_a_real_->DelegatedData(ex_policy)), probe_index_(encloser.probe_index_),
              basis_index_(encloser.basis_index_)
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            if (index_i == probe_index_)
            {
                a_real_[index_i] = Vecd(basis_index_ == 0 ? 1.0 : 0.0, basis_index_ == 1 ? 1.0 : 0.0,
                                        basis_index_ == 2 ? 1.0 : 0.0);
            }
            else
            {
                a_real_[index_i] = Vecd(0.0, 0.0, 0.0);
            }
        }

      protected:
        Vecd *a_real_;
        size_t probe_index_;
        UnsignedInt basis_index_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_a_real_;
    size_t probe_index_;
    UnsignedInt basis_index_;
};

inline void registerProbePenaltyFields(BaseParticles &particles, const AphiBlockNames &probe_input,
                                       const AphiBlockNames &probe_output)
{
    particles.registerStateVariable<Vecd>(probe_input.a_real, ZeroData<Vecd>::value);
    particles.registerStateVariable<Vecd>(probe_input.a_imag, ZeroData<Vecd>::value);
    particles.registerStateVariable<Real>(probe_input.phi_real, Real(0));
    particles.registerStateVariable<Real>(probe_input.phi_imag, Real(0));
    particles.registerStateVariable<Vecd>(probe_output.a_real, ZeroData<Vecd>::value);
    particles.registerStateVariable<Vecd>(probe_output.a_imag, ZeroData<Vecd>::value);
    particles.registerStateVariable<Real>(probe_output.phi_real, Real(0));
    particles.registerStateVariable<Real>(probe_output.phi_imag, Real(0));
}

template <class PenaltyPipeline>
inline Matd probeGradDivBlockByFiniteDifference(
    SPHBody &body, PenaltyPipeline &penalty_pipeline,
    StateDynamics<execution::MainExecutionPolicy, AphiZeroBlockCK> &zero_probe_input,
    StateDynamics<execution::MainExecutionPolicy, AphiZeroBlockCK> &zero_probe_output, BaseParticles &particles,
    const AphiBlockNames &probe_input, const AphiBlockNames &probe_output, size_t probe_index)
{
    Matd fd_block = Matd::Zero();
    for (UnsignedInt basis = 0; basis < 3; ++basis)
    {
        zero_probe_input.exec();
        zero_probe_output.exec();
        StateDynamics<execution::MainExecutionPolicy, AssignProbeUnitBasisAFieldCK> assign_probe(
            body, probe_input.a_real, probe_index, basis);
        assign_probe.exec();
        penalty_pipeline.exec();
        syncVariableToHost<Vecd>(particles, probe_output.a_real);
        const Vecd *lhs_a = particles.getVariableDataByName<Vecd>(probe_output.a_real);
        const Vecd response = -lhs_a[probe_index];
        for (UnsignedInt row = 0; row < 3; ++row)
        {
            fd_block(row, basis) = response[row];
        }
    }
    return fd_block;
}

template <class PenaltyPipeline>
inline void probeContactGradDivPcOnBody(
    SPHBody &body, PenaltyPipeline &penalty_pipeline,
    StateDynamics<execution::MainExecutionPolicy, AphiZeroBlockCK> &zero_probe_input,
    StateDynamics<execution::MainExecutionPolicy, AphiZeroBlockCK> &zero_probe_output, BaseParticles &particles,
    const AphiBlockNames &probe_input, const AphiBlockNames &probe_output,
    const AphiBlockJacobiDiagonalNames &diag_names, Real body_length, Real body_height, Real body_width,
    Real core_shell, Real x_interface, Real dp_0, size_t max_probe_per_zone,
    AphiContactGradDivPcConsistencyMetrics &metrics)
{
    syncVariableToHost<Matd>(particles, diag_names.graddiv_a_block);
    syncVariableToHost<Vecd>(particles, "Position");
    const Matd *graddiv_blocks = particles.getVariableDataByName<Matd>(diag_names.graddiv_a_block);
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    const size_t total_real_particles = particles.TotalRealParticles();

    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!isCoreParticle(positions[i], body_length, body_height, body_width, core_shell))
        {
            continue;
        }

        const Matd fd_block = probeGradDivBlockByFiniteDifference(body, penalty_pipeline, zero_probe_input,
                                                                  zero_probe_output, particles, probe_input,
                                                                  probe_output, i);
        const Real block_diff = hostGradDivBlockColumnMaxAbsDiff(graddiv_blocks[i], fd_block);

        if (isGradStencilSafeFromInterface(positions[i], x_interface, dp_0))
        {
            if (metrics.core_pc_probed >= max_probe_per_zone)
            {
                continue;
            }
            metrics.core_pc_probed += 1;
            metrics.core_pc_max_abs_diff = std::max(metrics.core_pc_max_abs_diff, block_diff);
            continue;
        }

        if (!isAwayFromInterface(positions[i], x_interface, dp_0))
        {
            if (metrics.interface_pc_probed >= max_probe_per_zone)
            {
                continue;
            }
            metrics.interface_pc_probed += 1;
            metrics.interface_pc_max_abs_diff = std::max(metrics.interface_pc_max_abs_diff, block_diff);
        }
    }
}

inline AphiContactGradDivPcConsistencyMetrics runContactGradDivPcConsistencyMetrics(
    int ac, char *av[], AphiContactADivergencePenaltyStencilMode penalty_stencil =
                              AphiContactADivergencePenaltyStencilMode::InnerContact)
{
    using MainExecutionPolicy = execution::MainExecutionPolicy;

    const Real dp_0 = 0.1;
    const Real body_length = 1.0;
    const Real body_height = 1.0;
    const Real body_width = 1.0;
    const Real boundary_width = 3.0 * dp_0;
    const Real core_shell = dp_0;
    const Real x_interface = 0.5;
    const Real sigma_left = 10.0;
    const Real sigma_right = 1.0;
    const Real nu = 1.5;
    const Real omega = 1.25;
    const size_t max_probe_per_zone = 8;

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
    const AphiBlockNames probe_input{"ProbeAReal", "ProbeAImag", "ProbePhiReal", "ProbePhiImag"};
    const AphiBlockNames probe_output{"ProbeLhsAReal", "ProbeLhsAImag", "ProbePhiReal", "ProbePhiImag"};
    const AphiBlockJacobiDiagonalNames diag_names;
    const AphiADivergencePenaltyScratchNames scratch = aphiDefaultADivergencePenaltyScratchNames();

    AphiLhsAssemblyOptions options;
    options.omega = omega;
    options.use_a_divergence_penalty = true;
    options.a_divergence_penalty_mode = AphiADivergencePenaltyMode::PairwiseUncorrected;
    options.contact_a_divergence_penalty_stencil = penalty_stencil;

    for (auto *body_ptr : {&left_body, &right_body})
    {
        registerProbePenaltyFields(body_ptr->getBaseParticles(), probe_input, probe_output);
    }

    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_left(left_body, sigma_left, nu, names);
    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_right(right_body, sigma_right, nu, names);
    StateDynamics<MainExecutionPolicy, SetAphiMaterialPropertiesCK> set_left_material(left_body, sigma_left, nu, names.material);
    StateDynamics<MainExecutionPolicy, SetAphiMaterialPropertiesCK> set_right_material(right_body, sigma_right, nu, names.material);
    StateDynamics<MainExecutionPolicy, AssignConstantMaterialSigmaCK> assign_left_sigma(left_body, sigma_left, names.material);
    StateDynamics<MainExecutionPolicy, AssignConstantMaterialSigmaCK> assign_right_sigma(right_body, sigma_right, names.material);
    StateDynamics<MainExecutionPolicy, AphiZeroBlockCK> zero_left_probe_input(left_body, probe_input);
    StateDynamics<MainExecutionPolicy, AphiZeroBlockCK> zero_left_probe_output(left_body, probe_output);
    StateDynamics<MainExecutionPolicy, AphiZeroBlockCK> zero_right_probe_input(right_body, probe_input);
    StateDynamics<MainExecutionPolicy, AphiZeroBlockCK> zero_right_probe_output(right_body, probe_output);

    UpdateCellLinkedList<MainExecutionPolicy, RealBody> update_left_cell_linked_list(left_body);
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> update_right_cell_linked_list(right_body);
    UpdateRelation<MainExecutionPolicy, Inner<>> update_left_inner(left_inner);
    UpdateRelation<MainExecutionPolicy, Inner<>> update_right_inner(right_inner);
    UpdateRelation<MainExecutionPolicy, Contact<>> update_left_contact(left_to_right);
    UpdateRelation<MainExecutionPolicy, Contact<>> update_right_contact(right_to_left);

    AphiComputeBlockJacobiContactDynamicsBundle<MainExecutionPolicy> compute_left_jacobi(
        left_body, left_inner, left_to_right, names.material, omega, options, diag_names);
    AphiComputeBlockJacobiContactDynamicsBundle<MainExecutionPolicy> compute_right_jacobi(
        right_body, right_inner, right_to_left, names.material, omega, options, diag_names);

    AphiContactPairwiseADivergencePenaltyPipelineBundle<MainExecutionPolicy> left_penalty_contact(
        left_body, left_inner, left_to_right, probe_input, probe_output, Real(1), scratch);
    AphiContactPairwiseADivergencePenaltyPipelineBundle<MainExecutionPolicy> right_penalty_contact(
        right_body, right_inner, right_to_left, probe_input, probe_output, Real(1), scratch);
    AphiPairwiseADivergencePenaltyPipelineBundle<MainExecutionPolicy> left_penalty_inner_only(
        left_body, left_inner, probe_input, probe_output, Real(1), scratch);
    AphiPairwiseADivergencePenaltyPipelineBundle<MainExecutionPolicy> right_penalty_inner_only(
        right_body, right_inner, probe_input, probe_output, Real(1), scratch);

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

    AphiContactGradDivPcConsistencyMetrics metrics;
    if (penalty_stencil == AphiContactADivergencePenaltyStencilMode::InnerOnly)
    {
        probeContactGradDivPcOnBody(left_body, left_penalty_inner_only, zero_left_probe_input, zero_left_probe_output,
                                    left_body.getBaseParticles(), probe_input, probe_output, diag_names, body_length,
                                    body_height, body_width, core_shell, x_interface, dp_0, max_probe_per_zone,
                                    metrics);
        probeContactGradDivPcOnBody(right_body, right_penalty_inner_only, zero_right_probe_input,
                                    zero_right_probe_output, right_body.getBaseParticles(), probe_input, probe_output,
                                    diag_names, body_length, body_height, body_width, core_shell, x_interface, dp_0,
                                    max_probe_per_zone, metrics);
    }
    else
    {
        probeContactGradDivPcOnBody(left_body, left_penalty_contact, zero_left_probe_input, zero_left_probe_output,
                                    left_body.getBaseParticles(), probe_input, probe_output, diag_names, body_length,
                                    body_height, body_width, core_shell, x_interface, dp_0, max_probe_per_zone,
                                    metrics);
        probeContactGradDivPcOnBody(right_body, right_penalty_contact, zero_right_probe_input,
                                    zero_right_probe_output, right_body.getBaseParticles(), probe_input, probe_output,
                                    diag_names, body_length, body_height, body_width, core_shell, x_interface, dp_0,
                                    max_probe_per_zone, metrics);
    }
    return metrics;
}

inline void printContactGradDivPcConsistencyMetrics(const std::string &prefix,
                                                    const AphiContactGradDivPcConsistencyMetrics &metrics)
{
    std::cout << prefix << " core_pc_probed=" << metrics.core_pc_probed
              << " core_pc_max_abs_diff=" << metrics.core_pc_max_abs_diff
              << " interface_pc_probed=" << metrics.interface_pc_probed
              << " interface_pc_max_abs_diff=" << metrics.interface_pc_max_abs_diff << std::endl;
}

} // namespace test
} // namespace electromagnetics
} // namespace SPH

#endif // APHI_CONTACT_GRADDIV_PC_DIAGNOSTIC_HELPERS_H
