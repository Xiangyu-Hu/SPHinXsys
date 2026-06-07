/**
 * Stage 10.8: Jacobi pairwise 3x3 grad(div A) block vs pairwise penalty apply FD diagnostic.
 */
#include "electromagnetic_dynamics/test_helpers/aphi_lhs_test_helpers.h"
#include "electromagnetic_dynamics/aphi_pairwise_a_divergence_penalty_pipeline.h"
#include "electromagnetic_dynamics/diagnostics/aphi_graddiv_block_diagnostic_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_test_device_sync.h"

#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics;
using namespace SPH::electromagnetics::test;

namespace
{

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

} // namespace

int main(int ac, char *av[])
{
    const Real dp_0 = 0.1;
    const Real body_length = 1.0;
    const Real body_height = 1.0;
    const Real body_width = 1.0;
    const Real boundary_width = 3.0 * dp_0;
    const Real core_shell = 2.0 * dp_0;
    const Real sigma = 2.0;
    const Real nu = 1.5;
    const Real omega = 1.25;
    const Real a_divergence_penalty = 50.0;
    const size_t max_probe_particles = 8;

    AphiLhsTestBody test_body(dp_0, body_length, body_height, body_width, boundary_width, ac, av);
    IOEnvironment io_environment(test_body.sph_system);

    AphiVariableNames names;
    const AphiBlockNames probe_input{"ProbeAReal", "ProbeAImag", "ProbePhiReal", "ProbePhiImag"};
    const AphiBlockNames probe_output{"ProbeLhsAReal", "ProbeLhsAImag", "ProbePhiReal", "ProbePhiImag"};
    const AphiBlockJacobiDiagonalNames diag_names;

    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_aphi(test_body.body, sigma, nu, names);
    StateDynamics<MainExecutionPolicy, SetAphiMaterialPropertiesCK> set_material(test_body.body, sigma, nu, names.material);
    test_body.body.getBaseParticles().registerStateVariable<Vecd>(probe_input.a_real, ZeroData<Vecd>::value);
    test_body.body.getBaseParticles().registerStateVariable<Vecd>(probe_input.a_imag, ZeroData<Vecd>::value);
    test_body.body.getBaseParticles().registerStateVariable<Real>(probe_input.phi_real, Real(0));
    test_body.body.getBaseParticles().registerStateVariable<Real>(probe_input.phi_imag, Real(0));
    test_body.body.getBaseParticles().registerStateVariable<Vecd>(probe_output.a_real, ZeroData<Vecd>::value);
    test_body.body.getBaseParticles().registerStateVariable<Vecd>(probe_output.a_imag, ZeroData<Vecd>::value);
    test_body.body.getBaseParticles().registerStateVariable<Real>(probe_output.phi_real, Real(0));
    test_body.body.getBaseParticles().registerStateVariable<Real>(probe_output.phi_imag, Real(0));

    AphiLhsAssemblyOptions options;
    options.omega = omega;
    options.use_a_divergence_penalty = true;
    options.a_divergence_penalty = a_divergence_penalty;
    options.a_divergence_penalty_mode = AphiADivergencePenaltyMode::PairwiseUncorrected;

    InteractionDynamicsCK<MainExecutionPolicy, AphiComputeBlockJacobiDiagonalCK<Inner<>>> compute_jacobi_diagonal(
        test_body.inner(), names.material, omega, options, diag_names);
    AphiPairwiseADivergencePenaltyPipelineBundle<MainExecutionPolicy> penalty_pipeline(
        test_body.body, test_body.inner(), probe_input, probe_output, Real(1));
    StateDynamics<MainExecutionPolicy, AphiZeroBlockCK> zero_probe_input(test_body.body, probe_input);
    StateDynamics<MainExecutionPolicy, AphiZeroBlockCK> zero_probe_output(test_body.body, probe_output);

    initialize_aphi.exec();
    set_material.exec();
    test_body.updateRelations();
    compute_jacobi_diagonal.exec();

    BaseParticles &particles = test_body.body.getBaseParticles();
    syncVariableToHost<Vecd>(particles, "Position");
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    const size_t total_real_particles = particles.TotalRealParticles();
    syncVariableToHost<Matd>(particles, diag_names.graddiv_a_block);
    const Matd *graddiv_blocks = particles.getVariableDataByName<Matd>(diag_names.graddiv_a_block);

    Real global_max_diff = 0.0;
    size_t probed = 0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!isCoreParticle(positions[i], body_length, body_height, body_width, core_shell))
        {
            continue;
        }
        if (probed >= max_probe_particles)
        {
            break;
        }
        probed += 1;

        Matd fd_block = Matd::Zero();
        for (UnsignedInt basis = 0; basis < 3; ++basis)
        {
            zero_probe_input.exec();
            zero_probe_output.exec();
            StateDynamics<MainExecutionPolicy, AssignProbeUnitBasisAFieldCK> assign_probe(
                test_body.body, probe_input.a_real, i, basis);
            assign_probe.exec();
            penalty_pipeline.exec();
            syncVariableToHost<Vecd>(particles, probe_output.a_real);
            const Vecd *lhs_a = particles.getVariableDataByName<Vecd>(probe_output.a_real);
            const Vecd response = -lhs_a[i];
            for (UnsignedInt row = 0; row < 3; ++row)
            {
                fd_block(row, basis) = response[row];
            }
        }
        global_max_diff = std::max(global_max_diff, hostGradDivBlockColumnMaxAbsDiff(graddiv_blocks[i], fd_block));
    }

    const Real penalty_to_laplace_ratio =
        hostCorePenaltyToLaplaceDiagRatio(particles, positions, total_real_particles, body_length, body_height,
                                          body_width, core_shell, diag_names, a_divergence_penalty);
    const bool passed = probed > 0 && global_max_diff < 0.05;

    std::cout << "test_3d_aphi_ck_graddiv_block_diagonal_diagnostic"
              << " probed=" << probed << " core_max_diff=" << global_max_diff
              << " penalty_to_laplace_diag_ratio=" << penalty_to_laplace_ratio << " passed=" << (passed ? 1 : 0)
              << std::endl;
    return passed ? 0 : 1;
}
