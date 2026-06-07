/**
 * Stage 10.8: compare B-corrected vs pairwise divA / grad(divA) discretizations and energy sign.
 */
#include "electromagnetic_dynamics/diagnostics/aphi_div_a_discretization_comparison_helpers.h"
#include "electromagnetic_dynamics/aphi_a_divergence_penalty_ck.h"
#include "electromagnetic_dynamics/test_helpers/aphi_lhs_test_helpers.h"

#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics;
using namespace SPH::electromagnetics::test;

namespace
{

class AssignDivergentAFieldCK : public LocalDynamics
{
  public:
    explicit AssignDivergentAFieldCK(SPHBody &sph_body, const AphiBlockNames &block_names)
        : LocalDynamics(sph_body),
          dv_position_(particles_->template getVariableByName<Vecd>("Position")),
          dv_a_real_(particles_->template getVariableByName<Vecd>(block_names.a_real)),
          dv_a_imag_(particles_->template getVariableByName<Vecd>(block_names.a_imag))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : position_(encloser.dv_position_->DelegatedData(ex_policy)),
              a_real_(encloser.dv_a_real_->DelegatedData(ex_policy)),
              a_imag_(encloser.dv_a_imag_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            const Vecd &position = position_[index_i];
            const Real x = position[0];
            const Real y = position[1];
            const Real z = position[2];
            const Real pi = Pi;

            a_real_[index_i] = Vecd(std::sin(pi * x), std::sin(pi * y), std::sin(pi * z));
            a_imag_[index_i] = Vecd(std::cos(pi * x), std::cos(pi * y), std::cos(pi * z));
        }

      protected:
        Vecd *position_;
        Vecd *a_real_;
        Vecd *a_imag_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_position_;
    DiscreteVariable<Vecd> *dv_a_real_;
    DiscreteVariable<Vecd> *dv_a_imag_;
};

class AssignManufacturedAFieldCK : public LocalDynamics
{
  public:
    explicit AssignManufacturedAFieldCK(SPHBody &sph_body, const AphiBlockNames &block_names)
        : LocalDynamics(sph_body),
          dv_position_(particles_->template getVariableByName<Vecd>("Position")),
          dv_a_real_(particles_->template getVariableByName<Vecd>(block_names.a_real)),
          dv_a_imag_(particles_->template getVariableByName<Vecd>(block_names.a_imag))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : position_(encloser.dv_position_->DelegatedData(ex_policy)),
              a_real_(encloser.dv_a_real_->DelegatedData(ex_policy)),
              a_imag_(encloser.dv_a_imag_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            const Vecd &position = position_[index_i];
            a_real_[index_i] = position;
            a_imag_[index_i] = Vecd(position[1], position[2], position[0]);
        }

      protected:
        Vecd *position_;
        Vecd *a_real_;
        Vecd *a_imag_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_position_;
    DiscreteVariable<Vecd> *dv_a_real_;
    DiscreteVariable<Vecd> *dv_a_imag_;
};

class CopyAphiBlockCK : public LocalDynamics
{
  public:
    CopyAphiBlockCK(SPHBody &sph_body, const AphiBlockNames &source_block, const AphiBlockNames &target_block)
        : LocalDynamics(sph_body),
          dv_source_a_real_(particles_->template getVariableByName<Vecd>(source_block.a_real)),
          dv_source_a_imag_(particles_->template getVariableByName<Vecd>(source_block.a_imag)),
          dv_target_a_real_(particles_->template getVariableByName<Vecd>(target_block.a_real)),
          dv_target_a_imag_(particles_->template getVariableByName<Vecd>(target_block.a_imag))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : source_a_real_(encloser.dv_source_a_real_->DelegatedData(ex_policy)),
              source_a_imag_(encloser.dv_source_a_imag_->DelegatedData(ex_policy)),
              target_a_real_(encloser.dv_target_a_real_->DelegatedData(ex_policy)),
              target_a_imag_(encloser.dv_target_a_imag_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            target_a_real_[index_i] = source_a_real_[index_i];
            target_a_imag_[index_i] = source_a_imag_[index_i];
        }

      protected:
        Vecd *source_a_real_;
        Vecd *source_a_imag_;
        Vecd *target_a_real_;
        Vecd *target_a_imag_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_source_a_real_;
    DiscreteVariable<Vecd> *dv_source_a_imag_;
    DiscreteVariable<Vecd> *dv_target_a_real_;
    DiscreteVariable<Vecd> *dv_target_a_imag_;
};

void registerRouteProbeBlock(BaseParticles &particles, const AphiBlockNames &block)
{
    particles.registerStateVariable<Vecd>(block.a_real, ZeroData<Vecd>::value);
    particles.registerStateVariable<Vecd>(block.a_imag, ZeroData<Vecd>::value);
    particles.registerStateVariable<Real>(aphiDivAFieldName(block.a_real), Real(0));
    particles.registerStateVariable<Real>(aphiDivAFieldName(block.a_imag), Real(0));
    particles.registerStateVariable<Vecd>(aphiGradDivAFieldName(block.a_real), ZeroData<Vecd>::value);
    particles.registerStateVariable<Vecd>(aphiGradDivAFieldName(block.a_imag), ZeroData<Vecd>::value);
    particles.registerStateVariable<Matd>(block.a_real + "Gradient", ZeroData<Matd>::value);
    particles.registerStateVariable<Matd>(block.a_imag + "Gradient", ZeroData<Matd>::value);
}

AphiDivADiscretizationComparisonMetrics runCase(const char *case_name, AphiLhsTestBody &test_body,
                                                const AphiBlockNames &source_block, const AphiBlockNames &block_b,
                                                const AphiBlockNames &block_p, Real body_length, Real body_height,
                                                Real body_width, Real core_shell, bool use_manufactured_field)
{
    (void)case_name;
    BaseParticles &particles = test_body.body.getBaseParticles();
    registerRouteProbeBlock(particles, block_b);
    registerRouteProbeBlock(particles, block_p);

    StateDynamics<MainExecutionPolicy, CopyAphiBlockCK> copy_to_b(test_body.body, source_block, block_b);
    StateDynamics<MainExecutionPolicy, CopyAphiBlockCK> copy_to_p(test_body.body, source_block, block_p);
    copy_to_b.exec();
    copy_to_p.exec();
    test_body.updateRelations();

    execBCorrectedDivAGradDivAPipeline(test_body.body, test_body.inner(), block_b);
    execPairwiseDivAGradDivAPipeline(test_body.body, test_body.inner(), block_p);

    const size_t total_real_particles = particles.TotalRealParticles();
    return hostDivADiscretizationComparisonMetrics(
        particles, source_block, aphiDivAFieldName(block_b.a_real), aphiDivAFieldName(block_b.a_imag),
        aphiGradDivAFieldName(block_b.a_real), aphiGradDivAFieldName(block_b.a_imag),
        aphiDivAFieldName(block_p.a_real), aphiDivAFieldName(block_p.a_imag),
        aphiGradDivAFieldName(block_p.a_real), aphiGradDivAFieldName(block_p.a_imag), total_real_particles,
        body_length, body_height, body_width, core_shell);
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
    const Real sigma = 2.0;
    const Real nu = 1.5;

    AphiLhsTestBody test_body(dp_0, body_length, body_height, body_width, boundary_width, ac, av);
    IOEnvironment io_environment(test_body.sph_system);

    AphiVariableNames names;
    const AphiBlockNames block_b{"ProbeBReal", "ProbeBImag", "ProbeBPhiReal", "ProbeBPhiImag"};
    const AphiBlockNames block_p{"ProbePReal", "ProbePImag", "ProbePPhiReal", "ProbePPhiImag"};

    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_aphi(test_body.body, sigma, nu, names);
    StateDynamics<MainExecutionPolicy, SetAphiMaterialPropertiesCK> set_material(test_body.body, sigma, nu, names.material);
    StateDynamics<MainExecutionPolicy, AssignDivergentAFieldCK> assign_divergent(test_body.body, names.solution);
    StateDynamics<MainExecutionPolicy, AssignManufacturedAFieldCK> assign_manufactured(test_body.body, names.solution);

    initialize_aphi.exec();
    set_material.exec();
    assign_divergent.exec();
    test_body.updateRelations();
    const AphiDivADiscretizationComparisonMetrics divergent_metrics =
        runCase("divergent", test_body, names.solution, block_b, block_p, body_length, body_height, body_width,
                core_shell, false);

    assign_manufactured.exec();
    test_body.updateRelations();
    const AphiDivADiscretizationComparisonMetrics manufactured_metrics =
        runCase("manufactured", test_body, names.solution, block_b, block_p, body_length, body_height, body_width,
                core_shell, true);

    const auto energy_sign_ok = [](Real sign) { return sign > -1.0e-5; };
    const bool divergent_sign_ok =
        energy_sign_ok(divergent_metrics.energy_sign_b) && energy_sign_ok(divergent_metrics.energy_sign_pairwise);
    const bool manufactured_sign_ok = energy_sign_ok(manufactured_metrics.energy_sign_b) &&
                                      energy_sign_ok(manufactured_metrics.energy_sign_pairwise);
    const bool routes_differ =
        divergent_metrics.div_a_b_vs_pairwise_l2_diff > TinyReal &&
        divergent_metrics.grad_div_a_b_vs_pairwise_l2_diff > TinyReal;
    const bool passed = divergent_metrics.core_particles > 0 && divergent_sign_ok && manufactured_sign_ok && routes_differ;

    std::cout << "test_3d_aphi_ck_div_a_discretization_comparison_diagnostic"
              << " divergent_divA_B_L2=" << divergent_metrics.div_a_b_l2
              << " divergent_divA_pairwise_L2=" << divergent_metrics.div_a_pairwise_l2
              << " divergent_divA_diff=" << divergent_metrics.div_a_b_vs_pairwise_l2_diff
              << " divergent_gradDivA_diff=" << divergent_metrics.grad_div_a_b_vs_pairwise_l2_diff
              << " divergent_energy_sign_B=" << divergent_metrics.energy_sign_b
              << " divergent_energy_sign_pairwise=" << divergent_metrics.energy_sign_pairwise
              << " manufactured_divA_diff=" << manufactured_metrics.div_a_b_vs_pairwise_l2_diff
              << " manufactured_gradDivA_diff=" << manufactured_metrics.grad_div_a_b_vs_pairwise_l2_diff
              << " manufactured_energy_sign_B=" << manufactured_metrics.energy_sign_b
              << " manufactured_energy_sign_pairwise=" << manufactured_metrics.energy_sign_pairwise
              << " core_particles=" << divergent_metrics.core_particles << " passed=" << (passed ? 1 : 0) << std::endl;

    return passed ? 0 : 1;
}
