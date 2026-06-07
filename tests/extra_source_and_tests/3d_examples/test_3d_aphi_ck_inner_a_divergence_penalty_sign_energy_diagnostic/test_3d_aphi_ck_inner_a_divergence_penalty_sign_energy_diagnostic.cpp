/**
 * Stage 10.7: sign/energy diagnostic for +grad(div A) vs -grad(div A).
 */
#include "electromagnetic_dynamics/aphi_a_divergence_penalty_ck.hpp"
#include "electromagnetic_dynamics/diagnostics/aphi_a_gauge_diagnostic_helpers.h"
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
          dv_a_imag_(particles_->template getVariableByName<Vecd>(block_names.a_imag)),
          dv_phi_real_(particles_->template getVariableByName<Real>(block_names.phi_real)),
          dv_phi_imag_(particles_->template getVariableByName<Real>(block_names.phi_imag))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : position_(encloser.dv_position_->DelegatedData(ex_policy)),
              a_real_(encloser.dv_a_real_->DelegatedData(ex_policy)),
              a_imag_(encloser.dv_a_imag_->DelegatedData(ex_policy)),
              phi_real_(encloser.dv_phi_real_->DelegatedData(ex_policy)),
              phi_imag_(encloser.dv_phi_imag_->DelegatedData(ex_policy))
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
            phi_real_[index_i] = 0.0;
            phi_imag_[index_i] = 0.0;
        }

      protected:
        Vecd *position_;
        Vecd *a_real_;
        Vecd *a_imag_;
        Real *phi_real_;
        Real *phi_imag_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_position_;
    DiscreteVariable<Vecd> *dv_a_real_;
    DiscreteVariable<Vecd> *dv_a_imag_;
    DiscreteVariable<Real> *dv_phi_real_;
    DiscreteVariable<Real> *dv_phi_imag_;
};

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
    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_aphi(test_body.body, sigma, nu, names);
    StateDynamics<MainExecutionPolicy, SetAphiMaterialPropertiesCK> set_material(test_body.body, sigma, nu, names.material);
    StateDynamics<MainExecutionPolicy, AssignDivergentAFieldCK> assign_field(test_body.body, names.solution);

    InteractionDynamicsCK<MainExecutionPolicy, LinearCorrectionMatrix<Inner<WithUpdate>>> linear_correction_matrix(
        test_body.inner());
    InteractionDynamicsCK<MainExecutionPolicy, LinearGradient<Inner<Vecd>>> a_real_gradient(test_body.inner(),
                                                                                            names.solution.a_real);
    InteractionDynamicsCK<MainExecutionPolicy, LinearGradient<Inner<Vecd>>> a_imag_gradient(test_body.inner(),
                                                                                            names.solution.a_imag);
    StateDynamics<MainExecutionPolicy, AphiVectorGradientDivergenceCK> div_a_real(
        test_body.body, names.solution.a_real + "Gradient", names.diagnostic.div_a_real);
    StateDynamics<MainExecutionPolicy, AphiVectorGradientDivergenceCK> div_a_imag(
        test_body.body, names.solution.a_imag + "Gradient", names.diagnostic.div_a_imag);
    InteractionDynamicsCK<MainExecutionPolicy, LinearGradient<Inner<Real>>> div_a_real_gradient(
        test_body.inner(), names.diagnostic.div_a_real);
    InteractionDynamicsCK<MainExecutionPolicy, LinearGradient<Inner<Real>>> div_a_imag_gradient(
        test_body.inner(), names.diagnostic.div_a_imag);
    StateDynamics<MainExecutionPolicy, AphiGradDivAPenaltyCK> grad_div_penalty(
        test_body.body, names.diagnostic.div_a_real + "Gradient", names.diagnostic.div_a_imag + "Gradient", names.v,
        Real(1.0));

    initialize_aphi.exec();
    set_material.exec();
    assign_field.exec();
    test_body.updateRelations();
    linear_correction_matrix.exec();
    a_real_gradient.exec();
    a_imag_gradient.exec();
    div_a_real.exec();
    div_a_imag.exec();
    div_a_real_gradient.exec();
    div_a_imag_gradient.exec();
    StateDynamics<MainExecutionPolicy, AphiZeroBlockCK> zero_v(test_body.body, names.v);
    zero_v.exec();
    grad_div_penalty.exec();

    BaseParticles &particles = test_body.body.getBaseParticles();
    const size_t total_real_particles = particles.TotalRealParticles();
    const AphiGradDivSignEnergyMetrics metrics = hostGradDivSignEnergyMetrics(
        particles, names, total_real_particles, body_length, body_height, body_width, core_shell);

    const bool sign_ok =
        metrics.div_a_norm_squared > TinyReal && metrics.inner_a_plus_graddiv_a > 0.0 && metrics.inner_a_minus_graddiv_a < 0.0;
    const bool passed = sign_ok;

    std::cout << "test_3d_aphi_ck_inner_a_divergence_penalty_sign_energy_diagnostic"
              << " inner_A_plus_graddiv_A=" << metrics.inner_a_plus_graddiv_a
              << " inner_A_minus_graddiv_A=" << metrics.inner_a_minus_graddiv_a
              << " div_A_norm2=" << metrics.div_a_norm_squared << " ratio_correct=" << metrics.ratio_minus
              << " ratio_wrong=" << metrics.ratio_plus << " sign_ok=" << (sign_ok ? 1 : 0)
              << " core_particles=" << metrics.particle_count << " passed=" << (passed ? 1 : 0) << std::endl;

    return passed ? 0 : 1;
}
