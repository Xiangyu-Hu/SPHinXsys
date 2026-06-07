/**
 * Stage 10.7: B=curl A manufactured diagnostic, A=(-y/2,x/2,0) -> B=(0,0,1), divA~0.
 */
#include "electromagnetic_dynamics/diagnostics/aphi_a_gauge_diagnostic_helpers.h"
#include "electromagnetic_dynamics/diagnostics/aphi_curl_a_diagnostic_ck.hpp"
#include "electromagnetic_dynamics/test_helpers/aphi_lhs_test_helpers.h"

#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics;
using namespace SPH::electromagnetics::test;

namespace
{

class AssignSoloidalAFieldCK : public LocalDynamics
{
  public:
    explicit AssignSoloidalAFieldCK(SPHBody &sph_body, const AphiBlockNames &block_names)
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

            a_real_[index_i] = Vecd(-0.5 * y, 0.5 * x, 0.0);
            a_imag_[index_i] = Vecd(0.0, 0.0, 0.0);
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
    const Real max_b_relative_error = 0.15;
    const Real max_div_a_relative = 0.2;
    const Vecd reference_b_real(0.0, 0.0, 1.0);
    const Vecd reference_b_imag(0.0, 0.0, 0.0);
    const std::string b_real_name = "BReal";
    const std::string b_imag_name = "BImag";

    AphiLhsTestBody test_body(dp_0, body_length, body_height, body_width, boundary_width, ac, av);
    IOEnvironment io_environment(test_body.sph_system);

    AphiVariableNames names;
    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_aphi(test_body.body, sigma, nu, names);
    StateDynamics<MainExecutionPolicy, SetAphiMaterialPropertiesCK> set_material(test_body.body, sigma, nu, names.material);
    StateDynamics<MainExecutionPolicy, AssignSoloidalAFieldCK> assign_field(test_body.body, names.solution);

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
    StateDynamics<MainExecutionPolicy, AphiVectorGradientCurlCK> curl_a_real(test_body.body,
                                                                             names.solution.a_real + "Gradient",
                                                                             b_real_name);
    StateDynamics<MainExecutionPolicy, AphiVectorGradientCurlCK> curl_a_imag(test_body.body,
                                                                             names.solution.a_imag + "Gradient",
                                                                             b_imag_name);

    initialize_aphi.exec();
    set_material.exec();
    assign_field.exec();
    test_body.updateRelations();
    linear_correction_matrix.exec();
    a_real_gradient.exec();
    a_imag_gradient.exec();
    div_a_real.exec();
    div_a_imag.exec();
    curl_a_real.exec();
    curl_a_imag.exec();

    BaseParticles &particles = test_body.body.getBaseParticles();
    const size_t total_real_particles = particles.TotalRealParticles();
    const AphiCurlAManufacturedErrorMetrics metrics = hostCurlAManufacturedErrorMetrics(
        particles, names, b_real_name, b_imag_name, reference_b_real, reference_b_imag, total_real_particles,
        body_length, body_height, body_width, core_shell);

    const bool passed =
        metrics.b_relative_error <= max_b_relative_error && metrics.div_a_relative <= max_div_a_relative;

    std::cout << "test_3d_aphi_ck_curl_a_manufactured_diagnostic"
              << " b_relative_error=" << metrics.b_relative_error << " b_error_L2=" << metrics.b_error_l2
              << " b_error_Linf=" << metrics.b_error_linf << " div_A_relative=" << metrics.div_a_relative
              << " core_particles=" << metrics.particle_count << " passed=" << (passed ? 1 : 0) << std::endl;

    return passed ? 0 : 1;
}
