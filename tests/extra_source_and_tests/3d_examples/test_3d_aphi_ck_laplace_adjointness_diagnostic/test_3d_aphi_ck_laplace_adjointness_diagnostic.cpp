/**
 * Stage 7: Vol-weighted adjointness diagnostic for pure Laplace operators.
 * Compares <x,Ky>_V and <Kx,y>_V without modifying pairwise Laplace weights.
 */
#include "electromagnetic_dynamics/test_helpers/aphi_lhs_test_helpers.h"

#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics;
using namespace SPH::electromagnetics::test;

namespace
{

class AssignTwoBlockFieldsCK : public LocalDynamics
{
  public:
    AssignTwoBlockFieldsCK(SPHBody &sph_body, const AphiBlockNames &block_x, const AphiBlockNames &block_y, bool vector_mode)
        : LocalDynamics(sph_body), vector_mode_(vector_mode),
          dv_position_(particles_->template getVariableByName<Vecd>("Position")),
          dv_x_a_real_(particles_->template getVariableByName<Vecd>(block_x.a_real)),
          dv_x_a_imag_(particles_->template getVariableByName<Vecd>(block_x.a_imag)),
          dv_x_phi_real_(particles_->template getVariableByName<Real>(block_x.phi_real)),
          dv_x_phi_imag_(particles_->template getVariableByName<Real>(block_x.phi_imag)),
          dv_y_a_real_(particles_->template getVariableByName<Vecd>(block_y.a_real)),
          dv_y_a_imag_(particles_->template getVariableByName<Vecd>(block_y.a_imag)),
          dv_y_phi_real_(particles_->template getVariableByName<Real>(block_y.phi_real)),
          dv_y_phi_imag_(particles_->template getVariableByName<Real>(block_y.phi_imag))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : vector_mode_(encloser.vector_mode_),
              position_(encloser.dv_position_->DelegatedData(ex_policy)),
              x_a_real_(encloser.dv_x_a_real_->DelegatedData(ex_policy)),
              x_a_imag_(encloser.dv_x_a_imag_->DelegatedData(ex_policy)),
              x_phi_real_(encloser.dv_x_phi_real_->DelegatedData(ex_policy)),
              x_phi_imag_(encloser.dv_x_phi_imag_->DelegatedData(ex_policy)),
              y_a_real_(encloser.dv_y_a_real_->DelegatedData(ex_policy)),
              y_a_imag_(encloser.dv_y_a_imag_->DelegatedData(ex_policy)),
              y_phi_real_(encloser.dv_y_phi_real_->DelegatedData(ex_policy)),
              y_phi_imag_(encloser.dv_y_phi_imag_->DelegatedData(ex_policy))
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

            if (vector_mode_)
            {
                x_a_real_[index_i] = Vecd(std::sin(pi * x), std::sin(pi * y), std::sin(pi * z));
                x_a_imag_[index_i] = Vecd(std::cos(pi * x), std::cos(pi * y), std::cos(pi * z));
                y_a_real_[index_i] = Vecd(0.5 * x, 0.25 * y, 0.75 * z);
                y_a_imag_[index_i] = Vecd(-0.2 * z, 0.3 * x, 0.1 * y);
                x_phi_real_[index_i] = 0.0;
                x_phi_imag_[index_i] = 0.0;
                y_phi_real_[index_i] = 0.0;
                y_phi_imag_[index_i] = 0.0;
            }
            else
            {
                x_a_real_[index_i] = Vecd::Zero();
                x_a_imag_[index_i] = Vecd::Zero();
                y_a_real_[index_i] = Vecd::Zero();
                y_a_imag_[index_i] = Vecd::Zero();
                x_phi_real_[index_i] = std::sin(pi * x) * std::cos(pi * y);
                x_phi_imag_[index_i] = std::cos(pi * x) * std::sin(pi * z);
                y_phi_real_[index_i] = 0.5 * x * y;
                y_phi_imag_[index_i] = 0.25 * y * z;
            }
        }

      protected:
        bool vector_mode_;
        Vecd *position_;
        Vecd *x_a_real_;
        Vecd *x_a_imag_;
        Real *x_phi_real_;
        Real *x_phi_imag_;
        Vecd *y_a_real_;
        Vecd *y_a_imag_;
        Real *y_phi_real_;
        Real *y_phi_imag_;
    };

  protected:
    bool vector_mode_;
    DiscreteVariable<Vecd> *dv_position_;
    DiscreteVariable<Vecd> *dv_x_a_real_;
    DiscreteVariable<Vecd> *dv_x_a_imag_;
    DiscreteVariable<Real> *dv_x_phi_real_;
    DiscreteVariable<Real> *dv_x_phi_imag_;
    DiscreteVariable<Vecd> *dv_y_a_real_;
    DiscreteVariable<Vecd> *dv_y_a_imag_;
    DiscreteVariable<Real> *dv_y_phi_real_;
    DiscreteVariable<Real> *dv_y_phi_imag_;
};

Real measureAdjointGap(SPHBody &body, Inner<> &inner, AphiVariableNames &names, const AphiLhsAssemblyOptions &options,
                       size_t total_real_particles)
{
    AphiApplyDynamicsBundle<MainExecutionPolicy> apply_x(body, inner, names.solution, names.lhs, names.material,
                                                       options.omega, options);
    AphiApplyDynamicsBundle<MainExecutionPolicy> apply_y(body, inner, names.v, names.t, names.material, options.omega,
                                                         options);
    ReduceDynamicsCK<MainExecutionPolicy, AphiBlockDotProductCK> dot_x_ky(body, names.solution, names.t);
    ReduceDynamicsCK<MainExecutionPolicy, AphiBlockDotProductCK> dot_kx_y(body, names.lhs, names.v);

    apply_x.exec();
    apply_y.exec();

    const Real x_ky = dot_x_ky.exec();
    const Real kx_y = dot_kx_y.exec();
    const Real scale = std::max({std::abs(x_ky), std::abs(kx_y), Real(1)});
    return std::abs(x_ky - kx_y) / scale;
}

} // namespace

int main(int ac, char *av[])
{
    const Real dp_0 = 0.1;
    const Real body_length = 1.0;
    const Real body_height = 1.0;
    const Real body_width = 1.0;
    const Real boundary_width = 3.0 * dp_0;
    const Real sigma = 2.0;
    const Real nu = 1.5;
    const Real validation_threshold = 2.0;

    AphiLhsTestBody test_body(dp_0, body_length, body_height, body_width, boundary_width, ac, av);
    IOEnvironment io_environment(test_body.sph_system);

    AphiVariableNames names;
    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_aphi_variables(test_body.body, sigma, nu, names);
    StateDynamics<MainExecutionPolicy, SetAphiMaterialPropertiesCK> set_material(test_body.body, sigma, nu, names.material);

    initialize_aphi_variables.exec();
    set_material.exec();
    test_body.updateRelations();

    BaseParticles &particles = test_body.body.getBaseParticles();
    const size_t total_real_particles = particles.TotalRealParticles();

    AphiLhsAssemblyOptions scalar_options;
    scalar_options.terms.laplace_a = false;
    scalar_options.terms.laplace_phi = true;
    scalar_options.terms.reaction = false;
    scalar_options.terms.grad_phi_coupling = false;
    scalar_options.terms.div_sigma_a_coupling = false;

    AphiLhsAssemblyOptions vector_options = scalar_options;
    vector_options.terms.laplace_a = true;
    vector_options.terms.laplace_phi = false;

    StateDynamics<MainExecutionPolicy, AssignTwoBlockFieldsCK> assign_scalar_fields(
        test_body.body, names.solution, names.v, false);
    assign_scalar_fields.exec();
    test_body.updateRelations();
    const Real scalar_gap = measureAdjointGap(test_body.body, test_body.inner(), names, scalar_options, total_real_particles);

    StateDynamics<MainExecutionPolicy, AssignTwoBlockFieldsCK> assign_vector_fields(
        test_body.body, names.solution, names.v, true);
    assign_vector_fields.exec();
    test_body.updateRelations();
    const Real vector_gap = measureAdjointGap(test_body.body, test_body.inner(), names, vector_options, total_real_particles);

    const bool passed = scalar_gap < validation_threshold && vector_gap < validation_threshold;

    std::cout << "test_3d_aphi_ck_laplace_adjointness_diagnostic"
              << " scalar_relative_gap=" << scalar_gap
              << " vector_relative_gap=" << vector_gap
              << " passed=" << (passed ? 1 : 0) << std::endl;

    return passed ? 0 : 1;
}
