#include "electromagnetic_dynamics/test_helpers/aphi_lhs_test_helpers.h"

#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics;
using namespace SPH::electromagnetics::test;

namespace
{

class AssignBlockVectorOpsFieldsCK : public LocalDynamics
{
  public:
    AssignBlockVectorOpsFieldsCK(SPHBody &sph_body, const AphiBlockNames &block_x, const AphiBlockNames &block_y)
        : LocalDynamics(sph_body),
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
            : position_(encloser.dv_position_->DelegatedData(ex_policy)),
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

            x_a_real_[index_i] = Vecd(std::sin(pi * x), std::cos(pi * y), std::sin(pi * z));
            x_a_imag_[index_i] = Vecd(std::cos(pi * x), std::sin(pi * y), std::cos(pi * z));
            x_phi_real_[index_i] = std::sin(pi * x) * std::cos(pi * y);
            x_phi_imag_[index_i] = std::cos(pi * x) * std::sin(pi * z);

            y_a_real_[index_i] = Vecd(0.5 * x, 0.25 * y, 0.75 * z);
            y_a_imag_[index_i] = Vecd(-0.2 * z, 0.3 * x, 0.1 * y);
            y_phi_real_[index_i] = 0.5 * x * y;
            y_phi_imag_[index_i] = 0.25 * y * z;
        }

      protected:
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
    const Real validation_threshold = 1.0e-5;
    const Real alpha = 0.37;
    const Real coeff_x = 1.25;
    const Real coeff_y = -0.42;

    AphiLhsTestBody test_body(dp_0, body_length, body_height, body_width, boundary_width, ac, av);
    IOEnvironment io_environment(test_body.sph_system);

    AphiVariableNames names;
    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_aphi_variables(test_body.body, sigma, nu, names);
    StateDynamics<MainExecutionPolicy, SetAphiMaterialPropertiesCK> set_material(test_body.body, sigma, nu, names.material);
    StateDynamics<MainExecutionPolicy, AssignBlockVectorOpsFieldsCK> assign_fields(test_body.body, names.solution, names.v);
    StateDynamics<MainExecutionPolicy, AphiCopyBlockCK> copy_x_to_s(test_body.body, names.s, names.solution);
    StateDynamics<MainExecutionPolicy, AphiBlockAXPYCK> axpy_s(test_body.body, names.s, alpha, names.v);
    StateDynamics<MainExecutionPolicy, AphiBlockLinearCombinationCK> linear_combination_t(
        test_body.body, names.t, coeff_x, coeff_y, names.solution, names.v);

    ReduceDynamicsCK<MainExecutionPolicy, AphiBlockDotProductCK> dot_xy(test_body.body, names.solution, names.v);
    ReduceDynamicsCK<MainExecutionPolicy, AphiBlockNormSquaredCK> norm_squared_x(test_body.body, names.solution);

    initialize_aphi_variables.exec();
    set_material.exec();
    assign_fields.exec();
    test_body.updateRelations();

    BaseParticles &particles = test_body.body.getBaseParticles();
    const size_t total_real_particles = particles.TotalRealParticles();

    const Real host_dot_xy = hostBlockDotProduct(particles, names.solution, names.v, total_real_particles);
    const Real host_norm_x = hostBlockNorm(particles, names.solution, total_real_particles);

    const Real device_dot_xy = dot_xy.exec();
    const Real device_norm_x = std::sqrt(norm_squared_x.exec());

    copy_x_to_s.exec();
    axpy_s.exec();
    linear_combination_t.exec();

    syncAphiBlockToHost(particles, names.s);
    syncAphiBlockToHost(particles, names.t);
    syncAphiBlockToHost(particles, names.solution);
    syncAphiBlockToHost(particles, names.v);

    const Vecd *x_a_real = particles.getVariableDataByName<Vecd>(names.solution.a_real);
    const Vecd *x_a_imag = particles.getVariableDataByName<Vecd>(names.solution.a_imag);
    const Real *x_phi_real = particles.getVariableDataByName<Real>(names.solution.phi_real);
    const Real *x_phi_imag = particles.getVariableDataByName<Real>(names.solution.phi_imag);
    const Vecd *y_a_real = particles.getVariableDataByName<Vecd>(names.v.a_real);
    const Vecd *y_a_imag = particles.getVariableDataByName<Vecd>(names.v.a_imag);
    const Real *y_phi_real = particles.getVariableDataByName<Real>(names.v.phi_real);
    const Real *y_phi_imag = particles.getVariableDataByName<Real>(names.v.phi_imag);
    const Vecd *s_a_real = particles.getVariableDataByName<Vecd>(names.s.a_real);
    const Vecd *s_a_imag = particles.getVariableDataByName<Vecd>(names.s.a_imag);
    const Real *s_phi_real = particles.getVariableDataByName<Real>(names.s.phi_real);
    const Real *s_phi_imag = particles.getVariableDataByName<Real>(names.s.phi_imag);
    const Vecd *t_a_real = particles.getVariableDataByName<Vecd>(names.t.a_real);
    const Vecd *t_a_imag = particles.getVariableDataByName<Vecd>(names.t.a_imag);
    const Real *t_phi_real = particles.getVariableDataByName<Real>(names.t.phi_real);
    const Real *t_phi_imag = particles.getVariableDataByName<Real>(names.t.phi_imag);

    Real copy_axpy_max_error = 0.0;
    Real linear_combination_max_error = 0.0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        copy_axpy_max_error = std::max(copy_axpy_max_error, (s_a_real[i] - (x_a_real[i] + alpha * y_a_real[i])).norm());
        copy_axpy_max_error = std::max(copy_axpy_max_error, (s_a_imag[i] - (x_a_imag[i] + alpha * y_a_imag[i])).norm());
        copy_axpy_max_error =
            std::max(copy_axpy_max_error, std::abs(s_phi_real[i] - (x_phi_real[i] + alpha * y_phi_real[i])));
        copy_axpy_max_error =
            std::max(copy_axpy_max_error, std::abs(s_phi_imag[i] - (x_phi_imag[i] + alpha * y_phi_imag[i])));

        linear_combination_max_error =
            std::max(linear_combination_max_error, (t_a_real[i] - (coeff_x * x_a_real[i] + coeff_y * y_a_real[i])).norm());
        linear_combination_max_error =
            std::max(linear_combination_max_error, (t_a_imag[i] - (coeff_x * x_a_imag[i] + coeff_y * y_a_imag[i])).norm());
        linear_combination_max_error = std::max(
            linear_combination_max_error, std::abs(t_phi_real[i] - (coeff_x * x_phi_real[i] + coeff_y * y_phi_real[i])));
        linear_combination_max_error = std::max(
            linear_combination_max_error, std::abs(t_phi_imag[i] - (coeff_x * x_phi_imag[i] + coeff_y * y_phi_imag[i])));
    }

    const Real dot_error = std::abs(device_dot_xy - host_dot_xy);
    const Real norm_error = std::abs(device_norm_x - host_norm_x);
    const bool passed = dot_error < validation_threshold && norm_error < validation_threshold &&
                        copy_axpy_max_error < validation_threshold && linear_combination_max_error < validation_threshold;

    std::cout << "test_3d_aphi_ck_block_vector_ops"
              << " dot_error=" << dot_error
              << " norm_error=" << norm_error
              << " copy_axpy_max_error=" << copy_axpy_max_error
              << " linear_combination_max_error=" << linear_combination_max_error
              << " passed=" << (passed ? 1 : 0) << std::endl;

    return passed ? 0 : 1;
}
