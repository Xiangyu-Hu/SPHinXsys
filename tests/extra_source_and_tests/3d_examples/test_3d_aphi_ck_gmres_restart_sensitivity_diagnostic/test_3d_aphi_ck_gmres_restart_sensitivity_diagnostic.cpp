/**
 * Stage 9A/9D-0 diagnostic: GMRES restart dimension sensitivity (m=10/20/30/50/80).
 */
#include "electromagnetic_dynamics/test_helpers/aphi_gmres_test_helpers.h"

#include <array>
#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics;
using namespace SPH::electromagnetics::test;

namespace
{

class AssignScalarPhiFieldsCK : public LocalDynamics
{
  public:
    explicit AssignScalarPhiFieldsCK(SPHBody &sph_body, const AphiBlockNames &block_names)
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

            a_real_[index_i] = Vecd::Zero();
            a_imag_[index_i] = Vecd::Zero();
            phi_real_[index_i] = std::sin(pi * x) * std::cos(pi * y);
            phi_imag_[index_i] = std::cos(pi * x) * std::sin(pi * z);
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
    const Real sigma = 2.0;
    const Real nu = 1.5;
    const Real phi_gauge_penalty = 10.0;
    const Real tolerance = 1.0e-5;
    const UnsignedInt max_workspace_dimension = 80;
    const UnsignedInt max_outer_iterations = 20;
    const std::array<UnsignedInt, 5> restart_values = {10, 20, 30, 50, 80};

    AphiLhsTestBody test_body(dp_0, body_length, body_height, body_width, boundary_width, ac, av);
    IOEnvironment io_environment(test_body.sph_system);

    AphiVariableNames names;
    const AphiGMRESWorkspaceNames workspace = buildAphiGMRESWorkspaceNames(max_workspace_dimension);
    RegisterAphiGMRESWorkspaceCK register_gmres_workspace(test_body.body, max_workspace_dimension);

    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_aphi_variables(test_body.body, sigma, nu, names);
    StateDynamics<MainExecutionPolicy, SetAphiMaterialPropertiesCK> set_material(test_body.body, sigma, nu, names.material);
    StateDynamics<MainExecutionPolicy, AssignScalarPhiFieldsCK> assign_exact(test_body.body, names.solution);
    StateDynamics<MainExecutionPolicy, AphiCopyBlockCK> copy_lhs_to_rhs(test_body.body, names.rhs, names.lhs);
    StateDynamics<MainExecutionPolicy, AphiZeroBlockCK> zero_solution(test_body.body, names.solution);

    AphiLhsAssemblyOptions options;
    options.terms.laplace_a = false;
    options.terms.laplace_phi = true;
    options.terms.reaction = false;
    options.terms.grad_phi_coupling = false;
    options.terms.div_sigma_a_coupling = false;
    options.use_phi_gauge_penalty = true;
    options.phi_gauge_penalty = phi_gauge_penalty;

    AphiApplyDynamicsBundle<MainExecutionPolicy> apply_exact(test_body.body, test_body.inner(), names.solution, names.lhs,
                                                             names.material, options.omega, options);

    (void)register_gmres_workspace;
    initialize_aphi_variables.exec();
    set_material.exec();
    assign_exact.exec();
    test_body.updateRelations();

    apply_exact.exec();
    copy_lhs_to_rhs.exec();
    test_body.updateRelations();

    std::cout << "test_3d_aphi_ck_gmres_restart_sensitivity_diagnostic";

    for (const UnsignedInt restart_dimension : restart_values)
    {
        zero_solution.exec();
        test_body.updateRelations();

        AphiGMRESSolverOptions gmres_options =
            defaultGMRESConvergenceOptions(tolerance, restart_dimension, max_outer_iterations);
        AphiGMRESSolverCK<MainExecutionPolicy> solver(test_body.body, test_body.inner(), names, workspace, options,
                                                        gmres_options);
        const AphiGMRESResult result = solver.solve();

        std::cout << " m=" << restart_dimension << " outer=" << result.outer_iteration_count
                  << " arnoldi=" << result.arnoldi_step_count << " rel_res=" << result.final_relative_residual
                  << " true_rel_res=" << result.final_true_relative_residual
                  << " recursive_true_gap=" << result.final_recursive_true_gap
                  << " converged=" << (result.converged ? 1 : 0)
                  << " breakdown_code=" << AphiGMRESBreakdownCodeName(result.breakdown_code);
    }

    std::cout << " passed=1" << std::endl;
    return 0;
}
