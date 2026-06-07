/**
 * Stage 8C: tiny host-assembled scalar phi_real matrix diagnostic (Vol-weighted symmetry / nullspace / eigenvalues).
 */
#include "electromagnetic_dynamics/test_helpers/aphi_lhs_test_helpers.h"
#include "electromagnetic_dynamics/diagnostics/aphi_scalar_phi_diagnostic_helpers.h"

#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics;
using namespace SPH::electromagnetics::test;

namespace
{

Eigen::MatrixXd assembleScalarPhiRealMatrix(SPHBody &body, Inner<> &inner, AphiVariableNames &names,
                                            BaseParticles &particles, size_t total_real_particles,
                                            const AphiLhsAssemblyOptions &options)
{
    AphiApplyDynamicsBundle<MainExecutionPolicy> apply_operator(body, inner, names.solution, names.lhs, names.material,
                                                                options.omega, options);
    Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(static_cast<Eigen::Index>(total_real_particles),
                                                   static_cast<Eigen::Index>(total_real_particles));

    for (size_t basis_index = 0; basis_index != total_real_particles; ++basis_index)
    {
        StateDynamics<MainExecutionPolicy, AssignScalarPhiBasisVectorCK> assign_basis(
            body, names.solution, basis_index, false);
        assign_basis.exec();
        apply_operator.exec();

        syncAphiBlockToHost(particles, names.lhs);
        const Real *lhs_phi_real = particles.getVariableDataByName<Real>(names.lhs.phi_real);
        for (size_t row = 0; row != total_real_particles; ++row)
        {
            matrix(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(basis_index)) = lhs_phi_real[row];
        }
    }
    return matrix;
}

void reportMatrixDiagnostics(BaseParticles &particles, size_t total_real_particles, const Eigen::MatrixXd &matrix,
                             Real penalty_label)
{
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *vol_data = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Eigen::VectorXd volumes(static_cast<Eigen::Index>(total_real_particles));
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        volumes(static_cast<Eigen::Index>(i)) = vol_data[i];
    }

    const HostWeightedScalarMatrixDiagnostics diagnostics = analyzeHostWeightedScalarMatrix(matrix, volumes);
    std::cout << " penalty=" << penalty_label << " weighted_sym_gap=" << diagnostics.weighted_sym_gap
              << " row_null_gap=" << diagnostics.row_null_gap << " left_null_gap=" << diagnostics.left_null_gap
              << " min_eig_sym=" << diagnostics.min_eig_sym << " max_eig_sym=" << diagnostics.max_eig_sym;
}

} // namespace

int main(int ac, char *av[])
{
    const Real dp_0 = 0.25;
    const Real body_length = 1.0;
    const Real body_height = 1.0;
    const Real body_width = 1.0;
    const Real boundary_width = 3.0 * dp_0;
    const Real sigma = 2.0;
    const Real nu = 1.5;

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

    std::cout << "test_3d_aphi_ck_scalar_phi_tiny_host_matrix_diagnostic"
              << " dp_0=" << dp_0 << " total_real_particles=" << total_real_particles;

    const Eigen::MatrixXd matrix_no_penalty =
        assembleScalarPhiRealMatrix(test_body.body, test_body.inner(), names, particles, total_real_particles,
                                    scalarPhiLaplacePenaltyOptions(Real(0), false));
    reportMatrixDiagnostics(particles, total_real_particles, matrix_no_penalty, Real(0));

    const Eigen::MatrixXd matrix_with_penalty =
        assembleScalarPhiRealMatrix(test_body.body, test_body.inner(), names, particles, total_real_particles,
                                    scalarPhiLaplacePenaltyOptions(Real(10), true));
    reportMatrixDiagnostics(particles, total_real_particles, matrix_with_penalty, Real(10));

    std::cout << " passed=1" << std::endl;
    return 0;
}
