#ifndef APHI_LHS_TEST_HELPERS_H
#define APHI_LHS_TEST_HELPERS_H

#include "sphinxsys.h"
#include "electromagnetic_dynamics/all_electromagnetic_dynamics_ck.h"
#include "electromagnetic_dynamics/test_helpers/aphi_test_device_sync.h"

#include <algorithm>
#include <cmath>
#include <cmath>
#include <vector>

#include <Eigen/Dense>

namespace SPH
{
namespace electromagnetics
{
namespace test
{

using MainExecutionPolicy = execution::MainExecutionPolicy;

class AphiLhsBoxShape : public ComplexShape
{
  public:
    AphiLhsBoxShape(const std::string &shape_name, const Vecd &center, const Vecd &halfsize) : ComplexShape(shape_name)
    {
        add<GeometricShapeBox>(Transform(center), halfsize);
    }
};

inline Real distanceToBoundary(const Vecd &position, Real length, Real height, Real width)
{
    return std::min({position[0], length - position[0], position[1], height - position[1], position[2], width - position[2]});
}

inline bool isCoreParticle(const Vecd &position, Real length, Real height, Real width, Real core_shell)
{
    return distanceToBoundary(position, length, height, width) > core_shell;
}

struct AphiLhsTestBody
{
    SPHSystem sph_system;
    SolidBody body;
    UniquePtr<Inner<>> inner_ck;

    AphiLhsTestBody(Real dp_0, Real body_length, Real body_height, Real body_width, Real boundary_width, int ac,
                    char *av[], Real dummy_shell_width = 0.0)
        : sph_system(BoundingBoxd(Vecd(-(boundary_width + dummy_shell_width), -(boundary_width + dummy_shell_width),
                                       -(boundary_width + dummy_shell_width)),
                                  Vecd(body_length + boundary_width + dummy_shell_width,
                                       body_height + boundary_width + dummy_shell_width,
                                       body_width + boundary_width + dummy_shell_width)),
                     dp_0),
          body(sph_system,
               makeShared<AphiLhsBoxShape>("AphiLhsTestBody",
                                           Vecd(0.5 * body_length, 0.5 * body_height, 0.5 * body_width),
                                           Vecd(0.5 * body_length, 0.5 * body_height, 0.5 * body_width)))
    {
        sph_system.handleCommandlineOptions(ac, av);
        body.defineAdaptation<SPHAdaptation>(1.15, 1.0);
        body.defineMaterial<Solid>();
        body.defineBodyLevelSetShape();
        body.generateParticles<BaseParticles, Lattice>();
        sph_system.initializeSystemCellLinkedLists();
        sph_system.initializeSystemConfigurations();
        inner_ck = makeUnique<Inner<>>(body);
    }

    Inner<> &inner() { return *inner_ck; }

    void updateRelations()
    {
        UpdateCellLinkedList<MainExecutionPolicy, RealBody> update_cell_linked_list(body);
        UpdateRelation<MainExecutionPolicy, Inner<>> update_inner_relation(inner());
        update_cell_linked_list.exec();
        update_inner_relation.exec();
    }
};

inline Real blockMaxAbsResidual(BaseParticles &particles, const AphiVariableNames &names, size_t total_real_particles)
{
    syncAphiLhsRegressionToHost(particles, names);
    const Vecd *res_a_real = particles.getVariableDataByName<Vecd>(names.residual.a_real);
    const Vecd *res_a_imag = particles.getVariableDataByName<Vecd>(names.residual.a_imag);
    const Real *res_phi_real = particles.getVariableDataByName<Real>(names.residual.phi_real);
    const Real *res_phi_imag = particles.getVariableDataByName<Real>(names.residual.phi_imag);
    Real max_error = 0.0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        max_error = std::max(max_error, res_a_real[i].norm());
        max_error = std::max(max_error, res_a_imag[i].norm());
        max_error = std::max(max_error, std::abs(res_phi_real[i]));
        max_error = std::max(max_error, std::abs(res_phi_imag[i]));
    }
    return max_error;
}

inline Real hostBlockDotProduct(BaseParticles &particles, const AphiBlockNames &block_x, const AphiBlockNames &block_y,
                                size_t total_real_particles)
{
    syncAphiBlockToHost(particles, block_x);
    syncAphiBlockToHost(particles, block_y);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const Vecd *x_a_real = particles.getVariableDataByName<Vecd>(block_x.a_real);
    const Vecd *x_a_imag = particles.getVariableDataByName<Vecd>(block_x.a_imag);
    const Real *x_phi_real = particles.getVariableDataByName<Real>(block_x.phi_real);
    const Real *x_phi_imag = particles.getVariableDataByName<Real>(block_x.phi_imag);
    const Vecd *y_a_real = particles.getVariableDataByName<Vecd>(block_y.a_real);
    const Vecd *y_a_imag = particles.getVariableDataByName<Vecd>(block_y.a_imag);
    const Real *y_phi_real = particles.getVariableDataByName<Real>(block_y.phi_real);
    const Real *y_phi_imag = particles.getVariableDataByName<Real>(block_y.phi_imag);
    Real dot = 0.0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        dot += vol[i] * (x_a_real[i].dot(y_a_real[i]) + x_a_imag[i].dot(y_a_imag[i]) +
                         x_phi_real[i] * y_phi_real[i] + x_phi_imag[i] * y_phi_imag[i]);
    }
    return dot;
}

inline Real hostBlockNorm(BaseParticles &particles, const AphiBlockNames &block_names, size_t total_real_particles)
{
    return std::sqrt(hostBlockDotProduct(particles, block_names, block_names, total_real_particles));
}

inline Real coreBlockOperatorNorm(BaseParticles &particles, const AphiBlockNames &block_names, const Vecd *positions,
                                 size_t total_real_particles, Real body_length, Real body_height, Real body_width,
                                 Real core_shell)
{
    syncAphiBlockToHost(particles, block_names);
    const Vecd *a_real = particles.getVariableDataByName<Vecd>(block_names.a_real);
    const Vecd *a_imag = particles.getVariableDataByName<Vecd>(block_names.a_imag);
    const Real *phi_real = particles.getVariableDataByName<Real>(block_names.phi_real);
    const Real *phi_imag = particles.getVariableDataByName<Real>(block_names.phi_imag);
    Real max_norm = 0.0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!isCoreParticle(positions[i], body_length, body_height, body_width, core_shell))
        {
            continue;
        }
        max_norm = std::max(max_norm, a_real[i].norm());
        max_norm = std::max(max_norm, a_imag[i].norm());
        max_norm = std::max(max_norm, std::abs(phi_real[i]));
        max_norm = std::max(max_norm, std::abs(phi_imag[i]));
    }
    return max_norm;
}

inline Real coreBlockMaxAbsDifference(BaseParticles &particles, const AphiBlockNames &lhs_block,
                                      const AphiBlockNames &rhs_block, const Vecd *positions,
                                      size_t total_real_particles, Real body_length, Real body_height, Real body_width,
                                      Real core_shell)
{
    syncAphiBlockToHost(particles, lhs_block);
    syncAphiBlockToHost(particles, rhs_block);
    const Vecd *lhs_a_real = particles.getVariableDataByName<Vecd>(lhs_block.a_real);
    const Vecd *lhs_a_imag = particles.getVariableDataByName<Vecd>(lhs_block.a_imag);
    const Real *lhs_phi_real = particles.getVariableDataByName<Real>(lhs_block.phi_real);
    const Real *lhs_phi_imag = particles.getVariableDataByName<Real>(lhs_block.phi_imag);
    const Vecd *rhs_a_real = particles.getVariableDataByName<Vecd>(rhs_block.a_real);
    const Vecd *rhs_a_imag = particles.getVariableDataByName<Vecd>(rhs_block.a_imag);
    const Real *rhs_phi_real = particles.getVariableDataByName<Real>(rhs_block.phi_real);
    const Real *rhs_phi_imag = particles.getVariableDataByName<Real>(rhs_block.phi_imag);
    Real max_error = 0.0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!isCoreParticle(positions[i], body_length, body_height, body_width, core_shell))
        {
            continue;
        }
        max_error = std::max(max_error, (lhs_a_real[i] - rhs_a_real[i]).norm());
        max_error = std::max(max_error, (lhs_a_imag[i] - rhs_a_imag[i]).norm());
        max_error = std::max(max_error, std::abs(lhs_phi_real[i] - rhs_phi_real[i]));
        max_error = std::max(max_error, std::abs(lhs_phi_imag[i] - rhs_phi_imag[i]));
    }
    return max_error;
}

inline Real coreBlockMaxAbsResidual(BaseParticles &particles, const AphiVariableNames &names, const Vecd *positions,
                                    size_t total_real_particles, Real body_length, Real body_height, Real body_width,
                                    Real core_shell)
{
    syncAphiLhsRegressionToHost(particles, names);
    const Vecd *res_a_real = particles.getVariableDataByName<Vecd>(names.residual.a_real);
    const Vecd *res_a_imag = particles.getVariableDataByName<Vecd>(names.residual.a_imag);
    const Real *res_phi_real = particles.getVariableDataByName<Real>(names.residual.phi_real);
    const Real *res_phi_imag = particles.getVariableDataByName<Real>(names.residual.phi_imag);
    Real max_error = 0.0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!isCoreParticle(positions[i], body_length, body_height, body_width, core_shell))
        {
            continue;
        }
        max_error = std::max(max_error, res_a_real[i].norm());
        max_error = std::max(max_error, res_a_imag[i].norm());
        max_error = std::max(max_error, std::abs(res_phi_real[i]));
        max_error = std::max(max_error, std::abs(res_phi_imag[i]));
    }
    return max_error;
}

inline Real hostVolSum(BaseParticles &particles, size_t total_real_particles)
{
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real sum_vol = Real(0);
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        sum_vol += vol[i];
    }
    return sum_vol;
}

inline Real hostScalarPhiRealVolWeightedMean(BaseParticles &particles, const AphiBlockNames &block_names,
                                             size_t total_real_particles)
{
    syncAphiBlockToHost(particles, block_names);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const Real *phi_real = particles.getVariableDataByName<Real>(block_names.phi_real);
    Real weighted_sum = Real(0);
    Real vol_sum = Real(0);
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        weighted_sum += vol[i] * phi_real[i];
        vol_sum += vol[i];
    }
    return weighted_sum / (vol_sum + TinyReal);
}

inline Real hostScalarPhiImagVolWeightedMean(BaseParticles &particles, const AphiBlockNames &block_names,
                                             size_t total_real_particles)
{
    syncAphiBlockToHost(particles, block_names);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const Real *phi_imag = particles.getVariableDataByName<Real>(block_names.phi_imag);
    Real weighted_sum = Real(0);
    Real vol_sum = Real(0);
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        weighted_sum += vol[i] * phi_imag[i];
        vol_sum += vol[i];
    }
    return weighted_sum / (vol_sum + TinyReal);
}

struct HostWeightedScalarMatrixDiagnostics
{
    Real weighted_sym_gap = Real(0);
    Real row_null_gap = Real(0);
    Real left_null_gap = Real(0);
    Real min_eig_sym = Real(0);
    Real max_eig_sym = Real(0);
};

inline HostWeightedScalarMatrixDiagnostics analyzeHostWeightedScalarMatrix(const Eigen::MatrixXd &matrix,
                                                                           const Eigen::VectorXd &volumes)
{
    HostWeightedScalarMatrixDiagnostics diagnostics;
    const Eigen::VectorXd vol = volumes.cwiseMax(TinyReal);
    const Eigen::MatrixXd weighted_matrix = vol.asDiagonal() * matrix;
    const Eigen::MatrixXd weighted_transpose = matrix.transpose() * vol.asDiagonal();
    const Real weighted_frobenius = weighted_matrix.norm();
    diagnostics.weighted_sym_gap =
        (weighted_matrix - weighted_transpose).norm() / (weighted_frobenius + TinyReal);

    const Eigen::VectorXd ones = Eigen::VectorXd::Ones(matrix.rows());
    diagnostics.row_null_gap = (matrix * ones).norm() / (matrix.rows() + TinyReal);
    diagnostics.left_null_gap = (ones.transpose() * weighted_matrix).norm() / (ones.norm() + TinyReal);

    const Eigen::VectorXd sqrt_vol = vol.cwiseSqrt();
    const Eigen::MatrixXd normalized_matrix =
        sqrt_vol.asDiagonal() * matrix * sqrt_vol.asDiagonal().inverse();
    const Eigen::MatrixXd sym_part = Real(0.5) * (normalized_matrix + normalized_matrix.transpose());
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen_solver(sym_part);
    if (eigen_solver.info() == Eigen::Success)
    {
        diagnostics.min_eig_sym = eigen_solver.eigenvalues().minCoeff();
        diagnostics.max_eig_sym = eigen_solver.eigenvalues().maxCoeff();
    }
    return diagnostics;
}

} // namespace test
} // namespace electromagnetics
} // namespace SPH

#endif // APHI_LHS_TEST_HELPERS_H
