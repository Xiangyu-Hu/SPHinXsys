#ifndef APHI_GRADDIV_BLOCK_DIAGNOSTIC_HELPERS_H
#define APHI_GRADDIV_BLOCK_DIAGNOSTIC_HELPERS_H

#include "electromagnetic_dynamics/aphi_block_jacobi_preconditioner_ck.hpp"
#include "electromagnetic_dynamics/aphi_coupling_modes_ck.h"
#include "electromagnetic_dynamics/test_helpers/aphi_lhs_test_helpers.h"

namespace SPH
{
namespace electromagnetics
{
namespace test
{

inline Real hostGradDivBlockColumnMaxAbsDiff(const Matd &pairwise_block, const Matd &fd_block)
{
    Real max_abs_diff = 0.0;
    for (UnsignedInt row = 0; row < 3; ++row)
    {
        for (UnsignedInt col = 0; col < 3; ++col)
        {
            max_abs_diff = std::max(max_abs_diff, std::abs(pairwise_block(row, col) - fd_block(row, col)));
        }
    }
    return max_abs_diff;
}

inline Real hostCorePenaltyToLaplaceDiagRatio(BaseParticles &particles, const Vecd *positions,
                                              size_t total_real_particles, Real body_length, Real body_height,
                                              Real body_width, Real core_shell,
                                              const AphiBlockJacobiDiagonalNames &diag_names, Real a_divergence_penalty)
{
    syncVariableToHost<Matd>(particles, diag_names.graddiv_a_block);
    syncVariableToHost<Real>(particles, diag_names.laplace_a_diag);
    const Matd *graddiv_blocks = particles.getVariableDataByName<Matd>(diag_names.graddiv_a_block);
    const Real *laplace_a_diag = particles.getVariableDataByName<Real>(diag_names.laplace_a_diag);

    Real ratio_sum = 0.0;
    size_t core_count = 0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!isCoreParticle(positions[i], body_length, body_height, body_width, core_shell))
        {
            continue;
        }
        core_count += 1;
        ratio_sum += a_divergence_penalty * graddiv_blocks[i].norm() / (std::abs(laplace_a_diag[i]) + TinyReal);
    }
    return core_count > 0 ? ratio_sum / static_cast<Real>(core_count) : 0.0;
}

} // namespace test
} // namespace electromagnetics
} // namespace SPH

#endif // APHI_GRADDIV_BLOCK_DIAGNOSTIC_HELPERS_H
