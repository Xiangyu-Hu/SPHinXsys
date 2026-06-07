#ifndef APHI_CONTACT_BODYWISE_RESIDUAL_HELPERS_H
#define APHI_CONTACT_BODYWISE_RESIDUAL_HELPERS_H

#include "electromagnetic_dynamics/test_helpers/aphi_lhs_test_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_test_device_sync.h"

#include <algorithm>
#include <cmath>

namespace SPH
{
namespace electromagnetics
{
namespace test
{

struct AphiBodywiseResidualEntry
{
    Real true_rel = 0.0;
    Real rhs_norm = 0.0;
    Real solution_norm = 0.0;
    Real absolute_true_residual_norm = 0.0;
};

struct AphiBodywiseTrueRelativeBreakdown
{
    Real global_true_rel = 0.0;
    Real max_bodywise_true_rel = 0.0;
    StdVec<AphiBodywiseResidualEntry> bodies{};
};

inline constexpr Real kMinRhsNormForBodywiseRelative = 1.0e-10;

inline AphiBodywiseResidualEntry bodywiseResidualEntry(SPHBody &body, const AphiVariableNames &names)
{
    BaseParticles &particles = body.getBaseParticles();
    const size_t total_real_particles = particles.TotalRealParticles();
    AphiBodywiseResidualEntry entry;
    entry.rhs_norm = hostBlockNorm(particles, names.rhs, total_real_particles);
    entry.solution_norm = hostBlockNorm(particles, names.solution, total_real_particles);
    entry.absolute_true_residual_norm = hostBlockNorm(particles, names.true_residual, total_real_particles);
    entry.true_rel = entry.rhs_norm >= kMinRhsNormForBodywiseRelative
                         ? entry.absolute_true_residual_norm / entry.rhs_norm
                         : entry.absolute_true_residual_norm;
    return entry;
}

inline AphiBodywiseTrueRelativeBreakdown buildBodywiseTrueRelativeBreakdown(
    const StdVec<SPHBody *> &bodies, const AphiVariableNames &names)
{
    AphiBodywiseTrueRelativeBreakdown breakdown;
    breakdown.bodies.reserve(bodies.size());

    Real sum_squared = 0.0;
    Real rhs_sum_squared = 0.0;
    for (SPHBody *body_ptr : bodies)
    {
        BaseParticles &particles = body_ptr->getBaseParticles();
        const size_t total_real_particles = particles.TotalRealParticles();
        sum_squared += std::pow(hostBlockNorm(particles, names.true_residual, total_real_particles), 2);
        rhs_sum_squared += std::pow(hostBlockNorm(particles, names.rhs, total_real_particles), 2);

        breakdown.bodies.push_back(bodywiseResidualEntry(*body_ptr, names));
    }

    breakdown.global_true_rel = std::sqrt(sum_squared) / (std::sqrt(rhs_sum_squared) + TinyReal);
    breakdown.max_bodywise_true_rel = 0.0;
    for (const AphiBodywiseResidualEntry &entry : breakdown.bodies)
    {
        const Real metric =
            entry.rhs_norm >= kMinRhsNormForBodywiseRelative ? entry.true_rel : entry.absolute_true_residual_norm;
        breakdown.max_bodywise_true_rel = std::max(breakdown.max_bodywise_true_rel, metric);
    }
    return breakdown;
}

} // namespace test
} // namespace electromagnetics
} // namespace SPH

#endif // APHI_CONTACT_BODYWISE_RESIDUAL_HELPERS_H
