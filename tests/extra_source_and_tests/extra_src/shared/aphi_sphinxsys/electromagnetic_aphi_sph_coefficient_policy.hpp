#ifndef ELECTROMAGNETIC_APHI_SPH_COEFFICIENT_POLICY_HPP
#define ELECTROMAGNETIC_APHI_SPH_COEFFICIENT_POLICY_HPP

#include "aphi_sphinxsys/electromagnetic_aphi_sph_coefficient_policy.h"

namespace SPH
{
namespace electromagnetics
{
namespace sph
{

inline Real harmonicMeanCoefficient(Real lhs, Real rhs)
{
    return 2.0 * lhs * rhs / (lhs + rhs + TinyReal);
}

inline Real pairwiseDiffusionWeight(Real dW_ijV_j, Real distance, Real pair_weight_regularization, Real reference_smoothing_length)
{
    return -2.0 * dW_ijV_j / (distance + pair_weight_regularization * reference_smoothing_length + TinyReal);
}

} // namespace sph
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_APHI_SPH_COEFFICIENT_POLICY_HPP
