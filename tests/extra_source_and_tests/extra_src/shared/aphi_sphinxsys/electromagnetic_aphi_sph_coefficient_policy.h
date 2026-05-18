#ifndef ELECTROMAGNETIC_APHI_SPH_COEFFICIENT_POLICY_H
#define ELECTROMAGNETIC_APHI_SPH_COEFFICIENT_POLICY_H

namespace SPH
{
namespace electromagnetics
{
namespace sph
{

/** Harmonic mean of nodal material coefficients across an interface (A-phi validated policy). */
inline Real harmonicMeanCoefficient(Real lhs, Real rhs);

/** Pairwise diffusion weight `-2 dW V / (r + reg h)` used by graph Laplace and native Laplace CK. */
inline Real pairwiseDiffusionWeight(Real dW_ijV_j, Real distance, Real pair_weight_regularization, Real reference_smoothing_length);

} // namespace sph
} // namespace electromagnetics
} // namespace SPH

#include "aphi_sphinxsys/electromagnetic_aphi_sph_coefficient_policy.hpp"

#endif // ELECTROMAGNETIC_APHI_SPH_COEFFICIENT_POLICY_H
