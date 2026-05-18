#ifndef ELECTROMAGNETIC_APHI_MATRIX_FREE_OPERATORS_HPP
#define ELECTROMAGNETIC_APHI_MATRIX_FREE_OPERATORS_HPP

#include "aphi_sphinxsys/electromagnetic_aphi_matrix_free_operators.h"

namespace SPH
{
namespace electromagnetics
{
namespace matrix_free
{

inline Real harmonicMean(Real lhs, Real rhs)
{
    return 2.0 * lhs * rhs / (lhs + rhs + TinyReal);
}

inline Real pairwiseDiffusionWeight(Real diffusion_i, Real diffusion_j, Real dW_ijV_j, Real distance)
{
    return 2.0 * harmonicMean(diffusion_i, diffusion_j) * dW_ijV_j / (distance + TinyReal);
}

inline Complex pairwiseScalarLaplaceContribution(const Complex &value_i, const Complex &value_j, Real weight)
{
    return weight * (value_i - value_j);
}

inline Vec3c pairwiseScalarGradientContribution(const Complex &value_i, const Complex &value_j, const Vecd &e_ij, Real dW_ijV_j)
{
    const Complex delta = value_j - value_i;
    Vec3c contribution = Vec3c::Zero();
    for (int axis = 0; axis != Dimensions; ++axis)
    {
        contribution[axis] = delta * e_ij[axis] * dW_ijV_j;
    }
    return contribution;
}

inline Real squaredMagnitude(const Complex &value)
{
    return std::norm(value);
}

inline Real squaredMagnitude(const Vec3c &value)
{
    Real norm = 0.0;
    for (int axis = 0; axis != Dimensions; ++axis)
    {
        norm += std::norm(value[axis]);
    }
    return norm;
}

} // namespace matrix_free
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_APHI_MATRIX_FREE_OPERATORS_HPP
