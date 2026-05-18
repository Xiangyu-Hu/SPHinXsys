#ifndef ELECTROMAGNETIC_APHI_MATRIX_FREE_OPERATORS_H
#define ELECTROMAGNETIC_APHI_MATRIX_FREE_OPERATORS_H

#include "aphi_sphinxsys/electromagnetic_aphi_matrix_free_fields.h"

namespace SPH
{
namespace electromagnetics
{
namespace matrix_free
{

Real harmonicMean(Real lhs, Real rhs);

Real pairwiseDiffusionWeight(Real diffusion_i, Real diffusion_j, Real dW_ijV_j, Real distance);

Complex pairwiseScalarLaplaceContribution(const Complex &value_i, const Complex &value_j, Real weight);

Vec3c pairwiseScalarGradientContribution(const Complex &value_i, const Complex &value_j, const Vecd &e_ij, Real dW_ijV_j);

Real squaredMagnitude(const Complex &value);

Real squaredMagnitude(const Vec3c &value);

} // namespace matrix_free
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_APHI_MATRIX_FREE_OPERATORS_H
