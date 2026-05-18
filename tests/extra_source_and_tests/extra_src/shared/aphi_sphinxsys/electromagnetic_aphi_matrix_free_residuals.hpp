#ifndef ELECTROMAGNETIC_APHI_MATRIX_FREE_RESIDUALS_HPP
#define ELECTROMAGNETIC_APHI_MATRIX_FREE_RESIDUALS_HPP

#include "aphi_sphinxsys/electromagnetic_aphi_matrix_free_residuals.h"
#include <algorithm>
#include <cmath>

namespace SPH
{
namespace electromagnetics
{
namespace matrix_free
{

inline ScalarLaplaceEdge::ScalarLaplaceEdge(size_t index_i, size_t index_j, Real pair_weight)
    : index_i_(index_i), index_j_(index_j), pair_weight_(pair_weight)
{
}

inline ScalarComplexHelmholtzResiduals::ScalarComplexHelmholtzResiduals()
    : l2_norm_(0.0), mean_abs_(0.0), max_abs_(0.0)
{
}

inline void ScalarComplexHelmholtzResiduals::resize(size_t size)
{
    laplace_term_.resize(size, Complex(0.0, 0.0));
    reaction_term_.resize(size, Complex(0.0, 0.0));
    rhs_term_.resize(size, Complex(0.0, 0.0));
    residual_.resize(size, Complex(0.0, 0.0));
    diagonal_scale_.resize(size, 0.0);
    squared_residual_.resize(size, 0.0);
    clear();
}

inline void ScalarComplexHelmholtzResiduals::clear()
{
    std::fill(laplace_term_.begin(), laplace_term_.end(), Complex(0.0, 0.0));
    std::fill(reaction_term_.begin(), reaction_term_.end(), Complex(0.0, 0.0));
    std::fill(rhs_term_.begin(), rhs_term_.end(), Complex(0.0, 0.0));
    std::fill(residual_.begin(), residual_.end(), Complex(0.0, 0.0));
    std::fill(diagonal_scale_.begin(), diagonal_scale_.end(), 0.0);
    std::fill(squared_residual_.begin(), squared_residual_.end(), 0.0);
    l2_norm_ = 0.0;
    mean_abs_ = 0.0;
    max_abs_ = 0.0;
}

inline void accumulateDirectionalScalarLaplaceContribution(size_t index_i, const Complex &value_i, const Complex &value_j,
                                                           Real pair_weight, ScalarComplexHelmholtzResiduals &residuals)
{
    residuals.laplace_term_[index_i] += pairwiseScalarLaplaceContribution(value_i, value_j, pair_weight);
    residuals.diagonal_scale_[index_i] += abs(pair_weight);
}

inline void accumulateSymmetricScalarLaplaceContribution(size_t index_i, size_t index_j, const StdVec<Complex> &field,
                                                         Real pair_weight, ScalarComplexHelmholtzResiduals &residuals)
{
    const Complex pair_contribution = pairwiseScalarLaplaceContribution(field[index_i], field[index_j], pair_weight);
    residuals.laplace_term_[index_i] += pair_contribution;
    residuals.laplace_term_[index_j] -= pair_contribution;
    residuals.diagonal_scale_[index_i] += abs(pair_weight);
    residuals.diagonal_scale_[index_j] += abs(pair_weight);
}

inline void accumulateScalarLaplaceResidualsFromEdges(const StdVec<Complex> &field, const StdVec<ScalarLaplaceEdge> &edges,
                                                      ScalarComplexHelmholtzResiduals &residuals)
{
    for (const ScalarLaplaceEdge &edge : edges)
    {
        accumulateSymmetricScalarLaplaceContribution(edge.index_i_, edge.index_j_, field, edge.pair_weight_, residuals);
    }
}

inline void finalizeScalarHelmholtzResiduals(const StdVec<Complex> &field, const StdVec<Complex> &rhs,
                                             const StdVec<Complex> &reaction_coefficient,
                                             ScalarComplexHelmholtzResiduals &residuals)
{
    const size_t number_of_particles = field.size();
    residuals.l2_norm_ = 0.0;
    residuals.mean_abs_ = 0.0;
    residuals.max_abs_ = 0.0;

    for (size_t index_i = 0; index_i != number_of_particles; ++index_i)
    {
        residuals.reaction_term_[index_i] = reaction_coefficient[index_i] * field[index_i];
        residuals.rhs_term_[index_i] = rhs[index_i];
        residuals.residual_[index_i] =
            residuals.laplace_term_[index_i] + residuals.reaction_term_[index_i] - residuals.rhs_term_[index_i];

        residuals.squared_residual_[index_i] = squaredMagnitude(residuals.residual_[index_i]);
        residuals.l2_norm_ += residuals.squared_residual_[index_i];

        const Real residual_abs = sqrt(residuals.squared_residual_[index_i]);
        residuals.mean_abs_ += residual_abs;
        residuals.max_abs_ = std::max(residuals.max_abs_, residual_abs);
    }

    if (number_of_particles != 0)
    {
        residuals.l2_norm_ = sqrt(residuals.l2_norm_ / static_cast<Real>(number_of_particles));
        residuals.mean_abs_ /= static_cast<Real>(number_of_particles);
    }
}

inline Real computeScalarHelmholtzResidualL2(const ScalarComplexHelmholtzResiduals &residuals)
{
    return residuals.l2_norm_;
}

} // namespace matrix_free
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_APHI_MATRIX_FREE_RESIDUALS_HPP
