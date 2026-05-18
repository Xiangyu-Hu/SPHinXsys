#ifndef ELECTROMAGNETIC_APHI_MATRIX_FREE_RESIDUALS_H
#define ELECTROMAGNETIC_APHI_MATRIX_FREE_RESIDUALS_H

#include "aphi_sphinxsys/electromagnetic_aphi_matrix_free_operators.hpp"

namespace SPH
{
namespace electromagnetics
{
namespace matrix_free
{

struct ScalarLaplaceEdge
{
    size_t index_i_;
    size_t index_j_;
    Real pair_weight_;

    ScalarLaplaceEdge(size_t index_i, size_t index_j, Real pair_weight);
};

struct ScalarComplexHelmholtzResiduals
{
    StdVec<Complex> laplace_term_;
    StdVec<Complex> reaction_term_;
    StdVec<Complex> rhs_term_;
    StdVec<Complex> residual_;
    StdVec<Real> diagonal_scale_;
    StdVec<Real> squared_residual_;

    Real l2_norm_;
    Real mean_abs_;
    Real max_abs_;

    ScalarComplexHelmholtzResiduals();

    void resize(size_t size);
    void clear();
};

void accumulateDirectionalScalarLaplaceContribution(size_t index_i, const Complex &value_i, const Complex &value_j,
                                                    Real pair_weight, ScalarComplexHelmholtzResiduals &residuals);

void accumulateSymmetricScalarLaplaceContribution(size_t index_i, size_t index_j, const StdVec<Complex> &field,
                                                  Real pair_weight, ScalarComplexHelmholtzResiduals &residuals);

void accumulateScalarLaplaceResidualsFromEdges(const StdVec<Complex> &field, const StdVec<ScalarLaplaceEdge> &edges,
                                               ScalarComplexHelmholtzResiduals &residuals);

void finalizeScalarHelmholtzResiduals(const StdVec<Complex> &field, const StdVec<Complex> &rhs,
                                      const StdVec<Complex> &reaction_coefficient,
                                      ScalarComplexHelmholtzResiduals &residuals);

Real computeScalarHelmholtzResidualL2(const ScalarComplexHelmholtzResiduals &residuals);

} // namespace matrix_free
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_APHI_MATRIX_FREE_RESIDUALS_H
