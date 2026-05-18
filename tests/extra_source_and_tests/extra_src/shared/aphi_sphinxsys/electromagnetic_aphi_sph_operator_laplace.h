#ifndef ELECTROMAGNETIC_APHI_SPH_OPERATOR_LAPLACE_H
#define ELECTROMAGNETIC_APHI_SPH_OPERATOR_LAPLACE_H

#include "inner_body_relation.h"
#include <complex>
#include <memory>

namespace SPH
{
template <typename...>
class Inner;

namespace electromagnetics
{
namespace matrix_free
{
struct ScalarComplexHelmholtzResiduals;
} // namespace matrix_free
namespace sph
{

struct SPHScalarNegativeLaplaceParameters
{
    Real pair_weight_regularization = 0.01;
    Real contact_diffusion_scale = 1.0;
};

struct SPHComplexNegativeLaplaceField
{
    StdVec<Complex> values_;
};

/**
 * Applies u -> -div(alpha grad u) using native Inner neighbor traversal (CK Interaction),
 * with the same pairwise weights as the matrix-free graph Laplace at uniform/low-contrast alpha.
 */
class SPHComplexScalarNegativeLaplaceOperator
{
  public:
    SPHComplexScalarNegativeLaplaceOperator(RealBody &body, Inner<> &ck_inner_relation,
                                            SPHScalarNegativeLaplaceParameters parameters = {});
    ~SPHComplexScalarNegativeLaplaceOperator();

    StdVec<Complex> apply(const StdVec<Complex> &field, const StdVec<Real> &diffusion_coefficient);

    /** Fills laplace_term_ and diagonal_scale_ for Helmholtz/Jacobi (same weights as graph Laplace residual). */
    void accumulateHelmholtzLaplaceResiduals(const StdVec<Complex> &field, const StdVec<Real> &diffusion_coefficient,
                                             matrix_free::ScalarComplexHelmholtzResiduals &residuals);

  private:
    void ensureInitialized();
    void ensureRelationUpdated();
    void copyFieldAndDiffusionToParticles(const StdVec<Complex> &field, const StdVec<Real> &diffusion_coefficient);

    RealBody &body_;
    SPHScalarNegativeLaplaceParameters parameters_;
    Inner<> &ck_inner_relation_;
    std::unique_ptr<class SPHComplexScalarNegativeLaplaceOperatorDynamics> dynamics_;
    bool initialized_ = false;
    bool relation_initialized_ = false;
    const Real *last_diffusion_data_ = nullptr;
    size_t last_diffusion_size_ = 0;
    bool last_diffusion_changed_ = true;
    bool diagonal_cache_valid_ = false;
    StdVec<Real> diagonal_scale_cache_;
};

} // namespace sph
} // namespace electromagnetics
} // namespace SPH

#include "aphi_sphinxsys/electromagnetic_aphi_sph_operator_laplace.hpp"

#endif // ELECTROMAGNETIC_APHI_SPH_OPERATOR_LAPLACE_H
