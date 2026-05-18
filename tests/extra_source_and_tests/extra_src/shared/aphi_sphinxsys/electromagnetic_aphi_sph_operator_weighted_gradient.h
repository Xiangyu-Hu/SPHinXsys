#ifndef ELECTROMAGNETIC_APHI_SPH_OPERATOR_WEIGHTED_GRADIENT_H
#define ELECTROMAGNETIC_APHI_SPH_OPERATOR_WEIGHTED_GRADIENT_H

#include "inner_body_relation.h"
#include <complex>
#include <memory>

namespace SPH
{
template <typename...>
class Inner;

namespace electromagnetics
{
namespace sph
{

/**
 * Harmonic-weighted scalar gradient: sum_j harmonicMean(w_i,w_j) (u_j-u_i) * dW_ij V_j e_ij.
 * Uses uncorrected kernel-gradient weights to match the validated matrix-free graph operator.
 */
class SPHComplexHarmonicWeightedGradientOperator
{
  public:
    explicit SPHComplexHarmonicWeightedGradientOperator(RealBody &body, Inner<> &ck_inner_relation);
    ~SPHComplexHarmonicWeightedGradientOperator();

    StdVec<Vec3c> computeFromField(const StdVec<Complex> &field, const StdVec<Real> &edge_weight_coefficient);

  private:
    void ensureInitialized();
    void copyFieldAndEdgeWeightsToParticles(const StdVec<Complex> &field, const StdVec<Real> &edge_weight_coefficient);

    RealBody &body_;
    Inner<> &ck_inner_relation_;
    std::unique_ptr<class SPHComplexHarmonicWeightedGradientOperatorDynamics> dynamics_;
    bool initialized_ = false;
};

} // namespace sph
} // namespace electromagnetics
} // namespace SPH

#include "aphi_sphinxsys/electromagnetic_aphi_sph_operator_weighted_gradient.hpp"

#endif // ELECTROMAGNETIC_APHI_SPH_OPERATOR_WEIGHTED_GRADIENT_H
