#ifndef ELECTROMAGNETIC_APHI_SPH_OPERATOR_DIVERGENCE_H
#define ELECTROMAGNETIC_APHI_SPH_OPERATOR_DIVERGENCE_H

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
 * div(F) = dFx/dx + dFy/dy + dFz/dz with uncorrected pairwise gradient weights
 * (same as validated matrix-free graph / SYCL standard-divergence path).
 */
class SPHComplexStandardVectorDivergenceOperator
{
  public:
    explicit SPHComplexStandardVectorDivergenceOperator(RealBody &body, Inner<> &ck_inner_relation);
    ~SPHComplexStandardVectorDivergenceOperator();

    StdVec<Complex> computeFromComponents(const StdVec<Complex> &field_x, const StdVec<Complex> &field_y,
                                          const StdVec<Complex> &field_z);

  private:
    void ensureInitialized();
    void copyComponentsToParticles(const StdVec<Complex> &field_x, const StdVec<Complex> &field_y,
                                 const StdVec<Complex> &field_z);

    RealBody &body_;
    Inner<> &ck_inner_relation_;
    std::unique_ptr<class SPHComplexStandardVectorDivergenceOperatorDynamics> dynamics_;
    bool initialized_ = false;
};

/**
 * div(sigma F) assembled as trace of harmonic-weighted gradients of Fx, Fy, Fz
 * (harmonic mean of sigma on each edge, uncorrected dW V e weights).
 */
class SPHComplexHarmonicVectorDivergenceOperator
{
  public:
    explicit SPHComplexHarmonicVectorDivergenceOperator(RealBody &body, Inner<> &ck_inner_relation);
    ~SPHComplexHarmonicVectorDivergenceOperator();

    StdVec<Complex> computeFromComponents(const StdVec<Complex> &field_x, const StdVec<Complex> &field_y,
                                          const StdVec<Complex> &field_z, const StdVec<Real> &edge_weight_coefficient);

  private:
    void ensureInitialized();
    void copyComponentsAndEdgeWeightsToParticles(const StdVec<Complex> &field_x, const StdVec<Complex> &field_y,
                                                   const StdVec<Complex> &field_z,
                                                   const StdVec<Real> &edge_weight_coefficient);

    RealBody &body_;
    Inner<> &ck_inner_relation_;
    std::unique_ptr<class SPHComplexHarmonicVectorDivergenceOperatorDynamics> dynamics_;
    bool initialized_ = false;
};

} // namespace sph
} // namespace electromagnetics
} // namespace SPH

#include "aphi_sphinxsys/electromagnetic_aphi_sph_operator_divergence.hpp"

#endif // ELECTROMAGNETIC_APHI_SPH_OPERATOR_DIVERGENCE_H
