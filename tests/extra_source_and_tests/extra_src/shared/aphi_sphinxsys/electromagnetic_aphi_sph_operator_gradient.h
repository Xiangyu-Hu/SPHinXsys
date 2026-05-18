#ifndef ELECTROMAGNETIC_APHI_SPH_OPERATOR_GRADIENT_H
#define ELECTROMAGNETIC_APHI_SPH_OPERATOR_GRADIENT_H

#include "inner_body_relation.h"
#include <complex>
#include <memory>
#include <string>

namespace SPH
{
template <typename...>
class Inner;

namespace electromagnetics
{
namespace sph
{

struct SPHComplexGradientField
{
    StdVec<Vec3c> values_;
};

struct SPHComplexScalarGradientOperatorDynamics;

/**
 * B-corrected scalar gradient via native shared_ck operators (LinearCorrectionMatrix + LinearGradient).
 * Complex input is split into real/imag particle variables; gradients are merged per component.
 */
class SPHComplexScalarGradientOperator
{
  public:
    explicit SPHComplexScalarGradientOperator(RealBody &body, Inner<> &ck_inner_relation);
    ~SPHComplexScalarGradientOperator();

    StdVec<Vec3c> computeFromField(const StdVec<Complex> &field);

  private:
    void ensureInitialized();
    void copyFieldToParticleScalars(const StdVec<Complex> &field);

    RealBody &body_;
    Inner<> &ck_inner_relation_;
    std::unique_ptr<SPHComplexScalarGradientOperatorDynamics> dynamics_;
    bool initialized_ = false;
};

} // namespace sph
} // namespace electromagnetics
} // namespace SPH

#include "aphi_sphinxsys/electromagnetic_aphi_sph_operator_gradient.hpp"

#endif // ELECTROMAGNETIC_APHI_SPH_OPERATOR_GRADIENT_H
