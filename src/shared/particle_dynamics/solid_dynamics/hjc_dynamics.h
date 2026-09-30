// SPDX-License-Identifier: Apache-2.0
#ifndef HJC_DYNAMICS_H
#define HJC_DYNAMICS_H
#include "elastic_dynamics.h"
#include "hjc_solid.h"

namespace SPH
{
namespace solid_dynamics
{
/** History-dependent Cauchy stress with the existing total Lagrangian SPH forces. */
class HJCIntegration1stHalf : public Integration1stHalf
{
  public:
    explicit HJCIntegration1stHalf(BaseInnerRelation &inner_relation);
    void initialization(size_t index_i, Real dt = 0.0);

  protected:
    HJCSolid &hjc_solid_;
};

/** Includes the current dense-EOS tangent in the longitudinal wave speed. */
class HJCAcousticTimeStep : public LocalDynamicsReduce<ReduceMin<Real>>
{
  public:
    explicit HJCAcousticTimeStep(SPHBody &body, Real acousticCFL = 0.2);
    Real reduce(size_t index_i, Real dt = 0.0);

  protected:
    Real *rho_, *mass_, *modulus_;
    Vecd *vel_, *force_, *force_prior_;
    Real smoothing_length_, acousticCFL_;
};
} // namespace solid_dynamics
} // namespace SPH
#endif
