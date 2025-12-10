/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for    *
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH   *
 * (smoothed particle hydrodynamics), a meshless computational method using  *
 * particle discretization.                                                  *
 *                                                                           *
 * SPHinXsys is partially funded by German Research Foundation               *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,            *
 *  HU1527/12-1 and HU1527/12-4.                                             *
 *                                                                           *
 * Portions copyright (c) 2017-2025 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file eulerian_compressible_fluid_integration.h
 * @brief Here, we define the common compressible eulerian classes for fluid dynamics.
 * @author Zhentong Wang and Xiangyu Hu
 */

#ifndef EULERIAN_COMPRESSIBLE_FLUID_INTEGRATION_H
#define EULERIAN_COMPRESSIBLE_FLUID_INTEGRATION_H

#include "base_general_dynamics.h"
#include "compressible_fluid.h"
#include "eulerian_riemann_solver.h"
#include "fluid_integration.hpp"
#include "fluid_time_step.h"
#include "viscous_dynamics.hpp"

namespace SPH
{
namespace fluid_dynamics
{
class BaseIntegrationInCompressible : public BaseIntegration<DataDelegateInner>
{
  public:
    explicit BaseIntegrationInCompressible(BaseInnerRelation &inner_relation);
    virtual ~BaseIntegrationInCompressible() {};

  protected:
    CompressibleFluid compressible_fluid_;
    Real *Vol_, *E_, *dE_dt_, *dmass_dt_;
    Vecd *mom_, *force_, *force_prior_;
};

template <class RiemannSolverType>
class EulerianCompressibleIntegration1stHalf : public BaseIntegrationInCompressible
{
  public:
    explicit EulerianCompressibleIntegration1stHalf(BaseInnerRelation &inner_relation, Real limiter_parameter = 5.0);
    virtual ~EulerianCompressibleIntegration1stHalf() {};
    RiemannSolverType riemann_solver_;
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);
};
using EulerianCompressibleIntegration1stHalfNoRiemann = EulerianCompressibleIntegration1stHalf<NoRiemannSolverInCompressibleEulerianMethod>;
using EulerianCompressibleIntegration1stHalfHLLCRiemann = EulerianCompressibleIntegration1stHalf<HLLCRiemannSolver>;
using EulerianCompressibleIntegration1stHalfHLLCWithLimiterRiemann = EulerianCompressibleIntegration1stHalf<HLLCWithLimiterRiemannSolver>;

/**
 * @class BaseIntegration2ndHalf
 * @brief  Template density relaxation scheme in HLLC Riemann solver with and without limiter
 */
template <class RiemannSolverType>
class EulerianCompressibleIntegration2ndHalf : public BaseIntegrationInCompressible
{
  public:
    explicit EulerianCompressibleIntegration2ndHalf(BaseInnerRelation &inner_relation, Real limiter_parameter = 5.0);
    virtual ~EulerianCompressibleIntegration2ndHalf() {};
    RiemannSolverType riemann_solver_;
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);
};
using EulerianCompressibleIntegration2ndHalfNoRiemann = EulerianCompressibleIntegration2ndHalf<NoRiemannSolverInCompressibleEulerianMethod>;
using EulerianCompressibleIntegration2ndHalfHLLCRiemann = EulerianCompressibleIntegration2ndHalf<HLLCRiemannSolver>;
using EulerianCompressibleIntegration2ndHalfHLLCWithLimiterRiemann = EulerianCompressibleIntegration2ndHalf<HLLCWithLimiterRiemannSolver>;

class CompressibleFluidInitialCondition : public FluidInitialCondition
{
  public:
    explicit CompressibleFluidInitialCondition(SPHBody &sph_body);

  protected:
    Vecd *mom_;
    Real *rho_, *Vol_, *mass_, *p_, *E_;
};

class EulerianCompressibleAcousticTimeStepSize : public AcousticTimeStep
{
  protected:
    Real *rho_, *p_;
    Vecd *vel_;
    Real smoothing_length_;

  public:
    explicit EulerianCompressibleAcousticTimeStepSize(SPHBody &sph_body, Real acousticCFL = 0.6);
    virtual ~EulerianCompressibleAcousticTimeStepSize() {};

    Real reduce(size_t index_i, Real dt = 0.0);
    virtual Real outputResult(Real reduced_value) override;
    CompressibleFluid compressible_fluid_;
};
} // namespace fluid_dynamics
} // namespace SPH
#endif // EULERIAN_COMPRESSIBLE_FLUID_INTEGRATION_H