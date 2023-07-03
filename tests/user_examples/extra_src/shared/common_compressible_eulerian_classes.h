/* -------------------------------------------------------------------------*
 *								SPHinXsys									*
 * -------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle*
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
 * physical accurate simulation and aims to model coupled industrial dynamic*
 * systems including fluid, solid, multi-body dynamics and beyond with SPH	*
 * (smoothed particle hydrodynamics), a meshless computational method using	*
 * particle discretization.													*
 *																			*
 * SPHinXsys is partially funded by German Research Foundation				*
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,			*
 *  HU1527/12-1 and HU1527/12-4													*
 *                                                                          *
 * Portions copyright (c) 2017-2023 Technical University of Munich and		*
 * the authors' affiliations.												*
 *                                                                          *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may  *
 * not use this file except in compliance with the License. You may obtain a*
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.       *
 *                                                                          *
 * ------------------------------------------------------------------------*/
/**
 * @file 	common_compressible_eulerian_classes.h
 * @brief 	Here, we define the common compressible eulerian classes for fluid dynamics.
 * @author	Zhentong Wang and Xiangyu Hu
 */

#ifndef COMMON_COMPRESSIBLE_EULERIAN_CLASSES_H
#define COMMON_COMPRESSIBLE_EULERIAN_CLASSES_H

#include "compressible_fluid.h"
#include "fluid_body.h"
#include "fluid_dynamics_inner.h"
#include "general_dynamics.h"
#include "riemann_solver.h"

namespace SPH
{
/**
 * @class EulerianCompressibleTimeStepInitialization
 * @brief initialize a time step for a body.
 * including initialize particle acceleration
 * induced by viscous, gravity and other forces,
 * set the number of ghost particles into zero.
 */
class EulerianCompressibleTimeStepInitialization : public TimeStepInitialization
{
  public:
    EulerianCompressibleTimeStepInitialization(SPHBody &sph_body, SharedPtr<Gravity> gravity_ptr = makeShared<Gravity>(Vecd::Zero()));
    virtual ~EulerianCompressibleTimeStepInitialization(){};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Real> &rho_;
    StdLargeVec<Vecd> &pos_, &vel_;
    StdLargeVec<Vecd> &dmom_dt_prior_;
    StdLargeVec<Real> &dE_dt_prior_;
};

/**
 * @class EulerianAcousticTimeStepSize
 * @brief Computing the acoustic time step size
 */
class EulerianCompressibleAcousticTimeStepSize : public fluid_dynamics::AcousticTimeStepSize
{
  protected:
    StdLargeVec<Real> &rho_, &p_;
    StdLargeVec<Vecd> &vel_;
    Real smoothing_length_;

  public:
    explicit EulerianCompressibleAcousticTimeStepSize(SPHBody &sph_body);
    virtual ~EulerianCompressibleAcousticTimeStepSize(){};

    Real reduce(size_t index_i, Real dt = 0.0);
    virtual Real outputResult(Real reduced_value) override;
    CompressibleFluid compressible_fluid_;
};

//----------------------------------------------------------------------
//	Remann Solver classes.
//----------------------------------------------------------------------
/**
 * @struct CompressibleFluidState
 * @brief  Struct for stored states of Riemann solver in compressible flow.
 */
struct CompressibleFluidState : FluidState
{
    Real &E_;
    CompressibleFluidState(Real &rho, Vecd &vel, Real &p, Real &E)
        : FluidState(rho, vel, p), E_(E){};
};
struct CompressibleFluidStarState : FluidStarState
{
    Real rho_;
    Real E_;
    CompressibleFluidStarState(Real rho, Vecd vel, Real p, Real E)
        : FluidStarState(vel, p), rho_(rho), E_(E){};
};

/**
 * @struct NoRiemannSolverInCompressibleEulerianMethod
 * @brief  NO RiemannSolver for weakly-compressible flow in Eulerian method for compressible flow.
 */
class NoRiemannSolverInCompressobleEulerianMethod
{
    CompressibleFluid &compressible_fluid_i_, &compressible_fluid_j_;

  public:
    NoRiemannSolverInCompressobleEulerianMethod(CompressibleFluid &fluid_i, CompressibleFluid &fluid_j);
    CompressibleFluidStarState getInterfaceState(const CompressibleFluidState &state_i, const CompressibleFluidState &state_j, const Vecd &e_ij);
};

/**
 * @struct HLLCRiemannSolver
 * @brief  HLLC Riemann solver.
 */
class HLLCRiemannSolver
{
    CompressibleFluid &compressible_fluid_i_, &compressible_fluid_j_;

  public:
    HLLCRiemannSolver(CompressibleFluid &compressible_fluid_i, CompressibleFluid &compressible_fluid_j, Real limiter_parameter = 0.0);
    CompressibleFluidStarState getInterfaceState(const CompressibleFluidState &state_i, const CompressibleFluidState &state_j, const Vecd &e_ij);
};
/**
 * @struct HLLCWithLimiterRiemannSolver
 * @brief  HLLC Riemann solver with dissipation limiter.
 */
class HLLCWithLimiterRiemannSolver
{
    CompressibleFluid &compressible_fluid_i_, &compressible_fluid_j_;
    Real limiter_parameter_;

  public:
    HLLCWithLimiterRiemannSolver(CompressibleFluid &compressible_fluid_i, CompressibleFluid &compressible_fluid_j, Real limiter_parameter = 5.0);
    CompressibleFluidStarState getInterfaceState(const CompressibleFluidState &state_i, const CompressibleFluidState &state_j, const Vecd &e_ij);
};

/**
 * @class EulerianViscousAccelerationInner
 * @brief  the viscosity force induced acceleration in Eulerian method
 */
class EulerianCompressibleViscousAccelerationInner : public fluid_dynamics::ViscousAccelerationInner
{
  public:
    explicit EulerianCompressibleViscousAccelerationInner(BaseInnerRelation &inner_relation);
    virtual ~EulerianCompressibleViscousAccelerationInner(){};
    void interaction(size_t index_i, Real dt = 0.0);
    StdLargeVec<Real> &dE_dt_prior_;
    StdLargeVec<Vecd> &dmom_dt_prior_;
};

//----------------------------------------------------------------------
//	Relaxation definition
//----------------------------------------------------------------------
/**
 * @class BaseIntegrationInCompressible
 * @brief Pure abstract base class for all fluid relaxation schemes in compressible flows
 */
class BaseIntegrationInCompressible : public fluid_dynamics::BaseIntegration
{
  public:
    explicit BaseIntegrationInCompressible(BaseInnerRelation &inner_relation);
    virtual ~BaseIntegrationInCompressible(){};

  protected:
    CompressibleFluid compressible_fluid_;
    StdLargeVec<Real> &Vol_, &E_, &dE_dt_, &dE_dt_prior_;
    StdLargeVec<Vecd> &mom_, &dmom_dt_, &dmom_dt_prior_;
};

/**
 * @class BaseIntegration1stHalf
 * @brief Template class for pressure relaxation scheme with the Riemann solver
 * as template variable
 */
template <class RiemannSolverType>
class BaseIntegration1stHalf : public BaseIntegrationInCompressible
{
  public:
    explicit BaseIntegration1stHalf(BaseInnerRelation &inner_relation, Real limiter_parameter = 5.0);
    virtual ~BaseIntegration1stHalf(){};
    RiemannSolverType riemann_solver_;
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);
};
using Integration1stHalf = BaseIntegration1stHalf<NoRiemannSolverInCompressobleEulerianMethod>;
using Integration1stHalfHLLCRiemann = BaseIntegration1stHalf<HLLCRiemannSolver>;
using Integration1stHalfHLLCWithLimiterRiemann = BaseIntegration1stHalf<HLLCWithLimiterRiemannSolver>;

/**
 * @class BaseIntegration2ndHalf
 * @brief  Template density relaxation scheme in HLLC Riemann solver with and without limiter
 */
template <class RiemannSolverType>
class BaseIntegration2ndHalf : public BaseIntegrationInCompressible
{
  public:
    explicit BaseIntegration2ndHalf(BaseInnerRelation &inner_relation, Real limiter_parameter = 5.0);
    virtual ~BaseIntegration2ndHalf(){};
    RiemannSolverType riemann_solver_;
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);
};
using Integration2ndHalf = BaseIntegration2ndHalf<NoRiemannSolverInCompressobleEulerianMethod>;
using Integration2ndHalfHLLCRiemann = BaseIntegration2ndHalf<HLLCRiemannSolver>;
using Integration2ndHalfHLLCWithLimiterRiemann = BaseIntegration2ndHalf<HLLCWithLimiterRiemannSolver>;

} // namespace SPH
#endif // COMMON_COMPRESSIBLE_EULERIAN_CLASSES_H