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
 * Portions copyright (c) 2017-2023 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	fluid_dynamics_inner.h
 * @brief 	Here, we define the algorithm classes for fluid dynamics within the body.
 * @details We consider here weakly compressible fluids.
 * 			Note that, as these are local dynamics which are combined with particle dynamics
 * 			algorithms as template, the name-hiding is used for functions in the derived classes.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef FLUID_DYNAMICS_INNER_H
#define FLUID_DYNAMICS_INNER_H

#include "base_fluid_dynamics.h"
#include "riemann_solver.h"
#include "weakly_compressible_fluid.h"

namespace SPH
{
namespace fluid_dynamics
{
/**
 * @class FluidInitialCondition
 * @brief  Set initial condition for a fluid body.
 * This is a abstract class to be override for case specific initial conditions
 */
class FluidInitialCondition : public LocalDynamics, public FluidDataSimple
{
  public:
    explicit FluidInitialCondition(SPHBody &sph_body);
    virtual ~FluidInitialCondition(){};

  protected:
    StdLargeVec<Vecd> &pos_, &vel_;
};

/**
 * @class BaseDensitySummation
 * @brief Base class for computing density by summation
 */
template <class DataDelegationType>
class BaseDensitySummation : public LocalDynamics, public DataDelegationType
{
  public:
    template <class BaseRelationType>
    explicit BaseDensitySummation(BaseRelationType &base_relation);
    virtual ~BaseDensitySummation(){};

  protected:
    StdLargeVec<Real> &rho_, rho_sum_, &mass_;
    Real rho0_, inv_sigma0_;
};

/**
 * @class BaseDensitySummationInner
 * @brief Base class for computing density by summation
 */
class BaseDensitySummationInner : public BaseDensitySummation<FluidDataInner>
{
  public:
    explicit BaseDensitySummationInner(BaseInnerRelation &inner_relation)
        : BaseDensitySummation<FluidDataInner>(inner_relation){};
    virtual ~BaseDensitySummationInner(){};
    void update(size_t index_i, Real dt = 0.0);
};

/**
 * @class DensitySummationInner
 * @brief  computing density by summation
 */
class DensitySummationInner : public BaseDensitySummationInner
{
  public:
    explicit DensitySummationInner(BaseInnerRelation &inner_relation);
    virtual ~DensitySummationInner(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Real W0_;
};

/**
 * @class DensitySummationInnerAdaptive
 * @brief  computing density by summation with variable smoothing length
 */
class DensitySummationInnerAdaptive : public BaseDensitySummationInner
{
  public:
    explicit DensitySummationInnerAdaptive(BaseInnerRelation &inner_relation);
    virtual ~DensitySummationInnerAdaptive(){};

    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    SPHAdaptation &sph_adaptation_;
    Kernel &kernel_;
    StdLargeVec<Real> &h_ratio_;
};

/**
 * @class BaseViscousAcceleration
 * @brief Base class for the viscosity force induced acceleration
 */
template <class DataDelegationType>
class BaseViscousAcceleration : public LocalDynamics, public DataDelegationType
{
  public:
    template <class BaseRelationType>
    explicit BaseViscousAcceleration(BaseRelationType &base_relation);
    virtual ~BaseViscousAcceleration(){};

  protected:
    StdLargeVec<Real> &rho_;
    StdLargeVec<Vecd> &vel_, &acc_prior_;
    Real mu_;
    Real smoothing_length_;
};

/**
 * @class ViscousAccelerationInner
 * @brief  the viscosity force induced acceleration
 */
class ViscousAccelerationInner : public BaseViscousAcceleration<FluidDataInner>
{
  public:
    explicit ViscousAccelerationInner(BaseInnerRelation &inner_relation)
        : BaseViscousAcceleration<FluidDataInner>(inner_relation){};
    virtual ~ViscousAccelerationInner(){};
    void interaction(size_t index_i, Real dt = 0.0);
};

/**
 * @class AngularConservativeViscousAccelerationInner
 * @brief the viscosity force induced acceleration, a formulation for conserving
 * angular momentum, to be tested for its practical applications.
 */
class AngularConservativeViscousAccelerationInner : public BaseViscousAcceleration<FluidDataInner>
{
  public:
    explicit AngularConservativeViscousAccelerationInner(BaseInnerRelation &inner_relation)
        : BaseViscousAcceleration<FluidDataInner>(inner_relation){};
    virtual ~AngularConservativeViscousAccelerationInner(){};
    void interaction(size_t index_i, Real dt = 0.0);
};

/**
 * @class AcousticTimeStepSize
 * @brief Computing the acoustic time step size
 */
class AcousticTimeStepSize : public LocalDynamicsReduce<Real, ReduceMax>, public FluidDataSimple
{
  public:
    explicit AcousticTimeStepSize(SPHBody &sph_body, Real acousticCFL = 0.6);
    virtual ~AcousticTimeStepSize(){};
    Real reduce(size_t index_i, Real dt = 0.0);
    virtual Real outputResult(Real reduced_value) override;

  protected:
    Fluid &fluid_;
    StdLargeVec<Real> &rho_, &p_;
    StdLargeVec<Vecd> &vel_;
    Real smoothing_length_min_;
    Real acousticCFL_;
};

/**
 * @class AdvectionTimeStepSizeForImplicitViscosity
 * @brief Computing the advection time step size when viscosity is handled implicitly
 */
class AdvectionTimeStepSizeForImplicitViscosity
    : public LocalDynamicsReduce<Real, ReduceMax>,
      public FluidDataSimple
{
  public:
    explicit AdvectionTimeStepSizeForImplicitViscosity(
        SPHBody &sph_body, Real U_ref, Real advectionCFL = 0.25);
    virtual ~AdvectionTimeStepSizeForImplicitViscosity(){};
    Real reduce(size_t index_i, Real dt = 0.0);
    virtual Real outputResult(Real reduced_value) override;

  protected:
    StdLargeVec<Vecd> &vel_;
    Real smoothing_length_min_;
    Real speed_ref_, advectionCFL_;
};

/**
 * @class AdvectionTimeStepSize
 * @brief Computing the advection time step size
 */
class AdvectionTimeStepSize : public AdvectionTimeStepSizeForImplicitViscosity
{
  public:
    explicit AdvectionTimeStepSize(SPHBody &sph_body, Real U_ref, Real advectionCFL = 0.25);
    virtual ~AdvectionTimeStepSize(){};
    Real reduce(size_t index_i, Real dt = 0.0);

  protected:
    Fluid &fluid_;
};

/**
 * @class VorticityInner
 * @brief  compute vorticity in the fluid field
 */
class VorticityInner : public LocalDynamics, public FluidDataInner
{
  public:
    explicit VorticityInner(BaseInnerRelation &inner_relation);
    virtual ~VorticityInner(){};

    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Vecd> &vel_;
    StdLargeVec<AngularVecd> vorticity_;
};

/**
 * @class BaseIntegration
 * @brief Base class for all fluid relaxation schemes
 */
template <class DataDelegationType>
class BaseIntegration : public LocalDynamics, public DataDelegationType
{
  public:
    template <class BaseRelationType>
    explicit BaseIntegration(BaseRelationType &base_relation);
    virtual ~BaseIntegration(){};

  protected:
    Fluid &fluid_;
    StdLargeVec<Real> &rho_, &p_, &drho_dt_;
    StdLargeVec<Vecd> &pos_, &vel_, &acc_, &acc_prior_;
};

/**
 * @class BaseIntegration1stHalfInner
 * @brief Template class for pressure relaxation scheme with the Riemann solver
 * as template variable
 */
template <class RiemannSolverType, class KernelCorrectionType>
class BaseIntegration1stHalfInner : public BaseIntegration<FluidDataInner>
{
  public:
    explicit BaseIntegration1stHalfInner(BaseInnerRelation &inner_relation);
    virtual ~BaseIntegration1stHalfInner(){};
    void initialization(size_t index_i, Real dt = 0.0);
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    KernelCorrectionType correction_;
    RiemannSolverType riemann_solver_;
};
using Integration1stHalfInner = BaseIntegration1stHalfInner<NoRiemannSolver, NoKernelCorrection>;
using Integration1stHalfInnerRiemann = BaseIntegration1stHalfInner<AcousticRiemannSolver, NoKernelCorrection>;
using Integration1stHalfInnerDissipativeRiemann = BaseIntegration1stHalfInner<DissipativeRiemannSolver, NoKernelCorrection>;

/**
 * @class BaseIntegration2ndHalfInner
 * @brief  Template density relaxation scheme with different Riemann solver
 */
template <class RiemannSolverType>
class BaseIntegration2ndHalfInner : public BaseIntegration<FluidDataInner>
{
  public:
    explicit BaseIntegration2ndHalfInner(BaseInnerRelation &inner_relation);
    virtual ~BaseIntegration2ndHalfInner(){};
    void initialization(size_t index_i, Real dt = 0.0);
    inline void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    RiemannSolverType riemann_solver_;
    StdLargeVec<Real> &Vol_, &mass_;
};
using Integration2ndHalfInner = BaseIntegration2ndHalfInner<NoRiemannSolver>;
/** define the mostly used density relaxation scheme using Riemann solver */
using Integration2ndHalfInnerRiemann = BaseIntegration2ndHalfInner<AcousticRiemannSolver>;
using Integration2ndHalfInnerDissipativeRiemann = BaseIntegration2ndHalfInner<DissipativeRiemannSolver>;

/**
 * @class Oldroyd_BIntegration1stHalf
 * @brief Pressure relaxation scheme with the mostly used Riemann solver.
 */
class Oldroyd_BIntegration1stHalf : public Integration1stHalfInnerDissipativeRiemann
{
  public:
    explicit Oldroyd_BIntegration1stHalf(BaseInnerRelation &inner_relation);
    virtual ~Oldroyd_BIntegration1stHalf(){};
    void initialization(size_t index_i, Real dt = 0.0);
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Matd> tau_, dtau_dt_;
};

/**
 * @class Oldroyd_BIntegration2ndHalf
 * @brief Density relaxation scheme with the mostly used Riemann solver.
 */
class Oldroyd_BIntegration2ndHalf : public Integration2ndHalfInnerDissipativeRiemann
{
  public:
    explicit Oldroyd_BIntegration2ndHalf(BaseInnerRelation &inner_relation);
    virtual ~Oldroyd_BIntegration2ndHalf(){};
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Oldroyd_B_Fluid &oldroyd_b_fluid_;
    StdLargeVec<Matd> &tau_, &dtau_dt_;
    Real mu_p_, lambda_;
};
} // namespace fluid_dynamics
} // namespace SPH
#endif // FLUID_DYNAMICS_INNER_H