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

#include "all_body_relations.h"
#include "all_particle_dynamics.h"
#include "base_kernel.h"
#include "base_particles.hpp"
#include "fluid_body.h"
#include "riemann_solver.h"
#include "weakly_compressible_fluid.h"

namespace SPH
{
namespace fluid_dynamics
{
typedef DataDelegateSimple<BaseParticles> FluidDataSimple;
typedef DataDelegateInner<BaseParticles> FluidDataInner;

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
 * @class BaseDensitySummationInner
 * @brief Base class for computing density by summation
 */
class BaseDensitySummationInner : public LocalDynamics, public FluidDataInner
{
  public:
    explicit BaseDensitySummationInner(BaseInnerRelation &inner_relation);
    virtual ~BaseDensitySummationInner(){};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Real> &rho_, rho_sum_, &mass_;
    Real rho0_, inv_sigma0_;
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

    inline void interaction(size_t index_i, Real dt = 0.0);

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

    inline void interaction(size_t index_i, Real dt = 0.0);

  protected:
    SPHAdaptation &sph_adaptation_;
    Kernel &kernel_;
    StdLargeVec<Real> &h_ratio_;
};

/**
 * @class BaseViscousAccelerationInner
 * @brief Base class for the viscosity force induced acceleration
 */
class BaseViscousAccelerationInner : public LocalDynamics, public FluidDataInner
{
  public:
    explicit BaseViscousAccelerationInner(BaseInnerRelation &inner_relation);
    virtual ~BaseViscousAccelerationInner(){};

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
class ViscousAccelerationInner : public BaseViscousAccelerationInner
{
  public:
    explicit ViscousAccelerationInner(BaseInnerRelation &inner_relation)
        : BaseViscousAccelerationInner(inner_relation){};
    virtual ~ViscousAccelerationInner(){};

    inline void interaction(size_t index_i, Real dt = 0.0);
};

/**
 * @class AngularConservativeViscousAccelerationInner
 * @brief the viscosity force induced acceleration, a formulation for conserving
 * angular momentum, to be tested for its practical applications.
 */
class AngularConservativeViscousAccelerationInner : public BaseViscousAccelerationInner
{
  public:
    explicit AngularConservativeViscousAccelerationInner(BaseInnerRelation &inner_relation)
        : BaseViscousAccelerationInner(inner_relation){};
    virtual ~AngularConservativeViscousAccelerationInner(){};

    inline void interaction(size_t index_i, Real dt = 0.0);
};

/**
 * @class TransportVelocityCorrectionInner
 * @brief The particle positions are corrected for more uniformed distribution
 * when there is negative pressure in the flow.
 * @details Note that the default coefficient is for using the dual time criteria method:
 * Dual-criteria time stepping for weakly compressible smoothed particle hydrodynamics.
 * C Zhang, M Rezavand, X Hu - Journal of Computational Physics,
 * Volume 404, 1 March 2020, 109135.
 * If single (acoustic) time step is used, the coefficient should be decrease
 * to about 1/4 of the default value.
 */
class TransportVelocityCorrectionInner : public LocalDynamics, public FluidDataInner
{
  public:
    explicit TransportVelocityCorrectionInner(BaseInnerRelation &inner_relation, Real coefficient = 0.2);
    virtual ~TransportVelocityCorrectionInner(){};

    inline void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Vecd> &pos_;
    StdLargeVec<int> &surface_indicator_;
    Real smoothing_length_sqr_;
    const Real coefficient_;
};

/**
 * @class TransportVelocityCorrectionInner
 * @brief transport velocity correction
 */
class TransportVelocityCorrectionInnerAdaptive : public LocalDynamics, public FluidDataInner
{
  public:
    explicit TransportVelocityCorrectionInnerAdaptive(BaseInnerRelation &inner_relation, Real coefficient = 0.2);
    virtual ~TransportVelocityCorrectionInnerAdaptive(){};

    inline void interaction(size_t index_i, Real dt = 0.0);

  protected:
    SPHAdaptation &sph_adaptation_;
    StdLargeVec<Vecd> &pos_;
    StdLargeVec<int> &surface_indicator_;
    Real smoothing_length_sqr_;
    const Real coefficient_;
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

    inline void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Vecd> &vel_;
    StdLargeVec<AngularVecd> vorticity_;
};

/**
 * @class BaseIntegration
 * @brief Pure abstract base class for all fluid relaxation schemes
 */
class BaseIntegration : public LocalDynamics, public FluidDataInner
{
  public:
    explicit BaseIntegration(BaseInnerRelation &inner_relation);
    virtual ~BaseIntegration(){};

  protected:
    Fluid &fluid_;
    StdLargeVec<Real> &rho_, &p_, &drho_dt_;
    StdLargeVec<Vecd> &pos_, &vel_, &acc_, &acc_prior_;
};

/**
 * @class BaseIntegration1stHalf
 * @brief Template class for pressure relaxation scheme with the Riemann solver
 * as template variable
 */
template <class RiemannSolverType>
class BaseIntegration1stHalf : public BaseIntegration
{
  public:
    explicit BaseIntegration1stHalf(BaseInnerRelation &inner_relation);
    virtual ~BaseIntegration1stHalf(){};
    RiemannSolverType riemann_solver_;
    void initialization(size_t index_i, Real dt = 0.0);

    inline void interaction(size_t index_i, Real dt = 0.0);

    void update(size_t index_i, Real dt = 0.0);

  protected:
    virtual Vecd computeNonConservativeAcceleration(size_t index_i);
};
using Integration1stHalf = BaseIntegration1stHalf<NoRiemannSolver>;
/** define the mostly used pressure relaxation scheme using Riemann solver */
using Integration1stHalfRiemann = BaseIntegration1stHalf<AcousticRiemannSolver>;
using Integration1stHalfDissipativeRiemann = BaseIntegration1stHalf<DissipativeRiemannSolver>;

/**
 * @class BaseIntegration2ndHalf
 * @brief  Template density relaxation scheme with different Riemann solver
 */
template <class RiemannSolverType>
class BaseIntegration2ndHalf : public BaseIntegration
{
  public:
    explicit BaseIntegration2ndHalf(BaseInnerRelation &inner_relation);
    virtual ~BaseIntegration2ndHalf(){};
    RiemannSolverType riemann_solver_;
    void initialization(size_t index_i, Real dt = 0.0);

    inline void interaction(size_t index_i, Real dt = 0.0);

    void update(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Real> &Vol_, &mass_;
};
using Integration2ndHalf = BaseIntegration2ndHalf<NoRiemannSolver>;
/** define the mostly used density relaxation scheme using Riemann solver */
using Integration2ndHalfRiemann = BaseIntegration2ndHalf<AcousticRiemannSolver>;
using Integration2ndHalfDissipativeRiemann = BaseIntegration2ndHalf<DissipativeRiemannSolver>;

/**
 * @class Oldroyd_BIntegration1stHalf
 * @brief Pressure relaxation scheme with the mostly used Riemann solver.
 */
class Oldroyd_BIntegration1stHalf : public Integration1stHalfDissipativeRiemann
{
  public:
    explicit Oldroyd_BIntegration1stHalf(BaseInnerRelation &inner_relation);
    virtual ~Oldroyd_BIntegration1stHalf(){};
    void initialization(size_t index_i, Real dt = 0.0);

    inline void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Matd> tau_, dtau_dt_;
};

/**
 * @class Oldroyd_BIntegration2ndHalf
 * @brief Density relaxation scheme with the mostly used Riemann solver.
 */
class Oldroyd_BIntegration2ndHalf : public Integration2ndHalfDissipativeRiemann
{
  public:
    explicit Oldroyd_BIntegration2ndHalf(BaseInnerRelation &inner_relation);
    virtual ~Oldroyd_BIntegration2ndHalf(){};

    inline void interaction(size_t index_i, Real dt = 0.0);

    void update(size_t index_i, Real dt = 0.0);

  protected:
    Oldroyd_B_Fluid &oldroyd_b_fluid_;
    StdLargeVec<Matd> &tau_, &dtau_dt_;
    Real mu_p_, lambda_;
};
} // namespace fluid_dynamics
} // namespace SPH
#endif // FLUID_DYNAMICS_INNER_H