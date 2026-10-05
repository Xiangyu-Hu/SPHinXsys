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
 * @file 	fluid_time_step.h
 * @brief 	Here, we define the algorithm classes for fluid dynamics within the body.
 * @details We consider here weakly compressible fluids.
 * 			Note that, as these are local dynamics which are combined with particle dynamics
 * 			algorithms as template, the name-hiding is used for functions in the derived classes.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef FLUID_TIME_STEP_H
#define FLUID_TIME_STEP_H

#include "base_fluid_dynamics.h"

namespace SPH
{
namespace fluid_dynamics
{
/**
 * @class AcousticTimeStep
 * @brief Computing the acoustic time step size
 */
class AcousticTimeStep : public LocalDynamicsReduce<ReduceMax<Real>>
{
  public:
    explicit AcousticTimeStep(SPHBody &sph_body, Real acousticCFL = 0.6);
    virtual ~AcousticTimeStep() {};
    Real reduce(size_t index_i, Real dt = 0.0);
    virtual Real outputResult(Real reduced_value) override;

  protected:
    Fluid &fluid_;
    Real *rho_, *p_;
    Vecd *vel_;
    Real h_min_;
    Real acousticCFL_;
};
/**
 * @class AcousticTimeStepWithAcceleration
 * @brief Opt-in acoustic bound including the current fluid acceleration.
 * @details The default CPU AcousticTimeStep and AdvectionTimeStep keep their
 * existing criteria. Moving the acceleration criterion to acoustic substeps
 * is a different dual-time-stepping scheme, not an interchangeable default:
 * the outer step also sets the neighbor-configuration and density-update
 * interval. Use the WithoutAcceleration advection variants only after
 * validating that interval for the application, and evaluate this bound at
 * every acoustic substep.
 */
class AcousticTimeStepWithAcceleration : public AcousticTimeStep
{
  public:
    explicit AcousticTimeStepWithAcceleration(SPHBody &sph_body, Real acousticCFL = 0.6);
    Real reduce(size_t index_i, Real dt = 0.0);

  protected:
    Real *mass_;
    Vecd *force_, *force_prior_;
};
/**
 * @class WallAccelerationTimeStep
 * @brief Computing a time step from the relative normal acceleration at fluid-wall contacts.
 * @details This bounds the acceleration term used by the wall Riemann integration when the
 * wall is moving or accelerating relative to the fluid.
 */
class WallAccelerationTimeStep
    : public LocalDynamicsReduce<ReduceMax<Real>>, public DataDelegateContact
{
  public:
    explicit WallAccelerationTimeStep(BaseContactRelation &wall_contact_relation,
                                      Real wallCFL = 0.25);
    virtual ~WallAccelerationTimeStep() {};
    Real reduce(size_t index_i, Real dt = 0.0);
    virtual Real outputResult(Real reduced_value) override;

  protected:
    Real *mass_;
    Vecd *force_prior_;
    StdVec<Vecd *> wall_acc_ave_;
    Real h_min_;
    Real wallCFL_;
};
/**
 * @class SurfaceTensionTimeStep
 * @brief Computing the acoustic time step size considering surface tension
 */
class SurfaceTensionTimeStep : public AcousticTimeStep
{
  public:
    explicit SurfaceTensionTimeStep(SPHBody &sph_body, Real acousticCFL = 0.6);
    virtual ~SurfaceTensionTimeStep() {};
    virtual Real outputResult(Real reduced_value) override;

  protected:
    Real rho0_, &surface_tension_coeff_;
};
/**
 * @class AdvectionTimeStep
 * @brief Computing the advection time step size when viscosity is handled implicitly
 */
class AdvectionTimeStep
    : public LocalDynamicsReduce<ReduceMax<Real>>
{
  protected:
    Real *mass_;
    Vecd *vel_, *force_, *force_prior_;
    Real h_min_;
    Real speed_ref_, advectionCFL_;

  public:
    AdvectionTimeStep(SPHBody &sph_body, Real U_ref, Real advectionCFL = 0.25);
    virtual ~AdvectionTimeStep() {};
    Real reduce(size_t index_i, Real dt = 0.0);
    virtual Real outputResult(Real reduced_value) override;
};

/**
 * @class AdvectionViscousTimeStep
 * @brief Computing the advection time step size
 */
class AdvectionViscousTimeStep : public AdvectionTimeStep
{
  public:
    AdvectionViscousTimeStep(SPHBody &sph_body, Real U_ref, Real advectionCFL = 0.25);
    virtual ~AdvectionViscousTimeStep() {};
    Real reduce(size_t index_i, Real dt = 0.0);
};
/**
 * @class AdvectionTimeStepWithoutAcceleration
 * @brief Opt-in velocity-only outer step, paired with AcousticTimeStepWithAcceleration.
 * @details Neighbor configurations are held over this interval. Existing CPU
 * cases should retain AdvectionTimeStep unless this alternative is validated.
 */
class AdvectionTimeStepWithoutAcceleration : public AdvectionTimeStep
{
  public:
    using AdvectionTimeStep::AdvectionTimeStep;
    Real reduce(size_t index_i, Real dt = 0.0);
};
/**
 * @class AdvectionViscousTimeStepWithoutAcceleration
 * @brief Opt-in velocity/viscosity outer step for the acceleration-aware acoustic scheme.
 */
class AdvectionViscousTimeStepWithoutAcceleration : public AdvectionViscousTimeStep
{
  public:
    using AdvectionViscousTimeStep::AdvectionViscousTimeStep;
    Real reduce(size_t index_i, Real dt = 0.0);
};
} // namespace fluid_dynamics
} // namespace SPH
#endif // FLUID_TIME_STEP_H
