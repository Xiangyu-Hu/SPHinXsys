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

#ifndef TIME_STEP_SIZE_MOVING_VELOCITY_H
#define TIME_STEP_SIZE_MOVING_VELOCITY_H

#include "all_body_relations.h"
#include "all_particle_dynamics.h"
#include "base_kernel.h"
#include "base_particles.hpp"
#include "fluid_body.h"
#include "weakly_compressible_fluid.h"

namespace SPH
{
namespace fluid_dynamics
{
typedef DataDelegateSimple<BaseParticles> FluidDataSimple;

/**
 * @class AcousticTimeStepSize
 * @brief Computing the acoustic time step size
 */
class AcousticTimeStepSizeMovingVelocity : public LocalDynamicsReduce<Real, ReduceMax>, public FluidDataSimple
{
  public:
    explicit AcousticTimeStepSizeMovingVelocity(SPHBody &sph_body, Real moving_velocity, Real acousticCFL = 0.6);
    virtual ~AcousticTimeStepSizeMovingVelocity(){};
    Real reduce(size_t index_i, Real dt = 0.0);
    virtual Real outputResult(Real reduced_value) override;

  protected:
    Fluid &fluid_;
    StdLargeVec<Real> &rho_, &p_;
    StdLargeVec<Vecd> &vel_;
    Real smoothing_length_min_;
    Real acousticCFL_;
    Real moving_velocity_;
};

/**
 * @class AdvectionTimeStepSize
 * @brief Computing the advection time step size
 */
class AdvectionTimeStepSizeMovingVelocity : public LocalDynamicsReduce<Real, ReduceMax>, public FluidDataSimple
{
  public:
    explicit AdvectionTimeStepSizeMovingVelocity(SPHBody &sph_body, Real U_ref, Real moving_velocity, Real advectionCFL = 0.25);
    virtual ~AdvectionTimeStepSizeMovingVelocity(){};
    Real reduce(size_t index_i, Real dt = 0.0);
    virtual Real outputResult(Real reduced_value) override;
  protected:
    Fluid &fluid_;
    Real moving_velocity_;
    StdLargeVec<Vecd> &vel_;
    Real smoothing_length_min_;
    Real speed_ref_, advectionCFL_;
};


} // namespace fluid_dynamics
} // namespace SPH
#endif // TIME_STEP_SIZE_MOVING_VELOCITY_H