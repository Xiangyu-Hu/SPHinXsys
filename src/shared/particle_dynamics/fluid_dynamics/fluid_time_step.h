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
 * @class AcousticTimeStepSize
 * @brief Computing the acoustic time step size
 */
class AcousticTimeStepSize : public LocalDynamicsReduce<ReduceMax>
{
  public:
    explicit AcousticTimeStepSize(SPHBody &sph_body, Real acousticCFL = 0.6);
    virtual ~AcousticTimeStepSize(){};
    Real reduce(size_t index_i, Real dt = 0.0);
    virtual Real outputResult(Real reduced_value) override;

  protected:
    Fluid &fluid_;
    Real *rho_, *p_, *mass_;
    Vecd *vel_, *force_, *force_prior_;
    Real smoothing_length_min_;
    Real acousticCFL_;
};

/**
 * @class AdvectionTimeStepSizeForImplicitViscosity
 * @brief Computing the advection time step size when viscosity is handled implicitly
 */
class AdvectionTimeStepSizeForImplicitViscosity
    : public LocalDynamicsReduce<ReduceMax>
{
  public:
    explicit AdvectionTimeStepSizeForImplicitViscosity(
        SPHBody &sph_body, Real U_ref, Real advectionCFL = 0.25);
    virtual ~AdvectionTimeStepSizeForImplicitViscosity(){};
    Real reduce(size_t index_i, Real dt = 0.0);
    virtual Real outputResult(Real reduced_value) override;

  protected:
    Real *mass_;
    Vecd *vel_, *force_, *force_prior_;
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
} // namespace fluid_dynamics
} // namespace SPH
#endif // FLUID_TIME_STEP_H