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
 * @file 	time_step_initialization.h
 * @brief 	TBD
 * @author	Chi Zhang and Xiangyu Hu
 */
#ifndef TIME_STEP_INITIALIZATION_H
#define TIME_STEP_INITIALIZATION_H

#include "base_general_dynamics.h"
#include <limits>

namespace SPH
{
/**
 * @class BaseTimeStepInitialization
 * @brief base class for time step initialization.
 */
class BaseTimeStepInitialization : public LocalDynamics
{
  private:
    SharedPtrKeeper<Gravity> gravity_ptr_keeper_;

  protected:
    Gravity *gravity_;

  public:
    BaseTimeStepInitialization(SPHBody &sph_body, SharedPtr<Gravity> &gravity_ptr)
        : LocalDynamics(sph_body), gravity_(gravity_ptr_keeper_.assignPtr(gravity_ptr)){};
    virtual ~BaseTimeStepInitialization(){};
};

/**
 * @class TimeStepInitialization
 * @brief initialize a time step for a body.
 */
class TimeStepInitialization
    : public BaseTimeStepInitialization,
      public GeneralDataDelegateSimple
{
  protected:
    StdLargeVec<Vecd> &pos_, &force_prior_;
    StdLargeVec<Real> &mass_;

  public:
    explicit TimeStepInitialization(SPHBody &sph_body, SharedPtr<Gravity> gravity_ptr = makeShared<Gravity>(Vecd::Zero()));
    virtual ~TimeStepInitialization(){};
    void update(size_t index_i, Real dt = 0.0);
};
} // namespace SPH
#endif // TIME_STEP_INITIALIZATION_H
