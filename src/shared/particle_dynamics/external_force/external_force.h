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
 * @file 	external_force.h
 * @brief 	Here, we define the base external force class.
 * @details The simple derived classes, such as gravity will be defined in applications.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef EXTERNAL_FORCE_H
#define EXTERNAL_FORCE_H

#include "base_data_package.h"

namespace SPH
{

class Gravity
{
  protected:
    Vecd reference_acceleration_;
    Vecd zero_potential_reference_;

  public:
    Gravity(Vecd gravity_vector, Vecd reference_position = Vecd::Zero());
    ~Gravity(){};

    Vecd InducedAcceleration(const Vecd &position = Vecd::Zero(), Real physical_time = 0.0);
    Real getPotential(const Vecd &position);
};

class StartupAcceleration : public Gravity
{
    Real target_time_;

  public:
    StartupAcceleration(Vecd target_velocity, Real target_time);
    ~StartupAcceleration(){};

    Vecd InducedAcceleration(const Vecd &position, Real physical_time);
};

class IncreaseToFullGravity : public Gravity
{
    Real time_to_full_gravity_;

  public:
    explicit IncreaseToFullGravity(Vecd gravity_vector, Real time_to_full_gravity)
        : Gravity(gravity_vector), time_to_full_gravity_(time_to_full_gravity) {}
    Vecd InducedAcceleration(const Vecd &position, Real physical_time);
};
} // namespace SPH
#endif // EXTERNAL_FORCE_H