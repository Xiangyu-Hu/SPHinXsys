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
/**
 * @class ExternalForce
 * @brief This class define external forces.
 */
class ExternalForce
{
  public:
    ExternalForce();
    virtual ~ExternalForce(){};
    /** This function can be used for runtime control of external force. */
    virtual Vecd InducedAcceleration(const Vecd &position = Vecd::Zero()) = 0;
};

/**
 * @class Gravity
 * @brief The gravity force, derived class of External force.
 */
class Gravity
{
  protected:
    Vecd global_acceleration_;
    Vecd zero_potential_reference_;
    DeviceVecd global_acceleration_device_;
    DeviceVecd zero_potential_reference_device_;

  public:
    Gravity(Vecd gravity_vector, Vecd reference_position = Vecd::Zero());

    /** This function can be used for runtime control of external force. */
    Vecd InducedAcceleration(const Vecd &position = VecdZero<Vecd>()) const
    {
        return global_acceleration_;
    }
    Real getPotential(const Vecd &position) const
    {
        return InducedAcceleration(position).dot(zero_potential_reference_ - position);
    }

    template<class Vec,
              typename = std::enable_if_t<std::conjunction_v<std::is_same<Vec, DeviceVecd>>,
                                          std::is_same<is_device_type_different_from_host<DeviceVecd>, std::true_type>>>
    DeviceVecd InducedAcceleration(const Vec &position = VecdZero<DeviceVecd>()) const {
        return global_acceleration_device_;
    }
    template<class Vec,
              typename = std::enable_if_t<std::conjunction_v<std::is_same<Vec, DeviceVecd>>,
                                          std::is_same<is_device_type_different_from_host<DeviceVecd>, std::true_type>>>
    DeviceReal getPotential(const Vec &position) const {
        return VecdDot(InducedAcceleration(position), DeviceVecd(zero_potential_reference_device_ - position));
    }
};
} // namespace SPH
#endif // EXTERNAL_FORCE_H