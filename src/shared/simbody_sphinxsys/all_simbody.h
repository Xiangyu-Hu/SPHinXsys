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
 * @file    all_simbody.h
 * @brief   headers for SimBody engine.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef ALL_SIMBODY_H
#define ALL_SIMBODY_H

#include "simbody_middle.h"
#include "state_engine.h"
#include "type_wrapper.h"
#include "xml_engine.h"

#include "base_data_type.h"

namespace SPH
{
template <>
struct ZeroData<SimTK::SpatialVec>
{
    static inline SimTK::SpatialVec value = SimTK::SpatialVec(SimTKVec3(0), SimTKVec3(0));
};

template <>
struct ZeroData<SimTKVec3>
{
    static inline SimTKVec3 value = SimTKVec3(0);
};

struct SimbodyState
{
    Vec3d initial_origin_location_;
    Vec3d origin_location_, origin_velocity_, origin_acceleration_;
    Vec3d angular_velocity_, angular_acceleration_;
    Mat3d rotation_;

    SimbodyState()
        : initial_origin_location_(ZeroData<Vec3d>::value),
          origin_location_(ZeroData<Vec3d>::value),
          origin_velocity_(ZeroData<Vec3d>::value),
          origin_acceleration_(ZeroData<Vec3d>::value),
          angular_velocity_(ZeroData<Vec3d>::value),
          angular_acceleration_(ZeroData<Vec3d>::value),
          rotation_(Mat3d::Identity()) {}

    // implemented according to the Simbody API function with the same name
    void findStationLocationVelocityAndAccelerationInGround(
        const Vec3d &initial_location,
        const Vec3d &initial_normal,
        Vec3d &locationOnGround,
        Vec3d &velocityInGround,
        Vec3d &accelerationInGround,
        Vec3d &normalInGround)
    {
        Vec3d temp_location = rotation_ * (initial_location - initial_origin_location_);
        locationOnGround = origin_location_ + temp_location;

        Vec3d temp_velocity = angular_velocity_.cross(temp_location);
        velocityInGround = origin_velocity_ + temp_velocity;
        accelerationInGround = origin_acceleration_ +
                               angular_acceleration_.cross(temp_location) +
                               angular_velocity_.cross(temp_velocity);
        normalInGround = rotation_ * initial_normal;
    };
};

template <>
struct ZeroData<SimbodyState>
{
    static inline SimbodyState value = SimbodyState();
};
} // namespace SPH
#endif // ALL_SIMBODY_H