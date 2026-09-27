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
 *  HU1527/12-1 and HU1527/12-4                                              *
 *                                                                           *
 * Portions copyright (c) 2017-2022 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	simtk_wrapper.h
 * @brief 	type wrapper between eigen and SimTK.
 * @author	Chi Zhang and Xiangyu Hu
 */
#ifndef SIMTK_WRAPPER_H
#define SIMTK_WRAPPER_H

#include "base_data_type.h"
#include "ownership.h"
#include "simbody_middle.h"

namespace SPH
{
using SimTKVec2 = SimTK::Vec2;
using SimTKVec3 = SimTK::Vec3;
using SimTKMat22 = SimTK::Mat22;
using SimTKMat33 = SimTK::Mat33;

inline SimTKVec2 EigenToSimTK(const Vec2d &eigen_vector)
{
    return SimTKVec2((double)eigen_vector[0], (double)eigen_vector[1]);
}

inline SimTKVec3 EigenToSimTK(const Vec3d &eigen_vector)
{
    return SimTKVec3((double)eigen_vector[0], (double)eigen_vector[1], (double)eigen_vector[2]);
}

inline Vec2d SimTKToEigen(const SimTKVec2 &simTK_vector)
{
    return Vec2d((Real)simTK_vector[0], (Real)simTK_vector[1]);
}

inline Vec3d SimTKToEigen(const SimTKVec3 &simTK_vector)
{
    return Vec3d((Real)simTK_vector[0], (Real)simTK_vector[1], (Real)simTK_vector[2]);
}

inline SimTKMat22 EigenToSimTK(const Mat2d &eigen_matrix)
{
    return SimTKMat22((double)eigen_matrix(0, 0), (double)eigen_matrix(0, 1),
                      (double)eigen_matrix(1, 0), (double)eigen_matrix(1, 1));
}

inline SimTKMat33 EigenToSimTK(const Mat3d &eigen_matrix)
{
    return SimTKMat33((double)eigen_matrix(0, 0), (double)eigen_matrix(0, 1), (double)eigen_matrix(0, 2),
                      (double)eigen_matrix(1, 0), (double)eigen_matrix(1, 1), (double)eigen_matrix(1, 2),
                      (double)eigen_matrix(2, 0), (double)eigen_matrix(2, 1), (double)eigen_matrix(2, 2));
}

inline Mat2d SimTKToEigen(const SimTKMat22 &simTK_matrix)
{
    return Mat2d{
        {(Real)simTK_matrix(0, 0), (Real)simTK_matrix(0, 1)},
        {(Real)simTK_matrix(1, 0), (Real)simTK_matrix(1, 1)}};
}

inline Mat3d SimTKToEigen(const SimTKMat33 &simTK_matrix)
{
    return Mat3d{
        {(Real)simTK_matrix(0, 0), (Real)simTK_matrix(0, 1), (Real)simTK_matrix(0, 2)},
        {(Real)simTK_matrix(1, 0), (Real)simTK_matrix(1, 1), (Real)simTK_matrix(1, 2)},
        {(Real)simTK_matrix(2, 0), (Real)simTK_matrix(2, 1), (Real)simTK_matrix(2, 2)}};
}

template <>
struct ZeroData<SimTK::SpatialVec>
{
    static inline const SimTK::SpatialVec value = SimTK::SpatialVec(SimTKVec3(0), SimTKVec3(0));
};

template <>
struct ZeroData<SimTKVec3>
{
    static inline const SimTKVec3 value = SimTKVec3(0);
};

struct SimbodyState
{
    Vec3d initial_origin_location_;
    Vec3d origin_location_, origin_velocity_, origin_acceleration_;
    Vec3d angular_velocity_, angular_acceleration_;
    Mat3d rotation_;

    SimbodyState();
    SimbodyState(
        const SimTKVec3 &sim_tk_initial_origin_location,
        SimTK::MobilizedBody &mobod, const SimTK::State &state);

    // implemented according to the Simbody API function with the same name
    void findStationLocationVelocityAndAccelerationInGround(
        const Vec3d &initial_location, const Vec3d &initial_normal,
        Vec3d &locationOnGround, Vec3d &velocityInGround,
        Vec3d &accelerationInGround, Vec3d &normalInGround);

    void printSimbodyState();
};

template <>
struct ZeroData<SimbodyState>
{
    static inline const SimbodyState value = SimbodyState();
};

class SimbodyStateEngine;

class SimbodySystem
{
    UniquePtrKeeper<SimbodyStateEngine> state_engine_keeper_;
    UniquePtrsKeeper<SimTK::Body::Rigid> rigid_bodies_keeper_;
    UniquePtrsKeeper<SimTK::MobilizedBody> mobilized_bodies_keeper_;

  public:
    SimbodySystem();
    SimTK::MultibodySystem &getMultibodySystem() { return MBsystem_; };
    SimTK::SimbodyMatterSubsystem &getSimbodyMatterSubsystem() { return simbody_matter_; };
    SimTK::RungeKuttaMersonIntegrator &getSimbodyIntegrator() { return integ_; };
    SimbodyStateEngine &getSimbodyStateEngine();

    SimTK::Body::Rigid &createRigidBody(
        const std::string &name, const SimTK::MassProperties &mass_properties);
    SimTK::Body::Rigid &getRigidBody(const std::string &name);

    template <class MobilizedBodyType, class ParentBodyType>
    MobilizedBodyType &createMobilizedBody(
        const std::string &name, ParentBodyType &parent_mobod,
        const SimTK::Transform &X_PF, const SimTK::Body::Rigid &body,
        const SimTK::Transform &X_BM)
    {
        MobilizedBodyType *mobilized_body =
            mobilized_bodies_keeper_.createPtr<MobilizedBodyType>(
                parent_mobod, X_PF, body, X_BM);
        mobilized_bodies_.push_back(std::make_pair(name, mobilized_body));
        return *mobilized_body;
    };

    template <class MobilizedBodyType>
    MobilizedBodyType &getMobilizedBody(const std::string &name)
    {
        for (size_t i = 0; i < mobilized_bodies_.size(); ++i)
        {
            if (mobilized_bodies_[i].first == name)
            {
                return *DynamicCast<MobilizedBodyType>(this, mobilized_bodies_[i].second);
            }
        }
        std::cerr << "\n Error: Mobilized body " << name << " not found!" << std::endl;
        throw std::runtime_error("Error: Mobilized body not found!");
    };

  protected:
    SimTK::MultibodySystem MBsystem_;
    SimTK::SimbodyMatterSubsystem simbody_matter_;
    SimTK::RungeKuttaMersonIntegrator integ_;
    StdVec<std::pair<std::string, SimTK::Body::Rigid *>> rigid_bodies_;
    StdVec<std::pair<std::string, SimTK::MobilizedBody *>> mobilized_bodies_;
};
} // namespace SPH
#endif // SIMTK_WRAPPER_H
