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
#include "vector_functions.h"

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

inline SimTK::Transform EigenToSimTK(const Transform3d &transform)
{
    Rotation3d aa = transform.getRotation();
    Vec3d translation = transform.getTranslation();
    SimTK::Rotation R(aa.angle(), SimTK::UnitVec3(aa.axis()[0], aa.axis()[1], aa.axis()[2]));
    SimTK::Vec3 t(translation(0), translation(1), translation(2));
    return SimTK::Transform(R, t);
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
} // namespace SPH
#endif // SIMTK_WRAPPER_H
