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
 * @file 	base_device_data_type.h
 * @brief 	This is the date type definition for SPHinXsys.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef BASE_DEVICE_DATA_TYPE_H
#define BASE_DEVICE_DATA_TYPE_H

// Specialize Eigen for device
#define EIGEN_RUNTIME_NO_MALLOC
#define EIGEN_NO_MALLOC
#define EIGEN_DONT_VECTORIZE
#define EIGEN_NO_DEBUG

#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>

namespace SPH
{

#if SPHINXSYS_DEVICE_USE_FLOAT
using DeviceReal = float;
#else
using DeviceReal = double;
#endif

/** Vector with float point number.*/
using DeviceVec2d = Eigen::Matrix<DeviceReal, 2, 1>;
using DeviceVec3d = Eigen::Matrix<DeviceReal, 3, 1>;
/** Small, 2*2 and 3*3, matrix with float point number. */
using DeviceMat2d = Eigen::Matrix<DeviceReal, 2, 2>;
using DeviceMat3d = Eigen::Matrix<DeviceReal, 3, 3>;
/** AlignedBox */
using DeviceAlignedBox2d = Eigen::AlignedBox<DeviceReal, 2>;
using DeviceAlignedBox3d = Eigen::AlignedBox<DeviceReal, 3>;
/** Rotation */
using DeviceRotation2d = Eigen::Rotation2D<DeviceReal>;
using DeviceRotation3d = Eigen::AngleAxis<DeviceReal>;

/** Type trait for device data type index, corresponding its host counterpart. */
template <typename T>
struct DeviceDataTypeIndex
{
    static constexpr int value = std::numeric_limits<int>::max();
};
template <>
struct DeviceDataTypeIndex<DeviceReal>
{
    static constexpr int value = 0;
};
template <>
struct DeviceDataTypeIndex<DeviceVec2d>
{
    static constexpr int value = 1;
};
template <>
struct DeviceDataTypeIndex<DeviceVec3d>
{
    static constexpr int value = 2;
};
template <>
struct DeviceDataTypeIndex<DeviceMat2d>
{
    static constexpr int value = 3;
};
template <>
struct DeviceDataTypeIndex<DeviceMat3d>
{
    static constexpr int value = 4;
};
template <>
struct DeviceDataTypeIndex<int>
{
    static constexpr int value = 5;
};
} // namespace SPH

#endif // BASE_DEVICE_DATA_TYPE_H
