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
 * @file 	base_data_type.h
 * @brief 	This is the date type definition for SPHinXsys.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef BASE_DATA_TYPE_H
#define BASE_DATA_TYPE_H

#include <algorithm>
#include <cassert>
#include <climits>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <typeinfo>
#include <vector>

#if SPHINXSYS_USE_SYCL
#include <CL/sycl.hpp>
#define SYCL_DEVICE_ONLY
#endif // SPHINXSYS_USE_SYCL

#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>

namespace SPH
{
#if SPHINXSYS_USE_FLOAT
using Real = float;
using UnsignedInt = u_int32_t;
#else
using Real = double;
using UnsignedInt = size_t;
#endif // SPHINXSYS_USE_FLOAT

/** Vector with integers. */
using Array2i = Eigen::Array<int, 2, 1>;
using Array3i = Eigen::Array<int, 3, 1>;
/** Vector with float point number.*/
using Vec2d = Eigen::Matrix<Real, 2, 1>;
using Vec3d = Eigen::Matrix<Real, 3, 1>;
/** Small, 2*2 and 3*3, matrix with float point number. */
using Mat2d = Eigen::Matrix<Real, 2, 2>;
using Mat3d = Eigen::Matrix<Real, 3, 3>;
/** Dynamic matrix*/
using MatXd = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>;

/** Unified initialize to zero for all data type. */
template <typename DataType>
struct ZeroData
{
    static inline DataType value = DataType::Zero();
};

template <>
struct ZeroData<bool>
{
    static inline bool value = false;
};
template <>
struct ZeroData<Real>
{
    static inline Real value = Real(0);
};
template <>
struct ZeroData<int>
{
    static inline int value = 0;
};

template <>
struct ZeroData<UnsignedInt>
{
    static inline UnsignedInt value = 0;
};

template <typename DataType>
struct IdentityMatrix
{
    static inline DataType value = DataType::Identity();
};

/** Type trait for data type index. */
template <typename T>
struct DataTypeIndex
{
    static constexpr int value = std::numeric_limits<int>::max();
};
template <>
struct DataTypeIndex<UnsignedInt>
{
    static constexpr int value = 0;
};
template <>
struct DataTypeIndex<int>
{
    static constexpr int value = 1;
};
template <>
struct DataTypeIndex<Real>
{
    static constexpr int value = 2;
};
template <>
struct DataTypeIndex<Vec2d>
{
    static constexpr int value = 3;
};
template <>
struct DataTypeIndex<Mat2d>
{
    static constexpr int value = 4;
};
template <>
struct DataTypeIndex<Vec3d>
{
    static constexpr int value = 5;
};
template <>
struct DataTypeIndex<Mat3d>
{
    static constexpr int value = 6;
};

/** Verbal boolean for positive and negative axis directions. */
const int xAxis = 0;
const int yAxis = 1;
const int zAxis = 2;
const bool positiveDirection = true;
const bool negativeDirection = false;
/** Constant parameters. */
constexpr Real Pi = Real(M_PI);
constexpr Real Eps = std::numeric_limits<Real>::epsilon();
constexpr Real SqrtEps = Real(1.0e-8);
constexpr Real TinyReal = Real(2.71051e-20);
constexpr Real MinReal = std::numeric_limits<Real>::min();
constexpr Real MaxReal = std::numeric_limits<Real>::max();
constexpr size_t MaxSize_t = std::numeric_limits<size_t>::max();
} // namespace SPH
#endif // BASE_DATA_TYPE_H
