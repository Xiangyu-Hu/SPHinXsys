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
#include <vector>

#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>

namespace SPH
{
/**
 * Matrix<T, 2, 1>::Identity  return {1,0}
 * Matrix<T, 2, 2>::Identity  return {{1,0},
 *								      {0,1},}
 * Matrix<T, n, n>::Ones  Set all element to One.
 * Enable Vectorization using CXX_FLAGS = -Ofast -march=native
 * -m :
 * These ‘-m’ options are defined for the x86 family of computers.
 * -march=cpu-type
 * Generate instructions for the machine type cpu-type. In contrast to -mtune=cpu-type,
 * which merely tunes the generated code for the specified cpu-type, -march=cpu-type allows GCC to generate code that
 * may not run at all on processors other than the one indicated.
 * Specifying -march=cpu-type implies -mtune=cpu-type, except where noted otherwise.
 * The choices for cpu-type are:
 * ‘native’ -> This selects the CPU to generate code for at compilation time by determining the processor type of the compiling machine.
 * 			Using -march=native enables all instruction subsets supported by the local machine (hence the result might not run on different machines).
 * 			Using -mtune=native produces code optimized for the local machine under the constraints of the selected instruction set.
 */

#if SPHINXSYS_USE_FLOAT
using Real = float;
using EigMat = Eigen::MatrixXf;
#else
using Real = double;
using EigMat = Eigen::MatrixXd;
#endif

/** Vector with integers. */
using Array2i = Eigen::Array<int, 2, 1>;
using Array3i = Eigen::Array<int, 3, 1>;
/** Vector with float point number.*/
using Vec2d = Eigen::Matrix<Real, 2, 1>;
using Vec3d = Eigen::Matrix<Real, 3, 1>;
/** Small, 2*2 and 3*3, matrix with float point number. */
using Mat2d = Eigen::Matrix<Real, 2, 2>;
using Mat3d = Eigen::Matrix<Real, 3, 3>;
/** AlignedBox */
using AlignedBox2d = Eigen::AlignedBox<Real, 2>;
using AlignedBox3d = Eigen::AlignedBox<Real, 3>;
/** Rotation */
using Rotation2d = Eigen::Rotation2D<Real>;
using Rotation3d = Eigen::AngleAxis<Real>;

/** Unified initialize to zero for all data type. */
/**
 * NOTE: Eigen::Matrix<> constexpr constructor?
 * Currently, there are no constexpr constructors/methods in Eigen.
 * And implementing this would be very complicated (for any non-trivial methods),
 * e.g., because SIMD functions are not easy to handle.
 */
template <typename DataType>
struct ZeroData
{
    static inline DataType value = DataType::Zero();
};
template <>
struct ZeroData<Real>
{
    static inline Real value = 0.0;
};
template <>
struct ZeroData<int>
{
    static inline int value = 0;
};
/** Type trait for data type index. */
template <typename T>
struct DataTypeIndex
{
    static constexpr int value = std::numeric_limits<int>::max();
};
template <>
struct DataTypeIndex<Real>
{
    static constexpr int value = 0;
};
template <>
struct DataTypeIndex<Vec2d>
{
    static constexpr int value = 1;
};
template <>
struct DataTypeIndex<Vec3d>
{
    static constexpr int value = 2;
};
template <>
struct DataTypeIndex<Mat2d>
{
    static constexpr int value = 3;
};
template <>
struct DataTypeIndex<Mat3d>
{
    static constexpr int value = 4;
};
template <>
struct DataTypeIndex<int>
{
    static constexpr int value = 5;
};
/** Useful float point constants. */
constexpr size_t MaxSize_t = std::numeric_limits<size_t>::max();
constexpr Real MinRealNumber = std::numeric_limits<Real>::min();
constexpr Real MaxRealNumber = std::numeric_limits<Real>::max();
/** Verbal boolean for positive and negative axis directions. */
const int xAxis = 0;
const int yAxis = 1;
const int zAxis = 2;
const bool positiveDirection = true;
const bool negativeDirection = false;
/** Constant parameters. */
constexpr Real Pi = Real(M_PI);
constexpr Real Eps = std::numeric_limits<Real>::epsilon();
constexpr Real TinyReal = Real(2.71051e-20);
constexpr Real Infinity = std::numeric_limits<Real>::max();
constexpr Real SqrtEps = Real(1.0e-8);
} // namespace SPH

#endif // BASE_DATA_TYPE_H
