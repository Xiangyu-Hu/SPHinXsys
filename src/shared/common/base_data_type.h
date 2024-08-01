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

#ifndef SPHINXSYS_USE_SYCL
#include <CL/sycl.hpp>
#define SYCL_DEVICE_ONLY // Eigen matrix and vector types degraded
#endif                   // SPHINXSYS_USE_SYCL

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
/** AlignedBox */
using AlignedBox2d = Eigen::AlignedBox<Real, 2>;
using AlignedBox3d = Eigen::AlignedBox<Real, 3>;
/** Rotation */
using Rotation2d = Eigen::Rotation2D<Real>;
using Rotation3d = Eigen::AngleAxis<Real>;

/** Unified initialize to zero for all data type. */
template <typename DataType>
struct ZeroData
{
    static inline DataType value = DataType::Zero();
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

/** Bounding box for system, body, body part and shape, first: lower bound, second: upper bound. */
template <typename VecType>
class BaseBoundingBox
{
  public:
    VecType first_, second_;
    int dimension_;

    BaseBoundingBox() : first_(VecType::Zero()), second_(VecType::Zero()), dimension_(VecType::Zero().size()){};
    BaseBoundingBox(const VecType &lower_bound, const VecType &upper_bound)
        : first_(lower_bound), second_(upper_bound), dimension_(lower_bound.size()){};
    /** Check the bounding box contain. */
    bool checkContain(const VecType &point)
    {
        bool is_contain = true;
        for (int i = 0; i < dimension_; ++i)
        {
            if (point[i] < first_[i] || point[i] > second_[i])
            {
                is_contain = false;
                break;
            }
        }
        return is_contain;
    };

    VecType getBoundSize()
    {
        return second_ - first_;
    };
};
/** Operator define. */
template <class T>
bool operator==(const BaseBoundingBox<T> &bb1, const BaseBoundingBox<T> &bb2)
{
    return bb1.first_ == bb2.first_ && bb1.second_ == bb2.second_ ? true : false;
};
/** Intersection fo bounding box.*/
template <class BoundingBoxType>
BoundingBoxType getIntersectionOfBoundingBoxes(const BoundingBoxType &bb1, const BoundingBoxType &bb2)
{
    /** Check that the inputs are correct. */
    int dimension = bb1.dimension_;
    /** Get the Bounding Box of the intersection of the two meshes. */
    BoundingBoxType bb(bb1);
    /** #1 check that there is overlap, if not, exception. */
    for (int i = 0; i < dimension; ++i)
        if (bb2.first_[i] > bb1.second_[i] || bb2.second_[i] < bb1.first_[i])
            std::runtime_error("getIntersectionOfBoundingBoxes: no overlap!");
    /** #2 otherwise modify the first one to get the intersection. */
    for (int i = 0; i < dimension; ++i)
    {
        /** If the lower limit is inside change the lower limit. */
        if (bb1.first_[i] < bb2.first_[i] && bb2.first_[i] < bb1.second_[i])
            bb.first_[i] = bb2.first_[i];
        /**  If the upper limit is inside, change the upper limit. */
        if (bb1.second_[i] > bb2.second_[i] && bb2.second_[i] > bb1.first_[i])
            bb.second_[i] = bb2.second_[i];
    }
    return bb;
}

/** obtain minimum dimension of a bounding box */
template <class BoundingBoxType>
auto MinimumDimension(const BoundingBoxType &bbox)
{
    return (bbox.second_ - bbox.first_).cwiseAbs().minCoeff();
};

template <typename RotationType, typename VecType>
class BaseTransform
{
  private:
    using MatType = decltype(RotationType().toRotationMatrix());
    MatType rotation_, inv_rotation_;
    VecType translation_;

  public:
    explicit BaseTransform(const RotationType &rotation, const VecType &translation = VecType::Zero())
        : rotation_(rotation.toRotationMatrix()), inv_rotation_(rotation_.transpose()), translation_(translation){};
    explicit BaseTransform(const VecType &translation)
        : rotation_(MatType::Identity()), inv_rotation_(rotation_.transpose()), translation_(translation){};
    BaseTransform() : BaseTransform(VecType::Zero()){};

    /** Forward rotation. */
    VecType xformFrameVecToBase(const VecType &origin)
    {
        return rotation_ * origin;
    };

    /** Forward transformation. Note that the rotation operation is carried out first. */
    VecType shiftFrameStationToBase(const VecType &origin)
    {
        return translation_ + xformFrameVecToBase(origin);
    };

    /** Inverse rotation. */
    VecType xformBaseVecToFrame(const VecType &target)
    {
        return inv_rotation_ * target;
    };

    /** Inverse transformation. Note that the inverse translation operation is carried out first. */
    VecType shiftBaseStationToFrame(const VecType &target)
    {
        return xformBaseVecToFrame(target - translation_);
    };
};
} // namespace SPH

#endif // BASE_DATA_TYPE_H
