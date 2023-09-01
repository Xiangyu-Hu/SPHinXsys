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
 * @file 	base_data_package.h
 * @brief 	Base data package for the library.
 * @author	Chi Zhang and Xiangyu Hu
 */
#ifndef BASE_DATA_PACKAGE_H
#define BASE_DATA_PACKAGE_H

#include "array_allocation.h"
#include "data_type.h"
#include "large_data_containers.h"
#include "ownership.h"

#define TBB_PARALLEL true

namespace SPH
{
/**
 * @class Transform
 * @brief Wrapper for SimTK::Transform
 * Note that the rotation is around the frame (or local) origin.
 */
class Transform
{
  private:
    Matd rotation_, inv_rotation_;
    Vecd translation_;

  public:
    explicit Transform(const Rotation &rotation, const Vecd &translation = Vecd::Zero())
        : rotation_(rotation.toRotationMatrix()), inv_rotation_(rotation_.transpose()), translation_(translation){};
    explicit Transform(const Vecd &translation)
        : rotation_(Matd::Identity()), inv_rotation_(rotation_.transpose()), translation_(translation){};
    Transform() : Transform(Vecd::Zero()){};

    /** Forward rotation. */
    Vecd xformFrameVecToBase(const Vecd &origin)
    {
        return rotation_ * origin;
    };

    /** Forward transformation. Note that the rotation operation is carried out first. */
    Vecd shiftFrameStationToBase(const Vecd &origin)
    {
        return translation_ + xformFrameVecToBase(origin);
    };

    /** Inverse rotation. */
    Vecd xformBaseVecToFrame(const Vecd &target)
    {
        return inv_rotation_ * target;
    };

    /** Inverse transformation. Note that the inverse translation operation is carried out first. */
    Vecd shiftBaseStationToFrame(const Vecd &target)
    {
        return xformBaseVecToFrame(target - translation_);
    };
};

constexpr Real OneOverDimensions = 1.0 / (Real)Dimensions;

/** Generalized data container assemble type */
template <template <typename DataType> typename DataContainerType>
using DataContainerAssemble =
    std::tuple<StdVec<DataContainerType<Real>>,
               StdVec<DataContainerType<Vec2d>>,
               StdVec<DataContainerType<Vec3d>>,
               StdVec<DataContainerType<Mat2d>>,
               StdVec<DataContainerType<Mat3d>>,
               StdVec<DataContainerType<int>>>;
/** Generalized data container address assemble type */
template <template <typename DataType> typename DataContainerType>
using DataContainerAddressAssemble =
    std::tuple<StdVec<DataContainerType<Real> *>,
               StdVec<DataContainerType<Vec2d> *>,
               StdVec<DataContainerType<Vec3d> *>,
               StdVec<DataContainerType<Mat2d> *>,
               StdVec<DataContainerType<Mat3d> *>,
               StdVec<DataContainerType<int> *>>;
/** Generalized data container unique pointer assemble type */
template <template <typename DataType> typename DataContainerType>
using DataContainerUniquePtrAssemble =
    std::tuple<UniquePtrsKeeper<DataContainerType<Real>>,
               UniquePtrsKeeper<DataContainerType<Vec2d>>,
               UniquePtrsKeeper<DataContainerType<Vec3d>>,
               UniquePtrsKeeper<DataContainerType<Mat2d>>,
               UniquePtrsKeeper<DataContainerType<Mat3d>>,
               UniquePtrsKeeper<DataContainerType<int>>>;

/** a type irrelevant operation on the data assembles  */
template <template <typename VariableType> typename OperationType>
struct DataAssembleOperation
{
    OperationType<Real> scalar_operation;
    OperationType<Vec2d> vector2d_operation;
    OperationType<Vec3d> vector3d_operation;
    OperationType<Mat2d> matrix2d_operation;
    OperationType<Mat3d> matrix3d_operation;
    OperationType<int> integer_operation;

    template <typename... OperationArgs>
    void operator()(OperationArgs &&...operation_args)
    {
        scalar_operation(std::forward<OperationArgs>(operation_args)...);
        vector2d_operation(std::forward<OperationArgs>(operation_args)...);
        vector3d_operation(std::forward<OperationArgs>(operation_args)...);
        matrix2d_operation(std::forward<OperationArgs>(operation_args)...);
        matrix3d_operation(std::forward<OperationArgs>(operation_args)...);
        integer_operation(std::forward<OperationArgs>(operation_args)...);
    }
};
} // namespace SPH

#endif // BASE_DATA_PACKAGE_H
