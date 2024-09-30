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
#include "common_functors.h"
#include "data_type.h"
#include "large_data_containers.h"
#include "ownership.h"
#include "vector_functions.h"

#define TBB_PARALLEL true

namespace SPH
{
constexpr Real OneOverDimensions = 1.0 / (Real)Dimensions;
constexpr int lastAxis = Dimensions - 1;
/** Generalized data container keeper */
template <typename ContainedDataType>
using DataContainerKeeper = StdVec<ContainedDataType>;
/** Generalized data address container keeper */
template <typename ContainedDataType>
using DataContainerAddressKeeper = StdVec<ContainedDataType *>;
/** Generalized data container unique pointer keeper */
template <typename ContainedDataType>
using DataContainerUniquePtrKeeper = UniquePtrsKeeper<ContainedDataType>;
/** Generalized data container allocated data keeper */
template <typename DataType>
using AllocatedData = DataType *;

template <template <typename> typename KeeperType, template <typename> typename ContainerType>
using DataAssemble = std::tuple<KeeperType<ContainerType<UnsignedInt>>,
                                KeeperType<ContainerType<int>>,
                                KeeperType<ContainerType<Real>>,
                                KeeperType<ContainerType<Vec2d>>,
                                KeeperType<ContainerType<Mat2d>>,
                                KeeperType<ContainerType<Vec3d>>,
                                KeeperType<ContainerType<Mat3d>>>;
/** Generalized data container assemble type */
template <template <typename> typename ContainerType>
using DataContainerAssemble = DataAssemble<DataContainerKeeper, ContainerType>;
/** Generalized data address container assemble type */
template <template <typename> typename ContainerType>
using DataContainerAddressAssemble = DataAssemble<DataContainerAddressKeeper, ContainerType>;
/** Generalized data container unique pointer assemble type */
template <template <typename> typename ContainerType>
using DataContainerUniquePtrAssemble = DataAssemble<DataContainerUniquePtrKeeper, ContainerType>;

/** a type irrelevant operation on the data assembles  */
template <template <typename> typename OperationType>
class DataAssembleOperation
{
    OperationType<int> integer_operation;
    OperationType<Real> scalar_operation;
    OperationType<Vec2d> vector2d_operation;
    OperationType<Mat2d> matrix2d_operation;
    OperationType<Vec3d> vector3d_operation;
    OperationType<Mat3d> matrix3d_operation;

  public:
    template <typename... Args>
    DataAssembleOperation(Args &&...args)
        : integer_operation(std::forward<Args>(args)...),
          scalar_operation(std::forward<Args>(args)...),
          vector2d_operation(std::forward<Args>(args)...),
          matrix2d_operation(std::forward<Args>(args)...),
          vector3d_operation(std::forward<Args>(args)...),
          matrix3d_operation(std::forward<Args>(args)...){};
    template <typename... OperationArgs>
    void operator()(OperationArgs &&...operation_args)
    {
        integer_operation(std::forward<OperationArgs>(operation_args)...);
        scalar_operation(std::forward<OperationArgs>(operation_args)...);
        vector2d_operation(std::forward<OperationArgs>(operation_args)...);
        matrix2d_operation(std::forward<OperationArgs>(operation_args)...);
        vector3d_operation(std::forward<OperationArgs>(operation_args)...);
        matrix3d_operation(std::forward<OperationArgs>(operation_args)...);
    }
};

// Please refer: https://www.cppstories.com/2022/tuple-iteration-basics/ for the following code
template <typename DataAssembleType, typename OperationType>
class OperationOnDataAssemble
{
    static constexpr std::size_t tuple_size_ = std::tuple_size_v<DataAssembleType>;
    DataAssembleType &data_assemble_;
    OperationType operation_;

    template <std::size_t... Is, typename... OperationArgs>
    void operationSequence(std::index_sequence<Is...>, OperationArgs &&...operation_args)
    {
        (operation_(std::get<Is>(data_assemble_), std::forward<OperationArgs>(operation_args)...), ...);
    }

  public:
    template <typename... Args>
    OperationOnDataAssemble(DataAssembleType &data_assemble, Args &&...args)
        : data_assemble_(data_assemble), operation_(std::forward<Args>(args)...){};

    template <typename... OperationArgs>
    void operator()(OperationArgs &&...operation_args)
    {
        operationSequence(std::make_index_sequence<tuple_size_>{}, std::forward<OperationArgs>(operation_args)...);
    }
};
} // namespace SPH
#endif // BASE_DATA_PACKAGE_H
