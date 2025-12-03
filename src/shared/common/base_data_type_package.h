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
 * Portions copyright (c) 2017-2025 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	base_data_type_package.h
 * @brief 	The collection of base data types for the library.
 * @author	Chi Zhang and Xiangyu Hu
 */
#ifndef BASE_DATA_TYPE_PACKAGE_H
#define BASE_DATA_TYPE_PACKAGE_H

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
                                KeeperType<ContainerType<Mat3d>>,
                                KeeperType<ContainerType<Vec6d>>,
                                KeeperType<ContainerType<Mat6d>>,
                                KeeperType<ContainerType<VecMatGrad2d>>,
                                KeeperType<ContainerType<VecMatGrad3d>>>;
/** Generalized data container assemble type */
template <template <typename> typename ContainerType>
using DataContainerAssemble = DataAssemble<DataContainerKeeper, ContainerType>;
/** Generalized data address container assemble type */
template <template <typename> typename ContainerType>
using DataContainerAddressAssemble = DataAssemble<DataContainerAddressKeeper, ContainerType>;
/** Generalized data container unique pointer assemble type */
template <template <typename> typename ContainerType>
using DataContainerUniquePtrAssemble = DataAssemble<DataContainerUniquePtrKeeper, ContainerType>;

// Please refer: https://www.cppstories.com/2022/tuple-iteration-basics/ for the following code
template <typename DataAssembleType, typename OperationType>
class OperationOnDataAssemble
{
    static constexpr std::size_t tuple_size_ = std::tuple_size_v<DataAssembleType>;
    OperationType operation_;

    template <std::size_t... Is, typename... OperationArgs>
    void operationSequence(DataAssembleType &data_assemble, std::index_sequence<Is...>, OperationArgs &&...operation_args)
    {
        (operation_(std::get<Is>(data_assemble), std::forward<OperationArgs>(operation_args)...), ...);
    }

  public:
    template <typename... Args>
    OperationOnDataAssemble(Args &&...args) : operation_(std::forward<Args>(args)...){};

    template <typename... OperationArgs>
    void operator()(DataAssembleType &data_assemble, OperationArgs &&...operation_args)
    {
        operationSequence(data_assemble, std::make_index_sequence<tuple_size_>{}, std::forward<OperationArgs>(operation_args)...);
    }
};

template <typename DataAssembleIn, typename DataAssembleOut, typename OperationType>
class OperationBetweenDataAssembles
{
    static constexpr std::size_t tuple_size_ = std::tuple_size_v<DataAssembleIn>;
    static constexpr std::size_t tuple_size_out_ = std::tuple_size_v<DataAssembleOut>;
    static_assert(tuple_size_ == tuple_size_out_, "The size of input and output data assembles must be the same.");
    OperationType operation_;

    template <std::size_t... Is, typename... OperationArgs>
    void operationSequence(DataAssembleIn &assemble_in, DataAssembleOut &assemble_out,
                           std::index_sequence<Is...>, OperationArgs &&...operation_args)
    {
        (operation_(std::get<Is>(assemble_in), std::get<Is>(assemble_out), std::forward<OperationArgs>(operation_args)...), ...);
    }

  public:
    template <typename... Args>
    OperationBetweenDataAssembles(Args &&...args) : operation_(std::forward<Args>(args)...){};

    template <typename... OperationArgs>
    void operator()(DataAssembleIn &assemble_in, DataAssembleOut &assemble_out, OperationArgs &&...operation_args)
    {
        operationSequence(assemble_in, assemble_out,
                          std::make_index_sequence<tuple_size_>{}, std::forward<OperationArgs>(operation_args)...);
    }
};
} // namespace SPH
#endif // BASE_DATA_TYPE_PACKAGE_H
