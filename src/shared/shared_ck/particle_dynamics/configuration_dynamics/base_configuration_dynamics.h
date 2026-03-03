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
 * @file    base_configuration_dynamics.h
 * @brief   TBD.
 * @author	Xiangyu Hu
 */

#ifndef BASE_CONFIGURATION_DYNAMICS_H
#define BASE_CONFIGURATION_DYNAMICS_H

#include "algorithm_primitive.h"
#include "base_data_type_package.h"
#include "execution_policy.h"

namespace SPH
{
using namespace execution;

template <typename... T>
struct SortMethod;

class QuickSort;
template <>
struct SortMethod<SequencedPolicy>
{
    typedef QuickSort type;
};

template <>
struct SortMethod<ParallelPolicy>
{
    typedef QuickSort type;
};

template <typename... T>
struct PlusUnsignedInt;

template <>
struct PlusUnsignedInt<SequencedPolicy>
{
    typedef std::plus<UnsignedInt> type;
};

template <>
struct PlusUnsignedInt<ParallelPolicy>
{
    typedef std::plus<UnsignedInt> type;
};

template <template <typename> class ContainerType>
class UpdateSortableVariables
{
    typedef DataAssemble<UniquePtr, ContainerType> TemporaryVariables;

    struct InitializeTemporaryVariables
    {
        template <typename DataType>
        void operator()(UniquePtr<ContainerType<DataType>> &variable_ptr, UnsignedInt data_size)
        {
            variable_ptr = makeUnique<ContainerType<DataType>>("Temporary", data_size);
        };
    };

    TemporaryVariables temp_variables_;
    OperationOnDataAssemble<TemporaryVariables, InitializeTemporaryVariables> initialize_temp_variables_;

  public:
    UpdateSortableVariables(UnsignedInt data_size) : initialize_temp_variables_()
    {
        initialize_temp_variables_(temp_variables_, data_size);
    };

    template <class ExecutionPolicy, typename DataType>
    void operator()(DataContainerAddressKeeper<ContainerType<DataType>> &variables,
                    ExecutionPolicy &ex_policy, UnsignedInt start_index, UnsignedInt end_index,
                    DiscreteVariable<UnsignedInt> *dv_index_permutation)
    {
        using ContainedDataType = typename ContainerType<DataType>::ContainedDataType;
        constexpr int type_index = DataTypeIndex<DataType>::value;
        ContainedDataType *temp_data_field = std::get<type_index>(temp_variables_)->DelegatedData(ex_policy);

        UnsignedInt *index_permutation = dv_index_permutation->DelegatedData(ex_policy);

        for (size_t k = 0; k != variables.size(); ++k)
        {
            ContainedDataType *sorted_data_field = variables[k]->DelegatedData(ex_policy);
            generic_for(ex_policy, IndexRange(start_index, end_index),
                        [=](size_t i)
                        { temp_data_field[i] = sorted_data_field[i]; });
            generic_for(ex_policy, IndexRange(start_index, end_index),
                        [=](size_t i)
                        { sorted_data_field[i] = temp_data_field[index_permutation[i]]; });
        }
    };
};
} // namespace SPH
#endif // BASE_CONFIGURATION_DYNAMICS_H
