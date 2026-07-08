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
#include "data_type.h"
#include "execution_policy.h"

#include <functional>
#include <tuple>

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
    DataAssemble<UniquePtr, ContainerType> temp_variables_;

  public:
    template <class ExecutionPolicy, typename DataType>
    void operator()(DataContainerAddressKeeper<ContainerType<DataType>> &variables,
                    ExecutionPolicy &ex_policy, UnsignedInt start_index, UnsignedInt end_index,
                    DiscreteVariable<UnsignedInt> *dv_index_permutation)
    {
        using ContainedDataType = typename ContainerType<DataType>::ContainedDataType;
        constexpr int type_index = DataTypeIndex<DataType>::value;
        UnsignedInt *index_permutation = dv_index_permutation->DelegatedData(ex_policy);

        if (!std::get<type_index>(temp_variables_))
        {
            UnsignedInt total_size = 0;
            for (UnsignedInt k = 0; k != variables.size(); ++k)
            {
                total_size = SMAX(total_size, variables[k]->getTotalSize());
            }

            std::get<type_index>(temp_variables_) =
                makeUnique<ContainerType<DataType>>("Temporary", total_size);
        }

        for (UnsignedInt k = 0; k != variables.size(); ++k)
        {
            MultiEntryView<ContainedDataType> sorted_view =
                variables[k]->DelegatedMultiEntryView(ex_policy);

            MultiEntryView<ContainedDataType> temp_view(
                std::get<type_index>(temp_variables_)->DelegatedData(ex_policy),
                variables[k]->getWidth());

            generic_for(ex_policy, IndexRange(start_index, end_index),
                        [=](UnsignedInt i)
                        {
                            for (UnsignedInt entry = 0; entry != sorted_view.Width(); ++entry)
                            {
                                temp_view[i][entry] = sorted_view[i][entry];
                            }
                        });
            generic_for(ex_policy, IndexRange(start_index, end_index),
                        [=](UnsignedInt i)
                        {
                            for (UnsignedInt entry = 0; entry != sorted_view.Width(); ++entry)
                            {
                                sorted_view[i][entry] = temp_view[index_permutation[i]][entry];
                            }
                        });
        }
    };
};
} // namespace SPH
#endif // BASE_CONFIGURATION_DYNAMICS_H
