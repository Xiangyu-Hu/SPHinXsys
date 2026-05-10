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
 * @file variable_functions.h
 * @brief Functions for the management of variables in SPHinXsys.
 * @details TBD.
 * @author Xiangyu Hu
 */

#ifndef VARIABLE_FUNCTIONS_H
#define VARIABLE_FUNCTIONS_H

#include "sphinxsys_variable.h"

namespace SPH
{
template <typename DataType, template <typename VariableDataType> class VariableType>
VariableType<DataType> *findVariableByName(DataContainerAddressAssemble<VariableType> &assemble,
                                           const std::string &name)
{
    constexpr int type_index = DataTypeIndex<DataType>::value;
    auto &variables = std::get<type_index>(assemble);
    auto result = std::find_if(variables.begin(), variables.end(),
                               [&](auto &variable) -> bool
                               { return variable->Name() == name; });

    return result != variables.end() ? *result : nullptr;
};

template <typename DataType, template <typename VariableDataType> class VariableType, typename... Args>
VariableType<DataType> *addVariableToAssemble(DataContainerAddressAssemble<VariableType> &assemble,
                                              DataContainerUniquePtrAssemble<VariableType> &ptr_assemble, Args &&...args)
{
    constexpr int type_index = DataTypeIndex<DataType>::value;
    UniquePtrsKeeper<VariableType<DataType>> &variable_ptrs = std::get<type_index>(ptr_assemble);
    VariableType<DataType> *new_variable =
        variable_ptrs.template createPtr<VariableType<DataType>>(std::forward<Args>(args)...);
    std::get<type_index>(assemble).push_back(new_variable);
    return new_variable;
};

template <template <typename> typename ContainerType, typename DataType, typename... Args>
ContainerType<DataType> *registerVariable(DataContainerAddressAssemble<ContainerType> &all_variable_set,
                                          DataContainerUniquePtrAssemble<ContainerType> &all_variable_ptrs_,
                                          const std::string &name, Args &&...args)
{
    ContainerType<DataType> *variable = findVariableByName<DataType, ContainerType>(all_variable_set, name);
    if (variable == nullptr)
    {
        return addVariableToAssemble<DataType, ContainerType>(
            all_variable_set, all_variable_ptrs_, name, std::forward<Args>(args)...);
    }
    return variable;
};

template <template <typename> typename ContainerType, typename DataType>
ContainerType<DataType> *addVariableToList(DataContainerAddressAssemble<ContainerType> &variable_set,
                                           ContainerType<DataType> *variable)
{
    ContainerType<DataType> *listed_variable = findVariableByName<DataType, ContainerType>(variable_set, variable->Name());
    if (listed_variable == nullptr)
    {
        constexpr int type_index = DataTypeIndex<DataType>::value;
        std::get<type_index>(variable_set).push_back(variable);
        return variable;
    }
    return nullptr; // no need to add
};

template <template <typename> class ContainerType>
struct PrepareVariablesToWrite
{
    template <class ExecutionPolicy, typename DataType>
    void operator()(DataContainerAddressKeeper<ContainerType<DataType>> &variables,
                    const ExecutionPolicy &ex_policy)
    {
        for (UnsignedInt i = 0; i != variables.size(); ++i)
        {
            variables[i]->prepareForOutput(ex_policy);
        }
    };
};

template <template <typename> class ContainerType>
struct FinalizeVariablesAfterRead
{
    template <class ExecutionPolicy, typename DataType>
    void operator()(DataContainerAddressKeeper<ContainerType<DataType>> &variables,
                    const ExecutionPolicy &ex_policy)
    {
        for (UnsignedInt i = 0; i != variables.size(); ++i)
        {
            variables[i]->finalizeLoadIn(ex_policy);
        }
    };
};
} // namespace SPH
#endif // VARIABLE_FUNCTIONS_H
