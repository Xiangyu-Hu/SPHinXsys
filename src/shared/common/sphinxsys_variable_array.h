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
 * @file sphinxsys_variable_array.h
 * @brief tbd
 * @author Xiangyu Hu
 */

#ifndef SPHINXSYS_VARIABLE_ARRAY_H
#define SPHINXSYS_VARIABLE_ARRAY_H

#include "sphinxsys_variable.h"

namespace SPH
{
template <typename DataType>
using DataArray = DataType *;

template <typename DataType, template <typename> class VariableType>
class VariableArray : public Entity
{
    UniquePtrKeeper<Entity> device_only_variable_array_keeper_;

  public:
    VariableArray(StdVec<VariableType<DataType> *> variables)
        : Entity("VariableArray"), variables_(variables),
          data_array_(nullptr), delegated_data_array_(nullptr)
    {
        data_array_ = new DataArray<DataType>[variables.size()];
        for (size_t i = 0; i != variables.size(); ++i)
        {
            data_array_[i] = variables[i]->Data();
        }
        delegated_data_array_ = data_array_;
    };
    ~VariableArray() { delete[] data_array_; };
    StdVec<VariableType<DataType> *> getVariables() { return variables_; };
    DataArray<DataType> *Data() { return data_array_; };

    template <class ExecutionPolicy>
    DataArray<DataType> *DelegatedDataArray(const ExecutionPolicy &ex_policy)
    {
        return data_array_;
    };

    template <class PolicyType>
    DataArray<DataType> *DelegatedOnDevice(const DeviceExecution<PolicyType> &ex_policy);

    template <class PolicyType>
    DataArray<DataType> *DelegatedDataArray(const DeviceExecution<PolicyType> &ex_policy);

    size_t getArraySize() { return variables_.size(); }

    bool isDataArrayDelegated() { return data_array_ != delegated_data_array_; };
    void setDelegateDataArray(DataArray<DataType> *data_array_)
    {
        delegated_data_array_ = data_array_;
    };

  private:
    StdVec<VariableType<DataType> *> variables_;
    DataArray<DataType> *data_array_;
    DataArray<DataType> *delegated_data_array_;
};

template <typename DataType, template <typename> class VariableType>
class DeviceOnlyVariableArray : public Entity
{
  public:
    template <class PolicyType>
    DeviceOnlyVariableArray(const DeviceExecution<PolicyType> &ex_policy,
                            VariableArray<DataType, VariableType> *host_variable_array);
    ~DeviceOnlyVariableArray();

  protected:
    DataArray<DataType> *device_only_data_array_;
};

template <typename DataType>
using DiscreteVariableArray = VariableArray<DataType, DiscreteVariable>;

template <typename DataType>
using AllocatedDataArray = DataArray<DataType> *;

template <typename AllocationType>
using VariableAllocationPair = std::pair<AllocationType, UnsignedInt>;

typedef DataAssemble<VariableAllocationPair, AllocatedDataArray> VariableDataArrays;
typedef DataAssemble<UniquePtr, DiscreteVariableArray> DiscreteVariableArrays;

struct DiscreteVariableArraysInitialization
{
    template <typename DataType>
    void operator()(const StdVec<DiscreteVariable<DataType> *> &variables,
                    UniquePtr<DiscreteVariableArray<DataType>> &variable_array_ptr)
    {
        variable_array_ptr = std::make_unique<DiscreteVariableArray<DataType>>(variables);
    }
};

struct VariableDataArraysInitialization
{
    template <typename DataType, class ExecutionPolicy>
    void operator()(const UniquePtr<DiscreteVariableArray<DataType>> &variable_array_ptr,
                    VariableAllocationPair<AllocatedDataArray<DataType>> &variable_allocation_size_pair,
                    const ExecutionPolicy &ex_policy)
    {
        variable_allocation_size_pair =
            std::make_pair(variable_array_ptr->DelegatedDataArray(ex_policy), variable_array_ptr->getArraySize());
    }
};
} // namespace SPH
#endif // SPHINXSYS_VARIABLE_ARRAY_H
