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
 * @file sphinxsys_variable_array.h
 * @brief Defines a collective view of several discrete variables
 * as an array for convenient data access and management.
 * @author Xiangyu Hu
 */

#ifndef SPHINXSYS_VARIABLE_ARRAY_H
#define SPHINXSYS_VARIABLE_ARRAY_H

#include "sphinxsys_variable.h"

namespace SPH
{
template <typename DataType>
class VariableArray;

template <typename DataType>
using DataPtr = DataType *;

template <typename DataType>
class DataArray
{
  public:
    DataArray() : data_ptr_(nullptr), array_size_(0) {};
    DataArray(DataPtr<DataType> *data_ptr, size_t array_size)
        : data_ptr_(data_ptr), array_size_(array_size) {};
    size_t ArraySize() { return array_size_; };

    DataType *operator[](size_t array_index) const
    {
        return data_ptr_[array_index];
    }

  protected:
    DataPtr<DataType> *data_ptr_;
    UnsignedInt array_size_;
};

template <typename DataType>
class DeviceOnlyVariableArray : public Quantity
{
  public:
    template <class PolicyType>
    DeviceOnlyVariableArray(const DeviceExecution<PolicyType> &ex_policy,
                            VariableArray<DataType> *host_variable_array);
    ~DeviceOnlyVariableArray();
    DataPtr<DataType> *DeviceOnlyDataPtr() { return device_only_data_ptr_; };

  protected:
    DataPtr<DataType> *device_only_data_ptr_;
};

template <typename DataType>
class VariableArray : public Quantity
{
    UniquePtrKeeper<Quantity> device_only_variable_array_keeper_;

  public:
    VariableArray(StdVec<DiscreteVariable<DataType> *> variables)
        : Quantity("VariableArray"), variables_(variables),
          array_size_(variables.size())
    {
        data_ptr_ = new DataPtr<DataType>[variables.size()];
        for (size_t i = 0; i != variables.size(); ++i)
        {
            data_ptr_[i] = variables[i]->Data();
        }
    };

    ~VariableArray() { delete[] data_ptr_; };
    StdVec<DiscreteVariable<DataType> *> getVariables() { return variables_; };
    size_t getArraySize() { return array_size_; }

    template <class ExecutionPolicy>
    DataArray<DataType> DelegatedDataArray(const ExecutionPolicy &ex_policy)
    {
        return DataArray<DataType>(data_ptr_, array_size_);
    };

    template <class PolicyType>
    DataArray<DataType> DelegatedDataArray(const DeviceExecution<PolicyType> &ex_policy)
    {
        return DataArray<DataType>(DelegatedOnDevice<PolicyType>(), array_size_);
    };

  protected:
    StdVec<DiscreteVariable<DataType> *> variables_;
    UnsignedInt array_size_;
    DataPtr<DataType> *data_ptr_ = nullptr;
    DeviceOnlyVariableArray<DataType> *device_only_variable_array_ = nullptr;
    friend class DeviceOnlyVariableArray<DataType>;

    template <class PolicyType>
    DataPtr<DataType> *DelegatedOnDevice();
    bool isDataArrayDelegated() { return device_only_variable_array_ != nullptr; };
};

typedef DataAssemble<TypeAlias, DataArray> VariableDataArrayAssemble;
typedef DataAssemble<UniquePtr, VariableArray> VariableArrayAssemble;

struct VariableArrayAssembleInitialization
{
    template <typename DataType>
    void operator()(const StdVec<DiscreteVariable<DataType> *> &variables,
                    UniquePtr<VariableArray<DataType>> &variable_array_ptr)
    {
        variable_array_ptr = std::make_unique<VariableArray<DataType>>(variables);
    }
};

struct VariableDataArrayAssembleInitialization
{
    template <typename DataType, class ExecutionPolicy>
    void operator()(const UniquePtr<VariableArray<DataType>> &variable_array_ptr,
                    DataArray<DataType> &variable_data_array,
                    const ExecutionPolicy &ex_policy)
    {
        variable_data_array = variable_array_ptr->DelegatedDataArray(ex_policy);
    }
};
} // namespace SPH
#endif // SPHINXSYS_VARIABLE_ARRAY_H
