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
class VariableArrayView
{
  public:
    VariableArrayView() : multi_entry_view_(nullptr), array_size_(0) {};
    VariableArrayView(MultiEntryView<DataType> *multi_entry_view, size_t array_size)
        : multi_entry_view_(multi_entry_view), array_size_(array_size) {};
    size_t ArraySize() { return array_size_; };

    MultiEntryView<DataType> &operator[](size_t array_index) const
    {
        return multi_entry_view_[array_index];
    }

  protected:
    MultiEntryView<DataType> *multi_entry_view_;
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
    MultiEntryView<DataType> *DeviceOnlyMultiEntryView() { return device_only_multi_entry_view_; };

  protected:
    MultiEntryView<DataType> *device_only_multi_entry_view_;
};

template <typename DataType>
class VariableArray : public Quantity
{
    UniquePtrKeeper<Quantity> device_only_variable_array_keeper_;

  public:
    VariableArray(StdVec<DiscreteVariable<DataType> *> variables)
        : Quantity("VariableArray"), variables_(variables),
          array_size_(variables.size()),
          multi_entry_view_(static_cast<MultiEntryView<DataType> *>(
              std::malloc(variables.size() * sizeof(MultiEntryView<DataType>))))
    {
        for (size_t i = 0; i != variables.size(); ++i)
        {
            multi_entry_view_[i] = MultiEntryView<DataType>(
                variables[i]->Data(), variables[i]->getWidth());
        }
    };

    ~VariableArray() { free(multi_entry_view_); };
    StdVec<DiscreteVariable<DataType> *> getVariables() { return variables_; };
    size_t getArraySize() { return array_size_; }
    MultiEntryView<DataType> *getArrayData() { return multi_entry_view_; };

    template <class ExecutionPolicy>
    VariableArrayView<DataType> DelegatedVariableArrayView(const ExecutionPolicy &ex_policy)
    {
        return VariableArrayView<DataType>(multi_entry_view_, array_size_);
    };

    template <class PolicyType>
    VariableArrayView<DataType> DelegatedVariableArrayView(const DeviceExecution<PolicyType> &ex_policy)
    {
        return VariableArrayView<DataType>(DelegatedOnDevice<PolicyType>(), array_size_);
    };

  protected:
    StdVec<DiscreteVariable<DataType> *> variables_;
    UnsignedInt array_size_;
    MultiEntryView<DataType> *multi_entry_view_ = nullptr;
    DeviceOnlyVariableArray<DataType> *device_only_variable_array_ = nullptr;
    friend class DeviceOnlyVariableArray<DataType>;

    template <class PolicyType>
    MultiEntryView<DataType> *DelegatedOnDevice();
    bool isVariableArrayViewDelegated() { return device_only_variable_array_ != nullptr; };
};

typedef DataAssemble<TypeAlias, VariableArrayView> VariableArrayViewAssemble;
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

struct VariableArrayViewAssembleInitialization
{
    template <typename DataType, class ExecutionPolicy>
    void operator()(const UniquePtr<VariableArray<DataType>> &variable_array_ptr,
                    VariableArrayView<DataType> &variable_array_view,
                    const ExecutionPolicy &ex_policy)
    {
        variable_array_view = variable_array_ptr->DelegatedVariableArrayView(ex_policy);
    }
};
} // namespace SPH
#endif // SPHINXSYS_VARIABLE_ARRAY_H
