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
 * @file sphinxsys_variable.h
 * @brief Here gives classes for the singular and discrete variables used in simulation.
 * @details These discrete variables are those discretized in spaces and time.
 * @author Xiangyu Hu
 */

#ifndef SPHINXSYS_VARIABLE_H
#define SPHINXSYS_VARIABLE_H

#include "base_data_type_package.h"
#include "execution_policy.h"
#include "ownership.h"

namespace SPH
{
using namespace execution;

template <typename DataType>
class SingularVariable;

template <typename DataType>
class DiscreteVariable;

class Entity
{
  public:
    explicit Entity(const std::string &name) : name_(name) {};
    ~Entity() {};
    std::string Name() const { return name_; };

  protected:
    const std::string name_;
};

template <typename DataType>
class DeviceSharedSingularVariable : public Entity
{
  public:
    DeviceSharedSingularVariable(SingularVariable<DataType> *host_variable);
    ~DeviceSharedSingularVariable();

  protected:
    DataType *device_shared_data_;
};

template <typename DataType>
class SingularVariable : public Entity
{
    UniquePtrKeeper<Entity> device_shared_singular_variable_keeper_;

  public:
    SingularVariable(const std::string &name, const DataType &value)
        : Entity(name), data_(new DataType(value)), delegated_(data_) {};

    template <typename... Args>
    SingularVariable(const std::string &name, Args &&...args)
        : Entity(name), data_(new DataType(std::forward<Args>(args)...)),
          delegated_(data_){};

    ~SingularVariable() { delete data_; };

    DataType *Data() { return delegated_; };
    void setValue(const DataType &value) { *delegated_ = value; };
    DataType getValue() const { return *delegated_; };
    void incrementValue(const DataType &value) { *delegated_ += value; };

    template <class ExecutionPolicy>
    DataType *DelegatedData(const ExecutionPolicy &ex_policy) { return delegated_; };

    template <class PolicyType>
    DataType *DelegatedData(const DeviceExecution<PolicyType> &ex_policy)
    {
        return DelegatedOnDevice();
    };

    DataType *DelegatedOnDevice()
    {
        if (!isDataDelegated())
        {
            device_shared_singular_variable_keeper_
                .createPtr<DeviceSharedSingularVariable<DataType>>(this);
        }
        return delegated_;
    };
    bool isDataDelegated() { return data_ != delegated_; };
    void setDelegateData(DataType *new_delegated) { delegated_ = new_delegated; };

  protected:
    DataType *data_;
    DataType *delegated_;
};

template <typename DataType>
class DeviceOnlyDiscreteVariable : public Entity
{
  public:
    DeviceOnlyDiscreteVariable(DiscreteVariable<DataType> *host_variable);
    ~DeviceOnlyDiscreteVariable();
    void reallocateData(DiscreteVariable<DataType> *host_variable);

  protected:
    DataType *device_only_data_field_;
};

template <typename DataType>
class DiscreteVariable : public Entity
{
    UniquePtrKeeper<Entity> device_only_variable_keeper_;

  public:
    typedef DataType ContainedDataType;
    template <class InitializationFunction>
    DiscreteVariable(const std::string &name, size_t data_size,
                     const InitializationFunction &initialization)
        : Entity(name), data_size_(data_size), data_field_(new DataType[data_size]),
          device_only_variable_(nullptr), device_data_field_(nullptr)
    {
        for (size_t i = 0; i < data_size; ++i)
        {
            data_field_[i] = initialization(i);
        }
    };

    DiscreteVariable(const std::string &name, size_t data_size,
                     DataType initial_value = ZeroData<DataType>::value)
        : DiscreteVariable(name, data_size, [&](UnsignedInt index)
                           { return initial_value; }) {};

    DiscreteVariable(const std::string &name, size_t data_size,
                     DiscreteVariable<DataType> *origin_variable)
        : DiscreteVariable(name, SMAX(origin_variable->getDataSize(), data_size),
                           [&](UnsignedInt index)
                           { return origin_variable->getValue(index); }) {};

    ~DiscreteVariable() { delete[] data_field_; };
    DataType *Data() { return data_field_; };
    void setValue(size_t index, const DataType &value) { data_field_[index] = value; };
    DataType getValue(size_t index) { return data_field_[index]; };

    template <class ExecutionPolicy>
    DataType *DelegatedData(const ExecutionPolicy &ex_policy) { return data_field_; };
    template <class PolicyType>
    DataType *DelegatedData(const DeviceExecution<PolicyType> &ex_policy) { return DelegatedOnDevice(); };
    DataType *DelegatedOnDevice();
    bool isDataDelegated() { return device_data_field_ != nullptr; };
    size_t getDataSize() { return data_size_; }
    void setDeviceData(DataType *data_field) { device_data_field_ = data_field; };

    template <class ExecutionPolicy>
    void reallocateData(const ExecutionPolicy &ex_policy, size_t tentative_size)
    {
        if (data_size_ < tentative_size)
        {
            reallocateData(tentative_size);
        }
    };

    void reallocateData(const ParallelDevicePolicy &par_device, size_t tentative_size)
    {
        if (data_size_ < tentative_size)
        {
            reallocateDataOnDevice(tentative_size);
        }
    };

    template <class ExecutionPolicy>
    void prepareForOutput(const ExecutionPolicy &ex_policy) {};
    void prepareForOutput(const ParallelDevicePolicy &ex_policy) { synchronizeWithDevice(); };

    template <class ExecutionPolicy>
    void finalizeLoadIn(const ExecutionPolicy &ex_policy) {};
    void finalizeLoadIn(const ParallelDevicePolicy &ex_policy) { synchronizeToDevice(); };

  private:
    size_t data_size_;
    DataType *data_field_;
    DeviceOnlyDiscreteVariable<DataType> *device_only_variable_;
    DataType *device_data_field_;

    void reallocateData(size_t tentative_size)
    {
        delete[] data_field_;
        data_size_ = tentative_size + tentative_size / 4;
        data_field_ = new DataType[data_size_];
    };

    void reallocateDataOnDevice(size_t tentative_size);
    void synchronizeWithDevice();
    void synchronizeToDevice();
};

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
        for (size_t i = 0; i != variables.size(); ++i)
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
        for (size_t i = 0; i != variables.size(); ++i)
        {
            variables[i]->finalizeLoadIn(ex_policy);
        }
    };
};
} // namespace SPH
#endif // SPHINXSYS_VARIABLE_H
