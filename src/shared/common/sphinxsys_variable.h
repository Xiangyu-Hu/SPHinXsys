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
class SingleVariable;

template <typename DataType>
class DiscreteVariable;

class Quantity
{
  public:
    explicit Quantity(const std::string &name) : name_(name) {};
    virtual ~Quantity() {};
    std::string Name() const { return name_; };
    Real getScalingRef() const { return scaling_ref_; };
    void setScalingRef(Real scaling_ref) { scaling_ref_ = scaling_ref; };

  protected:
    const std::string name_;
    Real scaling_ref_ = Real(1);
};

template <typename DataType>
class DeviceSharedSingleVariable : public Quantity
{
  public:
    DeviceSharedSingleVariable(SingleVariable<DataType> *host_variable);
    ~DeviceSharedSingleVariable();

  protected:
    DataType *device_shared_data_;
};

template <typename DataType>
class SingleVariable : public Quantity
{
    UniquePtrKeeper<Quantity> device_shared_singular_variable_keeper_;

  public:
    SingleVariable(const std::string &name, const DataType &value)
        : Quantity(name), data_(new DataType(value)), delegated_(data_) {};

    template <typename... Args>
    SingleVariable(const std::string &name, Args &&...args)
        : Quantity(name), data_(new DataType(std::forward<Args>(args)...)),
          delegated_(data_){};

    ~SingleVariable() { delete data_; };

    DataType *Data() { return delegated_; };
    void setValue(const DataType &value) { *delegated_ = value; };
    DataType getValue() const { return *delegated_; };
    DataType getValueWithScalingRef() const { return *delegated_ * scaling_ref_; };
    void incrementValue(const DataType &value) { *delegated_ += value; };

    template <class ExecutionPolicy>
    DataType *DelegatedData(const ExecutionPolicy &ex_policy) { return delegated_; };

    template <class PolicyType>
    DataType *DelegatedData(const DeviceExecution<PolicyType> &ex_policy)
    {
        return DelegatedOnDevice();
    };

  protected:
    DataType *data_;
    DataType *delegated_;
    friend class DeviceSharedSingleVariable<DataType>;

    DataType *DelegatedOnDevice()
    {
        if (!isDataDelegated())
        {
            device_shared_singular_variable_keeper_
                .createPtr<DeviceSharedSingleVariable<DataType>>(this);
        }
        return delegated_;
    };
    bool isDataDelegated() { return data_ != delegated_; };
    void setDelegateData(DataType *new_delegated) { delegated_ = new_delegated; };
};

template <typename DataType>
class DeviceOnlyDiscreteVariable : public Quantity
{
  public:
    DeviceOnlyDiscreteVariable(DiscreteVariable<DataType> *host_variable);
    ~DeviceOnlyDiscreteVariable();
    void reallocateData(DiscreteVariable<DataType> *host_variable);
    DataType *DeviceOnlyDataField() { return device_only_data_; };

  protected:
    DataType *device_only_data_;
};

template <typename DataType>
class DiscreteVariable : public Quantity
{
    UniquePtrKeeper<Quantity> device_only_variable_keeper_;

  public:
    typedef DataType ContainedDataType;
    template <class InitializationFunction>
    DiscreteVariable(const std::string &name, UnsignedInt size,
                     const InitializationFunction &initialization)
        : Quantity(name), size_(size), data_(new DataType[size]),
          device_only_variable_(nullptr)
    {
        fill(initialization, 0, size);
    };

    DiscreteVariable(const std::string &name, UnsignedInt size,
                     DataType initial_value = ZeroData<DataType>::value)
        : DiscreteVariable(name, size, [&](UnsignedInt index)
                           { return initial_value; }) {};

    ~DiscreteVariable() { delete[] data_; };
    DataType *Data() { return data_; };
    void setValue(UnsignedInt index, const DataType &value) { data_[index] = value; };
    DataType getValue(UnsignedInt index) { return data_[index]; };
    DataType getValueWithScalingRef(UnsignedInt index) const { return data_[index] * scaling_ref_; };
    UnsignedInt getSize() { return size_; }

    template <class FillFunction>
    void fill(const FillFunction &fill_function, UnsignedInt begin_index, UnsignedInt fill_size)
    {
        if (begin_index + fill_size > size_)
        {
            std::cout << "\n Error: trying to fill data out of range in DiscreteVariable '"
                      << this->name_ << "'!" << std::endl;
            exit(1);
        }

        for (UnsignedInt i = begin_index; i < begin_index + fill_size; ++i)
        {
            data_[i] = fill_function(i);
        }
    };

    template <class ExecutionPolicy>
    DataType *DelegatedData(const ExecutionPolicy &ex_policy) { return data_; };
    template <class PolicyType>
    DataType *DelegatedData(const DeviceExecution<PolicyType> &ex_policy) { return DelegatedOnDevice(); };

    template <class ExecutionPolicy>
    void reallocateData(const ExecutionPolicy &ex_policy, UnsignedInt tentative_size)
    {
        if (size_ < tentative_size)
        {
            reallocateData(tentative_size);
        }
    };

    void reallocateData(const ParallelDevicePolicy &par_device, UnsignedInt tentative_size)
    {
        if (size_ < tentative_size)
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
    UnsignedInt size_;
    DataType *data_;
    DeviceOnlyDiscreteVariable<DataType> *device_only_variable_ = nullptr;
    friend class DeviceOnlyDiscreteVariable<DataType>;

    DataType *DelegatedOnDevice();
    bool isDataDelegated() { return device_only_variable_ != nullptr; };

    void reallocateData(UnsignedInt tentative_size)
    {
        delete[] data_;
        size_ = tentative_size + tentative_size / 4;
        data_ = new DataType[size_];
    };

    void reallocateDataOnDevice(UnsignedInt tentative_size);
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

/** Generalized particle variable type*/
typedef DataContainerAddressAssemble<DiscreteVariable> DiscreteVariables;
/** Generalized particle variable type*/
typedef DataContainerAddressAssemble<SingleVariable> SingleVariables;
} // namespace SPH
#endif // SPHINXSYS_VARIABLE_H
