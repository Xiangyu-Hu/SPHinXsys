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
 * @file sphinxsys_variable.h
 * @brief Here gives classes for the constants and variables used in simulation.
 * @details These variables are those discretized in spaces and time.
 * @author Xiangyu Hu
 */

#ifndef SPHINXSYS_VARIABLE_H
#define SPHINXSYS_VARIABLE_H

#include "base_data_package.h"
#include "execution_policy.h"
#include "ownership.h"

namespace SPH
{
using namespace execution;

template <typename DataType>
class SingularVariable;

template <typename DataType>
class DiscreteVariable;

class BaseVariable
{
  public:
    explicit BaseVariable(const std::string &name)
        : name_(name){};
    virtual ~BaseVariable(){};
    std::string Name() const { return name_; };

  protected:
    const std::string name_;
};

template <typename DataType>
class DeviceSharedSingularVariable : public BaseVariable

{
  public:
    DeviceSharedSingularVariable(SingularVariable<DataType> *host_variable);
    virtual ~DeviceSharedSingularVariable();

  protected:
    DataType *device_shared_value_;
};
template <typename DataType>
class SingularVariable : public BaseVariable
{
    UniquePtrKeeper<BaseVariable> device_shared_singular_variable_keeper_;

  public:
    SingularVariable(const std::string &name, const DataType &value)
        : BaseVariable(name), value_(new DataType(value)), delegated_(value_){};
    virtual ~SingularVariable() { delete value_; };

    DataType *ValueAddress() { return delegated_; };

    void setValue(const DataType &value) { *delegated_ = value; };
    DataType getValue() { return *delegated_; };

    void incrementValue(const DataType &value) { *delegated_ += value; };

    template <class ExecutionPolicy>
    DataType *DelegatedData(const ExecutionPolicy &ex_policy) { return delegated_; };

    DataType *DelegatedData(const ParallelDevicePolicy &par_device)
    {
        if (!isValueDelegated())
        {
            device_shared_singular_variable_keeper_
                .createPtr<DeviceSharedSingularVariable<DataType>>(this);
        }
        return delegated_;
    };
    bool isValueDelegated() { return value_ != delegated_; };
    void setDelegateValueAddress(DataType *new_delegated) { delegated_ = new_delegated; };

  protected:
    DataType *value_;
    DataType *delegated_;
};

template <typename DataType>
class DeviceOnlyDiscreteVariable : public BaseVariable
{
  public:
    DeviceOnlyDiscreteVariable(DiscreteVariable<DataType> *host_variable);
    virtual ~DeviceOnlyDiscreteVariable();
    void reallocateDataField(DiscreteVariable<DataType> *host_variable);

  protected:
    DataType *device_only_data_field_;
};

template <typename DataType>
class DiscreteVariable : public BaseVariable
{
    UniquePtrKeeper<BaseVariable> device_only_variable_keeper_;

  public:
    DiscreteVariable(const std::string &name, size_t data_size)
        : BaseVariable(name), data_size_(data_size),
          data_field_(nullptr), device_only_variable_(nullptr),
          device_data_field_(nullptr)
    {
        data_field_ = new DataType[data_size];
    };
    virtual ~DiscreteVariable() { delete[] data_field_; };
    DataType *DataField() { return data_field_; };

    template <class ExecutionPolicy>
    DataType *DelegatedDataField(const ExecutionPolicy &ex_policy) { return data_field_; };
    DataType *DelegatedDataField(const ParallelDevicePolicy &par_device);

    bool existDeviceDataField() { return device_data_field_ != nullptr; };
    size_t getDataFieldSize() { return data_size_; }
    void setDeviceDataField(DataType *data_field) { device_data_field_ = data_field; };

    template <class ExecutionPolicy>
    void reallocateDataField(const ExecutionPolicy &ex_policy, size_t tentative_size)
    {
        if (data_size_ < tentative_size)
        {
            reallocateDataField(tentative_size);
        }
    };

    void reallocateDataField(const ParallelDevicePolicy &par_device, size_t tentative_size);
    void synchronizeWithDevice();
    void synchronizeToDevice();

  private:
    size_t data_size_;
    DataType *data_field_;
    DeviceOnlyDiscreteVariable<DataType> *device_only_variable_;
    DataType *device_data_field_;

    void reallocateDataField(size_t tentative_size)
    {
        delete[] data_field_;
        data_size_ = tentative_size + tentative_size / 4;
        data_field_ = new DataType[data_size_];
    };
};

template <typename DataType>
class MeshVariable : public BaseVariable
{
  public:
    using PackageData = PackageDataMatrix<DataType, 4>;
    MeshVariable(const std::string &name, size_t data_size)
        : BaseVariable(name), data_field_(nullptr){};
    virtual ~MeshVariable() { delete[] data_field_; };

    // void setDataField(PackageData* mesh_data){ data_field_ = mesh_data; };
    PackageData *DataField() { return data_field_; };
    void allocateAllMeshVariableData(const size_t size)
    {
        data_field_ = new PackageData[size];
    }

  private:
    PackageData *data_field_;
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
} // namespace SPH
#endif // SPHINXSYS_VARIABLE_H
