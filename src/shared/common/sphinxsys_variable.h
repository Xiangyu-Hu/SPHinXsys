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

#include <algorithm>
#include <iterator>

namespace SPH
{
using namespace execution;

template <typename DataType>
class SingleVariable;

template <typename DataType>
class DiscreteVariable;

template <typename DataType>
class DataView
{
  public:
    DataView(DataType *data) : data_(data) {};

    DataType &operator[](UnsignedInt index) const
    {
        return *(data_ + index);
    }

  protected:
    DataType *data_;
};

template <typename DataType>
class EntryView
{
  public:
    EntryView(DataType *data, UnsignedInt entry, UnsignedInt width)
        : data_(data), entry_(entry), width_(width) {};

    UnsignedInt Entry() const { return entry_; };
    UnsignedInt Width() const { return width_; };

    DataType &operator[](UnsignedInt index) const
    {
        return *(data_ + entry_ + index * width_);
    }

  protected:
    DataType *data_;
    UnsignedInt entry_, width_;
};

template <typename DataType>
class MultiEntryView
{
  public:
    MultiEntryView() : data_(nullptr), width_(0) {};
    MultiEntryView(DataType *data, UnsignedInt width)
        : data_(data), width_(width) {};

    void setData(DataType *data) { data_ = data; };
    void setWidth(UnsignedInt width) { width_ = width; };
    UnsignedInt Width() const { return width_; };

    DataType *operator[](UnsignedInt index) const
    {
        return data_ + index * width_;
    }

  protected:
    DataType *data_;
    UnsignedInt width_;
};

class Quantity
{
  public:
    explicit Quantity(const std::string &name) : name_(name) {};
    virtual ~Quantity() {};
    std::string Name() const { return name_; };
    void setName(const std::string &name) { name_ = name; };
    Real getScalingRef() const { return scaling_ref_; };
    void setScalingRef(Real scaling_ref) { scaling_ref_ = scaling_ref; };

  protected:
    std::string name_;
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

struct MultiEntryTag
{
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
        : Quantity(name), size_(size), width_(1), data_(new DataType[size]),
          device_only_variable_(nullptr)
    {
        fill(initialization, 0, size);
    };

    DiscreteVariable(const std::string &name, UnsignedInt size,
                     DataType initial_value = ZeroData<DataType>::value)
        : DiscreteVariable(name, size, [&](UnsignedInt index)
                           { return initial_value; }) {};

    DiscreteVariable(const std::string &name, UnsignedInt size,
                     const MultiEntryTag &tag, UnsignedInt width)
        : Quantity(name), size_(size), width_(width), data_(new DataType[size * width]),
          device_only_variable_(nullptr)
    {
        for (size_t i = 0; i < width; i++)
        {
            entry_names_.push_back(std::to_string(i));
            fill([&](UnsignedInt index) // zero initialization
                 { return ZeroData<DataType>::value; }, 0, size, i);
        }
    };

    DiscreteVariable(const std::string &name, UnsignedInt size, StdVec<std::string> entry_names)
        : DiscreteVariable(name, size, MultiEntryTag{}, entry_names.size())
    {
        entry_names_ = entry_names;
    };

    ~DiscreteVariable() { delete[] data_; };
    DataType *Data() { return data_; };
    DataType getValue(UnsignedInt index) { return data_[index]; };
    DataType getValueWithScalingRef(UnsignedInt index) const { return data_[index] * scaling_ref_; };
    UnsignedInt getSize() { return size_; }
    UnsignedInt getWidth() { return width_; }
    UnsignedInt getTotalSize() { return size_ * width_; }
    std::string getEntryName(UnsignedInt entry) { return !entry_names_.empty() ? entry_names_[entry] : ""; }

    UnsignedInt getEntryIndexByName(std::string entry_name)
    {
        auto iter = std::find(entry_names_.begin(), entry_names_.end(), entry_name);
        if (iter != entry_names_.end())
        {
            return std::distance(entry_names_.begin(), iter);
        }
        else
        {
            std::cout << "\n Error: the variable '" << this->name_
                      << "' does not have a entry named '" << entry_name << "'!" << std::endl;
            exit(1);
        }
    };

    DataType getEntryValueWithScalingRef(UnsignedInt index, UnsignedInt entry) const
    {
        return data_[index * width_ + entry] * scaling_ref_;
    };

    template <class FillFunction>
    void fill(const FillFunction &fill_function, UnsignedInt begin_index,
              UnsignedInt fill_size, UnsignedInt entry = 0)
    {
        if (begin_index + fill_size > size_)
        {
            std::cout << "\n Error: trying to fill data out of range in DiscreteVariable '"
                      << this->name_ << "'!" << std::endl;
            exit(1);
        }

        for (UnsignedInt i = begin_index; i < begin_index + fill_size; ++i)
        {
            data_[i * width_ + entry] = fill_function(i);
        }
    };

    template <class ExecutionPolicy>
    DataType *DelegatedData(const ExecutionPolicy &ex_policy) { return data_; };

    template <class PolicyType>
    DataType *DelegatedData(const DeviceExecution<PolicyType> &ex_policy)
    {
        return DelegatedOnDevice();
    };

    template <class ExecutionPolicy>
    DataView<DataType> DelegatedDataView(const ExecutionPolicy &ex_policy)
    {
        if (width_ != 1)
        {
            std::cout << "\n Error: the variable '" << this->name_
                      << "' is not a single entry variable!" << std::endl;
            exit(1);
        }

        return DataView<DataType>(DelegatedData(ex_policy));
    };

    template <class ExecutionPolicy>
    EntryView<DataType> DelegatedEntryView(const ExecutionPolicy &ex_policy, UnsignedInt entry)
    {
        if (entry >= width_)
        {
            std::cout << "\n Error: entry index out of range in variable '"
                      << this->name_ << "'!" << std::endl;
            exit(1);
        }

        return EntryView<DataType>(DelegatedData(ex_policy), entry, width_);
    };

    template <class ExecutionPolicy>
    EntryView<DataType> DelegatedEntryView(const ExecutionPolicy &ex_policy, std::string entry_name)
    {
        return DelegatedEntryView(ex_policy, getEntryIndexByName(entry_name));
    };

    template <class ExecutionPolicy>
    MultiEntryView<DataType> DelegatedMultiEntryView(const ExecutionPolicy &ex_policy)
    {
        return MultiEntryView<DataType>(DelegatedData(ex_policy), width_);
    };

    DataView<DataType> getDataView() { return DelegatedDataView(ParallelPolicy{}); };
    EntryView<DataType> getEntryView(UnsignedInt entry) { return DelegatedEntryView(ParallelPolicy{}, entry); };
    MultiEntryView<DataType> getMultiEntryView() { return DelegatedMultiEntryView(ParallelPolicy{}); };

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
    UnsignedInt size_, width_;
    StdVec<std::string> entry_names_;
    DataType *data_;
    DeviceOnlyDiscreteVariable<DataType> *device_only_variable_ = nullptr;
    friend class DeviceOnlyDiscreteVariable<DataType>;

    DataType *DelegatedOnDevice();
    bool isDataDelegated() { return device_only_variable_ != nullptr; };

    void reallocateData(UnsignedInt tentative_size)
    {
        delete[] data_;
        size_ = tentative_size + tentative_size / 4;
        data_ = new DataType[size_ * width_];
    };

    void reallocateDataOnDevice(UnsignedInt tentative_size);
    void synchronizeWithDevice();
    void synchronizeToDevice();
};

/** Generalized particle variable type*/
typedef DataContainerAddressAssemble<DiscreteVariable> DiscreteVariables;
/** Generalized particle variable type*/
typedef DataContainerAddressAssemble<SingleVariable> SingleVariables;
} // namespace SPH
#endif // SPHINXSYS_VARIABLE_H
