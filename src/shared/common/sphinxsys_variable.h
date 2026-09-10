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
#include <array>
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
    UniquePtrsKeeper<Quantity> device_shared_singular_variable_keeper_;

  public:
    SingleVariable(const std::string &name, const DataType &value)
        : Quantity(name), data_(new DataType(value)) { delegated_.fill(data_); };

    template <typename... Args>
    SingleVariable(const std::string &name, Args &&...args)
        : Quantity(name), data_(new DataType(std::forward<Args>(args)...)) { delegated_.fill(data_); };

    ~SingleVariable()
    {
        for (DataType *replica : host_replica_)
        {
            delete replica;
        }
        delete data_;
    };
    //----------------------------------------------------------------------
    // The delegate is resolved per device. In a multi-device run each device holds
    // its own replica, since quantities such as the number of local particles
    // differ between subdomains. When a single device is used, currentSubdomainID()
    // is always 0 and the behavior is that of a single delegate.
    //----------------------------------------------------------------------
    DataType *Data() { return delegated_[currentSubdomainID()]; };
    DataType *Data(int device_id) { return delegated_[device_id]; };
    void setValue(const DataType &value) { *Data() = value; };
    void setValue(int device_id, const DataType &value) { *Data(device_id) = value; };
    DataType getValue() const { return *delegated_[currentSubdomainID()]; };
    DataType getValue(int device_id) const { return *delegated_[device_id]; };
    DataType getValueWithScalingRef() const { return *delegated_[currentSubdomainID()] * scaling_ref_; };
    void incrementValue(const DataType &value) { *Data() += value; };

    template <class ExecutionPolicy>
    DataType *DelegatedData(const ExecutionPolicy &ex_policy) { return delegated_[currentSubdomainID()]; };

    DataType *DelegatedData(const SYCLDevicePolicy &ex_policy)
    {
        return DelegatedOnDevice();
    };

    /** Host side decomposition: each subdomain owns its own value, since counters such
     *  as the number of local particles differ between subdomains. */
    template <class PolicyType>
    DataType *DelegatedData(const DecomposedExecution<PolicyType> &ex_policy)
    {
        return DelegatedOnHostSubdomain();
    };
    /** The device replica of the current subdomain. A non-template overload, so that it
     *  beats the host template above, which is an exact match for MultiDevicePolicy too. */
    DataType *DelegatedData(const MultiDevicePolicy &ex_policy)
    {
        return DelegatedOnDevice();
    };

  protected:
    DataType *data_;
    /** One delegate per subdomain, all initially aliasing the host data. */
    std::array<DataType *, MaxSubdomains> delegated_;
    /** Host replicas owned by this variable. The device replicas are owned by
     *  DeviceSharedSingleVariable instead, hence the two separate arrays. */
    std::array<DataType *, MaxSubdomains> host_replica_{};
    friend class DeviceSharedSingleVariable<DataType>;

    DataType *DelegatedOnDevice()
    {
        if (!isDataDelegated())
        {
            device_shared_singular_variable_keeper_
                .template createPtr<DeviceSharedSingleVariable<DataType>>(this);
        }
        return delegated_[currentSubdomainID()];
    };
    /** Mirrors DelegatedOnDevice(), with the replica in ordinary host memory. The
     *  replicas are separate allocations rather than aliases of data_, so that the
     *  host path reproduces the aliasing rules of the device path exactly. */
    DataType *DelegatedOnHostSubdomain()
    {
        const int subdomain_id = currentSubdomainID();
        if (!isDataDelegated())
        {
            host_replica_[subdomain_id] = new DataType(*data_);
            delegated_[subdomain_id] = host_replica_[subdomain_id];
        }
        return delegated_[subdomain_id];
    };

    bool isDataDelegated() { return data_ != delegated_[currentSubdomainID()]; };
    bool isDataDelegated(int subdomain_id) { return data_ != delegated_[subdomain_id]; };
    void setDelegateData(DataType *new_delegated) { delegated_[currentSubdomainID()] = new_delegated; };
};

/**
 * @class HostOnlyDiscreteVariable
 * @brief Host memory replica of a discrete variable, one per subdomain.
 * @details The host counterpart of DeviceOnlyDiscreteVariable, deliberately given the
 *          same shape and the same lifetime rules. Under a host side decomposition the
 *          replicas hold disjoint particle sets, exactly as the device replicas do, and
 *          the variable's own data_ array is only used as staging for scatter, gather
 *          and I/O. Keeping the two classes symmetric is what makes a bug found on the
 *          host path the same bug as on the device path.
 */
template <typename DataType>
class HostOnlyDiscreteVariable : public Quantity
{
  public:
    explicit HostOnlyDiscreteVariable(DiscreteVariable<DataType> *host_variable)
        : Quantity(host_variable->Name()), host_only_data_(nullptr),
          total_size_(host_variable->getTotalSize())
    {
        host_only_data_ = new DataType[total_size_];
        std::copy(host_variable->Data(), host_variable->Data() + total_size_, host_only_data_);
    };
    ~HostOnlyDiscreteVariable() { delete[] host_only_data_; };

    /** Grow to the host variable's new size, keeping the existing contents. A growth
     *  triggered by one subdomain reallocates every replica, and the other subdomains
     *  may already have written theirs in the same step (for instance the neighbor
     *  lists built one subdomain after another), so those must survive. */
    void reallocateData(DiscreteVariable<DataType> *host_variable)
    {
        const UnsignedInt new_total_size = host_variable->getTotalSize();
        DataType *new_data = new DataType[new_total_size];
        std::copy(host_only_data_, host_only_data_ + std::min(total_size_, new_total_size), new_data);
        delete[] host_only_data_;
        host_only_data_ = new_data;
        total_size_ = new_total_size;
    };
    DataType *HostOnlyDataField() { return host_only_data_; };

  protected:
    DataType *host_only_data_;
    UnsignedInt total_size_;
};

template <typename DataType>
class DeviceOnlyDiscreteVariable : public Quantity
{
  public:
    DeviceOnlyDiscreteVariable(DiscreteVariable<DataType> *host_variable);
    ~DeviceOnlyDiscreteVariable();
    /** Grow to the host variable's new size, keeping the existing contents; see the
     *  host replica for why. */
    void reallocateData(DiscreteVariable<DataType> *host_variable);
    DataType *DeviceOnlyDataField() { return device_only_data_; };

  protected:
    DataType *device_only_data_;
    UnsignedInt total_size_;
};

struct MultiEntryTag
{
};

template <typename DataType>
class DiscreteVariable : public Quantity
{
    UniquePtrsKeeper<Quantity> subdomain_replica_keeper_;

  public:
    typedef DataType ContainedDataType;
    template <class InitializationFunction>
    DiscreteVariable(const std::string &name, UnsignedInt size,
                     const InitializationFunction &initialization)
        : Quantity(name), size_(size), width_(1), data_(new DataType[size])
    {
        fill(initialization, 0, size);
    };

    DiscreteVariable(const std::string &name, UnsignedInt size,
                     DataType initial_value = ZeroData<DataType>::value)
        : DiscreteVariable(name, size, [&](UnsignedInt index)
                           { return initial_value; }) {};

    DiscreteVariable(const std::string &name, UnsignedInt size,
                     const MultiEntryTag &tag, UnsignedInt width)
        : Quantity(name), size_(size), width_(width), data_(new DataType[size * width])
    {
        for (UnsignedInt i = 0; i < width; i++)
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
    
    DataType *DelegatedData(const SYCLDevicePolicy &ex_policy)
    {
        return DelegatedOnDevice();
    };

    /** Host side decomposition: the replica of the subdomain bound to this thread. */
    template <class PolicyType>
    DataType *DelegatedData(const DecomposedExecution<PolicyType> &ex_policy)
    {
        return DelegatedOnHostSubdomain();
    };
    /** The device replica of the current subdomain; non-template so that it wins over the
     *  host template above for MultiDevicePolicy. */
    DataType *DelegatedData(const MultiDevicePolicy &ex_policy)
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

    void reallocateData(const SYCLDevicePolicy &sycl_device, UnsignedInt tentative_size)
    {
        if (size_ < tentative_size)
        {
            reallocateDataOnDevice(tentative_size);
        }
    };

    void reallocateData(const MultiDevicePolicy &ex_policy, UnsignedInt tentative_size)
    {
        reallocateData(SYCLDevicePolicy{}, tentative_size);
    };

    template <class PolicyType>
    void reallocateData(const DecomposedExecution<PolicyType> &ex_policy, UnsignedInt tentative_size)
    {
        if (size_ < tentative_size)
        {
            reallocateData(tentative_size);
            for (int subdomain_id = 0; subdomain_id < numberOfSubdomains(); ++subdomain_id)
            {
                if (host_only_variable_[subdomain_id] != nullptr)
                {
                    host_only_variable_[subdomain_id]->reallocateData(this);
                }
            }
        }
    };

    /** Stage the replica of the given subdomain into, or out of, the host data array.
     *  The device path does the same over PCIe; here it is a copy, which is precisely
     *  why scatter and gather are cheap to debug on this path. */
    void synchronizeWithHostSubdomain(int subdomain_id)
    {
        if (host_only_variable_[subdomain_id] != nullptr)
        {
            DataType *replica = host_only_variable_[subdomain_id]->HostOnlyDataField();
            std::copy(replica, replica + getTotalSize(), data_);
        }
    };

    void synchronizeToHostSubdomain(int subdomain_id)
    {
        if (host_only_variable_[subdomain_id] != nullptr)
        {
            std::copy(data_, data_ + getTotalSize(),
                      host_only_variable_[subdomain_id]->HostOnlyDataField());
        }
    };

    template <class ExecutionPolicy>
    void prepareForOutput(const ExecutionPolicy &ex_policy) {};
    void prepareForOutput(const SYCLDevicePolicy &ex_policy) { synchronizeWithDevice(); };
    template <class PolicyType>
    void prepareForOutput(const DecomposedExecution<PolicyType> &ex_policy) { prepareForOutput(PolicyType{}); };

    template <class ExecutionPolicy>
    void finalizeLoadIn(const ExecutionPolicy &ex_policy) {};
    void finalizeLoadIn(const SYCLDevicePolicy &ex_policy) { synchronizeToDevice(); };
    template <class PolicyType>
    void finalizeLoadIn(const DecomposedExecution<PolicyType> &ex_policy) { finalizeLoadIn(PolicyType{}); };

  private:
    UnsignedInt size_, width_;
    StdVec<std::string> entry_names_;
    DataType *data_;
    /** One device-only allocation per device. In a multi-device run the replicas hold
     *  disjoint particle sets (the subdomain owned by that device plus its halo), not
     *  copies of one global array; the host data_ is only used as I/O staging. */
    std::array<DeviceOnlyDiscreteVariable<DataType> *, MaxSubdomains> device_only_variable_{};
    /** The host side counterpart, used by the host decomposed policies. */
    std::array<HostOnlyDiscreteVariable<DataType> *, MaxSubdomains> host_only_variable_{};
    friend class DeviceOnlyDiscreteVariable<DataType>;
    friend class HostOnlyDiscreteVariable<DataType>;

    DataType *DelegatedOnDevice();
    bool isDataDelegated() { return device_only_variable_[currentSubdomainID()] != nullptr; };
    bool isDataDelegated(int subdomain_id) { return device_only_variable_[subdomain_id] != nullptr; };

    DataType *DelegatedOnHostSubdomain()
    {
        const int subdomain_id = currentSubdomainID();
        if (host_only_variable_[subdomain_id] == nullptr)
        {
            std::lock_guard<std::mutex> lock(execution::replicaCreationMutex());
            if (host_only_variable_[subdomain_id] == nullptr)
            {
                host_only_variable_[subdomain_id] =
                    subdomain_replica_keeper_
                        .template createPtr<HostOnlyDiscreteVariable<DataType>>(this);
            }
        }
        return host_only_variable_[subdomain_id]->HostOnlyDataField();
    };

    void reallocateData(UnsignedInt tentative_size)
    {
        delete[] data_;
        size_ = tentative_size + tentative_size / 4;
        data_ = new DataType[size_ * width_];
    };

    void reallocateDataOnDevice(UnsignedInt tentative_size);

  public:
    /** Host staging of the replica held by the current device. In a multi-device run
     *  these only make sense inside a SubdomainScope, or through the gather routines of
     *  the domain decomposition, which assemble the global field for I/O. */
    void synchronizeWithDevice();
    void synchronizeToDevice();
};

/** Generalized particle variable type*/
typedef DataContainerAddressAssemble<DiscreteVariable> DiscreteVariables;
/** Generalized particle variable type*/
typedef DataContainerAddressAssemble<SingleVariable> SingleVariables;
} // namespace SPH
#endif // SPHINXSYS_VARIABLE_H
