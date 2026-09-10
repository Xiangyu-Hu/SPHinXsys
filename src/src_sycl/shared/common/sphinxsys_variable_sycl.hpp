#ifndef SPHINXSYS_VARIABLE_SYCL_HPP
#define SPHINXSYS_VARIABLE_SYCL_HPP

#include "implementation_sycl.h"
#include "sphinxsys_variable.h"

namespace SPH
{
//=================================================================================================//
template <typename DataType>
DeviceSharedSingleVariable<DataType>::
    DeviceSharedSingleVariable(SingleVariable<DataType> *host_variable)
    : Quantity(host_variable->Name()),
      device_shared_data_(allocateDeviceShared<DataType>(1))
{
    copyToDevice(host_variable->Data(), device_shared_data_, 1);
    host_variable->setDelegateData(device_shared_data_);
}
//=================================================================================================//
template <typename DataType>
DeviceSharedSingleVariable<DataType>::~DeviceSharedSingleVariable()
{
    freeDeviceData(device_shared_data_);
}
//=================================================================================================//
template <typename DataType>
void DiscreteVariable<DataType>::synchronizeWithDevice()
{
    if (isDataDelegated())
    {
        copyFromDevice(data_, device_only_variable_[currentSubdomainID()]->DeviceOnlyDataField(),
                       getTotalSize());
    }
}
//=================================================================================================//
template <typename DataType>
void DiscreteVariable<DataType>::synchronizeToDevice()
{
    if (isDataDelegated())
    {
        copyToDevice(data_, device_only_variable_[currentSubdomainID()]->DeviceOnlyDataField(),
                     getTotalSize());
    }
}
//=================================================================================================//
template <typename DataType>
DeviceOnlyDiscreteVariable<DataType>::
    DeviceOnlyDiscreteVariable(DiscreteVariable<DataType> *host_variable)
    : Quantity(host_variable->Name()), device_only_data_(nullptr),
      total_size_(host_variable->getTotalSize())
{
    device_only_data_ = allocateDeviceOnly<DataType>(total_size_);
    copyToDevice(host_variable->Data(), device_only_data_, total_size_);
}
//=================================================================================================//
template <typename DataType>
DeviceOnlyDiscreteVariable<DataType>::~DeviceOnlyDiscreteVariable()
{
    freeDeviceData(device_only_data_);
}
//=================================================================================================//
template <typename DataType>
void DeviceOnlyDiscreteVariable<DataType>::
    reallocateData(DiscreteVariable<DataType> *host_variable)
{
    // Staged through the host: the old contents are copied back, the new allocation
    // made, and the contents copied in again. Rare enough (a neighbor list growth)
    // for the round trip not to matter.
    const UnsignedInt new_total_size = host_variable->getTotalSize();
    const UnsignedInt kept_size = std::min(total_size_, new_total_size);
    StdVec<DataType> kept(kept_size);
    copyFromDevice(kept.data(), device_only_data_, kept_size);
    freeDeviceData(device_only_data_);
    device_only_data_ = allocateDeviceOnly<DataType>(new_total_size);
    copyToDevice(kept.data(), device_only_data_, kept_size);
    total_size_ = new_total_size;
}
//=================================================================================================//
template <typename DataType>
DataType *DiscreteVariable<DataType>::DelegatedOnDevice()
{
    const int device_id = currentSubdomainID();
    if (!isDataDelegated())
    { // the allocation lands on the device bound to the calling thread, because
      // allocateDeviceOnly() resolves the queue through currentSubdomainID() as well
        std::lock_guard<std::mutex> lock(execution::replicaCreationMutex());
        if (!isDataDelegated())
        {
            device_only_variable_[device_id] =
                subdomain_replica_keeper_
                    .template createPtr<DeviceOnlyDiscreteVariable<DataType>>(this);
        }
    }
    return device_only_variable_[device_id]->DeviceOnlyDataField();
}
//=================================================================================================//
template <typename DataType>
void DiscreteVariable<DataType>::reallocateDataOnDevice(UnsignedInt tentative_size)
{
    reallocateData(tentative_size);
    // Every replica is grown, so that the host staging buffer stays a valid
    // destination for any of them.
    for (int device_id = 0; device_id < numberOfSubdomains(); ++device_id)
    {
        if (device_only_variable_[device_id] != nullptr)
        {
            execution::SubdomainScope scope(device_id);
            device_only_variable_[device_id]->reallocateData(this);
        }
    }
}
//=================================================================================================//
} // namespace SPH

#endif // SPHINXSYS_VARIABLE_SYCL_HPP