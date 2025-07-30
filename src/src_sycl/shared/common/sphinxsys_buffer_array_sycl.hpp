#ifndef SPHINXSYS_BUFFER_ARRAY_SYCL_HPP
#define SPHINXSYS_BUFFER_ARRAY_SYCL_HPP

#include "sphinxsys_buffer_array.h"

#include "sphinxsys_variable_array_sycl.hpp"

namespace SPH
{
//=================================================================================================//
template <typename DataType>
template <class PolicyType>
DataArray<DataType> *VariableBufferArray<DataType>::
    DelegatedBufferArray(const DeviceExecution<PolicyType> &ex_policy)
{
    DiscreteVariableArray<DataType>::DelegatedDataArray(ex_policy);
    if (!isBufferArrayDelegated())
    {
        device_buffer_array_keeper_
            .template createPtr<DeviceVariableBufferArray<DataType>>(ex_policy, this);
    }
    return delegated_buffer_array_;
}
//=================================================================================================//
template <typename DataType>
void VariableBufferArray<DataType>::setHostStagingBufferArray()
{
    if (!isBufferArrayHostStaged())
    {
        host_staging_buffer_array_ =
            device_buffer_array_keeper_.getPtr()->allocateHostStagingBufferArray();
    }
}
//=================================================================================================//
template <typename DataType>
template <class PolicyType>
DataArray<DataType> *VariableBufferArray<DataType>::
    getHostStagingBufferArray(const DeviceExecution<PolicyType> &ex_policy, UnsignedInt data_size)
{
    for (size_t i = 0; i != this->array_size_; ++i)
    {
        copyFromDevice(host_staging_buffer_array_[i], delegated_buffer_array_[i], data_size);
    }
    return host_staging_buffer_array_;
}
//=================================================================================================//
template <typename DataType>
template <class PolicyType>
DeviceVariableBufferArray<DataType>::
    DeviceVariableBufferArray(const DeviceExecution<PolicyType> &ex_policy,
                              VariableBufferArray<DataType> *host_buffer_array)
    : Entity(host_buffer_array->Name()),
      array_size_(host_buffer_array->getArraySize()),
      buffer_size_(host_buffer_array->getBufferSize()),
      device_only_buffer_array_(allocateDeviceOnly<DataArray<DataType>>(array_size_)),
      host_staging_buffer_array_(nullptr)
{
    for (size_t i = 0; i != array_size_; ++i)
    {
        device_only_buffer_array_[i] = allocateDeviceOnly<DataType>(buffer_size_);
    }
    host_buffer_array->setDelegateBufferArray(device_only_buffer_array_);
}
//=================================================================================================//
template <typename DataType>
DataArray<DataType> *DeviceVariableBufferArray<DataType>::allocateHostStagingBufferArray()
{
    if (host_staging_buffer_array_ == nullptr)
    {
        host_staging_buffer_array_ = allocateHostStaging<DataArray<DataType>>(array_size_);
        for (size_t i = 0; i != array_size_; ++i)
        {
            host_staging_buffer_array_[i] = allocateHostStaging<DataType>(buffer_size_);
        }
    }
    return host_staging_buffer_array_;
}
//=================================================================================================//
template <typename DataType>
DeviceVariableBufferArray<DataType>::~DeviceVariableBufferArray()
{
    for (size_t i = 0; i != array_size_; ++i)
    {
        freeDeviceData(device_only_buffer_array_[i]);
        if (host_staging_buffer_array_[i] != nullptr)
        {
            freeDeviceData(host_staging_buffer_array_[i]);
        }
    }
    freeDeviceData(device_only_buffer_array_);
}
//=================================================================================================//
} // namespace SPH
#endif // SPHINXSYS_BUFFER_ARRAY_SYCL_HPP
