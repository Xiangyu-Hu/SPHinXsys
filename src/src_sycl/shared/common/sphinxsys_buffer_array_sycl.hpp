#ifndef SPHINXSYS_BUFFER_ARRAY_SYCL_HPP
#define SPHINXSYS_BUFFER_ARRAY_SYCL_HPP

#include "sphinxsys_buffer_array.h"

#include "sphinxsys_variable_array_sycl.hpp"

namespace SPH
{
//=================================================================================================//
template <typename DataType>
template <class PolicyType>
DataArray<DataType> *VariableBufferArray<DataType>::DelegatedBufferArrayOnDevice()
{
    DiscreteVariableArray<DataType>::DelegatedDataArray(DeviceExecution<PolicyType>{}); // check variable array first
    if (!isBufferArrayDelegated())
    {
        device_buffer_array_keeper_
            .template createPtr<DeviceVariableBufferArray<DataType>>(DeviceExecution<PolicyType>{}, this);
    }
    return delegated_buffer_array_;
}
//=================================================================================================//
template <typename DataType>
DataArray<DataType> *VariableBufferArray<DataType>::synchronizeBufferArrayFromDevice(UnsignedInt data_size)
{
    if (!isBufferArrayHostStaged())
    {
        device_buffer_array_keeper_.getPtr()->allocateHostStagingBufferArray(this);
    }
    device_buffer_array_keeper_.getPtr()->synchronizeHostStagingBufferArray(data_size);
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
      device_only_ptr_(nullptr), host_staging_ptr_(nullptr),
      host_staging_buffer_array_(nullptr)
{
    device_only_ptr_ = allocateDeviceOnly<DataType>(buffer_size_ * array_size_);
    for (size_t i = 0; i != array_size_; ++i)
    {
        copyToDevice(device_only_ptr_ + i * buffer_size_, device_only_buffer_array_ + i, 1);
    }
    host_buffer_array->setDelegateBufferArray(device_only_buffer_array_);
}
//=================================================================================================//
template <typename DataType>
void DeviceVariableBufferArray<DataType>::
    allocateHostStagingBufferArray(VariableBufferArray<DataType> *host_buffer_array)
{
    if (host_staging_buffer_array_ == nullptr)
    {
        host_staging_buffer_array_ = allocateHostStaging<DataArray<DataType>>(array_size_);
        host_staging_ptr_ = allocateHostStaging<DataType>(buffer_size_ * array_size_);
        for (size_t i = 0; i != array_size_; ++i)
        {
            host_staging_buffer_array_[i] = host_staging_ptr_ + i * buffer_size_;
        }
        host_buffer_array->setHostStagingBufferArray(host_staging_buffer_array_);
    }
}
//=================================================================================================//
template <typename DataType>
void DeviceVariableBufferArray<DataType>::synchronizeHostStagingBufferArray(UnsignedInt data_size)
{
    for (size_t i = 0; i != array_size_; ++i)
    {
        UnsignedInt offset = i * buffer_size_;
        copyFromDevice(host_staging_ptr_ + offset, device_only_ptr_ + offset, data_size);
    }
}
//=================================================================================================//
template <typename DataType>
DeviceVariableBufferArray<DataType>::~DeviceVariableBufferArray()
{
    freeDeviceData(device_only_buffer_array_);
    freeDeviceData(device_only_ptr_);

    if (host_staging_buffer_array_ != nullptr)
    {
        freeDeviceData(host_staging_buffer_array_);
        freeDeviceData(host_staging_ptr_);
    }
}
//=================================================================================================//
} // namespace SPH
#endif // SPHINXSYS_BUFFER_ARRAY_SYCL_HPP
