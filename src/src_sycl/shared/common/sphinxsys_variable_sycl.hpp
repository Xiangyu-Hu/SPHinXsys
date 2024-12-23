#ifndef SPHINXSYS_VARIABLE_SYCL_HPP
#define SPHINXSYS_VARIABLE_SYCL_HPP

#include "implementation_sycl.h"
#include "sphinxsys_variable.h"

namespace SPH
{
//=================================================================================================//
template <typename DataType>
DeviceSharedSingularVariable<DataType>::
    DeviceSharedSingularVariable(SingularVariable<DataType> *host_variable)
    : Entity(host_variable->Name()),
      device_shared_data_(allocateDeviceShared<DataType>(1))
{
    *device_shared_data_ = *host_variable->Data();
    host_variable->setDelegateData(device_shared_data_);
}
//=================================================================================================//
template <typename DataType>
DeviceSharedSingularVariable<DataType>::~DeviceSharedSingularVariable()
{
    freeDeviceData(device_shared_data_);
}
//=================================================================================================//
template <typename DataType>
void DiscreteVariable<DataType>::synchronizeWithDevice()
{
    if (isDataDelegated())
    {
        copyFromDevice(data_field_, device_data_field_, data_size_);
    }
}
//=================================================================================================//
template <typename DataType>
void DiscreteVariable<DataType>::synchronizeToDevice()
{
    if (isDataDelegated())
    {
        copyToDevice(data_field_, device_data_field_, data_size_);
    }
}
//=================================================================================================//
template <typename DataType>
DeviceOnlyDiscreteVariable<DataType>::
    DeviceOnlyDiscreteVariable(DiscreteVariable<DataType> *host_variable)
    : Entity(host_variable->Name()), device_only_data_field_(nullptr)
{
    size_t data_size = host_variable->getDataSize();
    device_only_data_field_ = allocateDeviceOnly<DataType>(data_size);
    copyToDevice(host_variable->Data(), device_only_data_field_, data_size);
    host_variable->setDeviceData(device_only_data_field_);
}
//=================================================================================================//
template <typename DataType>
DeviceOnlyDiscreteVariable<DataType>::~DeviceOnlyDiscreteVariable()
{
    freeDeviceData(device_only_data_field_);
}
//=================================================================================================//
template <typename DataType>
void DeviceOnlyDiscreteVariable<DataType>::
    reallocateData(DiscreteVariable<DataType> *host_variable)
{
    freeDeviceData(device_only_data_field_);
    size_t new_host_variable_size = host_variable->getDataSize();
    device_only_data_field_ = allocateDeviceOnly<DataType>(new_host_variable_size);
    host_variable->setDeviceData(device_only_data_field_);
}
//=================================================================================================//
template <typename DataType>
DataType *DiscreteVariable<DataType>::DelegatedData(const ParallelDevicePolicy &par_device)
{
    if (!isDataDelegated())
    {
        device_only_variable_ =
            device_only_variable_keeper_
                .createPtr<DeviceOnlyDiscreteVariable<DataType>>(this);
    }
    return device_data_field_;
}
//=================================================================================================//
template <typename DataType>
void DiscreteVariable<DataType>::reallocateData(
    const ParallelDevicePolicy &par_device, size_t tentative_size)
{
    if (data_size_ < tentative_size)
    {
        reallocateData(tentative_size);
        device_only_variable_->reallocateData(this);
    }
}
//=================================================================================================//
} // namespace SPH

#endif // SPHINXSYS_VARIABLE_SYCL_HPP