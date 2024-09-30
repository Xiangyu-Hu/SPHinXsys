#ifndef SPHINXSYS_CONSTANT_SYCL_HPP
#define SPHINXSYS_CONSTANT_SYCL_HPP

#include "execution_sycl.h"
#include "sphinxsys_constant.h"

namespace SPH
{
//=================================================================================================//
template <typename DataType>
DeviceSharedDiscreteConstant<DataType>::
    DeviceSharedDiscreteConstant(DiscreteConstant<DataType> *host_variable)
    : Entity(host_variable->Name()), device_shared_data_field_(nullptr)
{
    size_t data_size = host_variable->getDataFieldSize();
    device_shared_data_field_ = allocateDeviceShared<DataType>(data_size);
    copyToDevice(host_variable->DataField(), device_shared_data_field_, data_size);
    host_variable->setDeviceDataField(device_shared_data_field_);
}
//=================================================================================================//
template <typename DataType>
DeviceSharedDiscreteConstant<DataType>::~DeviceSharedDiscreteConstant()
{
    freeDeviceData(device_shared_data_field_);
}
//=================================================================================================//
} // namespace SPH

#endif // SPHINXSYS_CONSTANT_SYCL_HPP