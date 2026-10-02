#ifndef SPHINXSYS_CONSTANT_SYCL_HPP
#define SPHINXSYS_CONSTANT_SYCL_HPP

#include "implementation_sycl.h"
#include "sphinxsys_constant.h"

namespace SPH
{
//=================================================================================================//
template <typename DataType>
DataType *ConstantArray<DataType>::DelegatedOnDevice(const SYCLDevicePolicy &ex_policy)
{
    if (!isDataDelegated())
    {
        device_only_constant_array_keeper_
            .createPtr<DeviceOnlyConstantArray<DataType>>(SYCLDevicePolicy{}, this);
    }
    return delegated_;
};
//=================================================================================================//
template <typename DataType>
DeviceOnlyConstantArray<DataType>::DeviceOnlyConstantArray(
    const SYCLDevicePolicy &ex_policy, ConstantArray<DataType> *host_constant)
    : Quantity(host_constant->Name()), device_only_data_(nullptr)
{
    size_t data_size = host_constant->getSize();
    DataType *host_data = host_constant->Data();
    device_only_data_ = allocateDeviceOnly<DataType>(data_size);
    copyToDevice(host_data, device_only_data_, data_size);
    host_constant->setDelegateData(device_only_data_);
}
//=================================================================================================//
template <typename DataType>
DeviceOnlyConstantArray<DataType>::~DeviceOnlyConstantArray()
{
    freeDeviceData(device_only_data_);
}
//=================================================================================================//
template <typename GeneratorType, typename ComputingKernelType>
ComputingKernelType *ComputingKernelArray<GeneratorType, ComputingKernelType>::DelegatedOnDevice(
    const SYCLDevicePolicy &ex_policy)
{
    if (!isDataDelegated())
    {
        device_only_kernel_array_keeper_.createPtr<DeviceOnlyComputingKernelArray<
            GeneratorType, ComputingKernelType>>(SYCLDevicePolicy{}, this);
    }
    return delegated_;
}
//=================================================================================================//
template <typename GeneratorType, typename ComputingKernelType>
DeviceOnlyComputingKernelArray<GeneratorType, ComputingKernelType>::DeviceOnlyComputingKernelArray(
    const SYCLDevicePolicy &ex_policy,
    ComputingKernelArray<GeneratorType, ComputingKernelType> *host_constant)
    : Quantity(host_constant->Name()), device_only_data_(nullptr)
{
    size_t data_size = host_constant->getSize();
    StdVec<GeneratorType *> generators = host_constant->getGenerators();
    ComputingKernelType *host_data = host_constant->Data();
    for (size_t i = 0; i != data_size; ++i)
    {
        host_data[i] = ComputingKernelType(ex_policy, *generators[i]);
    }
    device_only_data_ = allocateDeviceOnly<ComputingKernelType>(data_size);
    copyToDevice(host_data, device_only_data_, data_size);
    host_constant->setDelegateData(device_only_data_);
}
//=================================================================================================//
template <typename GeneratorType, typename ComputingKernelType>
DeviceOnlyComputingKernelArray<GeneratorType, ComputingKernelType>::~DeviceOnlyComputingKernelArray()
{
    freeDeviceData(device_only_data_);
}
//=================================================================================================//
} // namespace SPH
#endif // SPHINXSYS_CONSTANT_SYCL_HPP
