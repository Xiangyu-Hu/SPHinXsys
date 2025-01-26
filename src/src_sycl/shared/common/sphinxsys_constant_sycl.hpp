#ifndef SPHINXSYS_CONSTANT_SYCL_HPP
#define SPHINXSYS_CONSTANT_SYCL_HPP

#include "implementation_sycl.h"
#include "sphinxsys_constant.h"

namespace SPH
{
//=================================================================================================//
template <typename GeneratorType, typename DataType>
DataType *ConstantArray<GeneratorType, DataType>::
    DelegatedData(const ParallelDevicePolicy &par_device)
{
    if (!isDataDelegated())
    {
        device_only_constant_array_keeper_
            .createPtr<DeviceOnlyConstantArray<GeneratorType, DataType>>(this);
    }
    return delegated_;
};
//=================================================================================================//
template <typename GeneratorType, typename DataType>
DeviceOnlyConstantArray<GeneratorType, DataType>::
    DeviceOnlyConstantArray(ConstantArray<GeneratorType, DataType> *host_constant)
    : Entity(host_constant->Name()), device_only_data_(nullptr)
{
    StdVec<GeneratorType *> generators = host_constant->getGenerators();
    size_t data_size = host_constant->getDataSize();
    DataType *host_data = host_constant->Data();
    for (size_t i = 0; i != data_size; ++i)
    {
        host_data[i] = DataType(ParallelDevicePolicy{}, *generators[i]);
    }
    device_only_data_ = allocateDeviceOnly<DataType>(data_size);
    copyToDevice(host_data, device_only_data_, data_size);
    host_constant->setDelegateData(device_only_data_);
}
//=================================================================================================//
template <typename GeneratorType, typename DataType>
DeviceOnlyConstantArray<GeneratorType, DataType>::~DeviceOnlyConstantArray()
{
    freeDeviceData(device_only_data_);
}
//=================================================================================================//
} // namespace SPH
#endif // SPHINXSYS_CONSTANT_SYCL_HPP
