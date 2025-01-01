#ifndef SPHINXSYS_CONSTANT_SYCL_HPP
#define SPHINXSYS_CONSTANT_SYCL_HPP

#include "implementation_sycl.h"
#include "sphinxsys_constant.h"

namespace SPH
{
//=================================================================================================//
template <typename GeneratorType, typename FunctionType>
DataType *Constant<GeneratorType, FunctionType>::
    DelegatedData(const ParallelDevicePolicy &par_device)
{
    if (!isDataDelegated())
    {
        device_shared_singular_constant_keeper_
            .createPtr<DeviceOnlyConstant<GeneratorType, FunctionType>>(this);
    }
    return delegated_;
};
//=================================================================================================//
template <typename GeneratorType, typename FunctionType>
DeviceOnlyConstant<GeneratorType, FunctionType>::
    DeviceOnlyConstant(DiscreteConstant<GeneratorType, FunctionType> *host_constant)
    : Entity(host_constant->Name()), device_only_data_(nullptr)
{
    GeneratorType *generator = host_constant->getGenerator();
    size_t data_size = host_constant->getDataSize();
    device_only_data_ = allocateDeviceShared<DataType>(data_size);
    for (size_t i = 0; i != data_size; ++i)
    {
        DataType data = generator[i]->getFunction<DataType>(ParallelDevicePolicy{});
        copyToDevice(data, device_only_data_ + i, 1);
    }

    host_constant->setDeviceData(device_only_data_);
}
//=================================================================================================//
template <typename GeneratorType, typename FunctionType>
DeviceOnlyConstant<GeneratorType, FunctionType>::~DeviceOnlyConstant()
{
    freeDeviceData(device_only_data_);
}
//=================================================================================================//
} // namespace SPH

#endif // SPHINXSYS_CONSTANT_SYCL_HPP