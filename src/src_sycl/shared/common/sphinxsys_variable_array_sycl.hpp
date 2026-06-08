#ifndef SPHINXSYS_VARIABLE_ARRAY_SYCL_HPP
#define SPHINXSYS_VARIABLE_ARRAY_SYCL_HPP

#include "implementation_sycl.h"
#include "sphinxsys_variable_array.h"

namespace SPH
{
//=================================================================================================//
template <typename DataType>
template <class PolicyType>
MultiEntryView<DataType> *VariableArray<DataType>::DelegatedOnDevice()
{
    if (!isVariableArrayViewDelegated())
    {
        device_only_variable_array_ = device_only_variable_array_keeper_.createPtr<
            DeviceOnlyVariableArray<DataType>>(DeviceExecution<PolicyType>{}, this);
    }
    return device_only_variable_array_->DeviceOnlyMultiEntryView();
}
//=================================================================================================//
template <typename DataType>
template <class PolicyType>
DeviceOnlyVariableArray<DataType>::
    DeviceOnlyVariableArray(const DeviceExecution<PolicyType> &ex_policy,
                            VariableArray<DataType> *host_variable_array)
    : Quantity(host_variable_array->Name()), device_only_multi_entry_view_(nullptr)
{
    StdVec<DiscreteVariable<DataType> *> host_variables = host_variable_array->getVariables();
    size_t data_size = host_variable_array->getArraySize();
    MultiEntryView<DataType> *multi_entry_view = host_variable_array->getArrayData();
    device_only_multi_entry_view_ = allocateDeviceOnly<MultiEntryView<DataType>>(data_size);
    for (size_t i = 0; i != data_size; ++i)
    {
        multi_entry_view[i].setData(host_variables[i]->DelegatedData(ex_policy));
        copyToDevice(multi_entry_view, device_only_multi_entry_view_ + i, 1);
    }
}
//=================================================================================================//
template <typename DataType>
DeviceOnlyVariableArray<DataType>::~DeviceOnlyVariableArray()
{
    freeDeviceData(device_only_multi_entry_view_);
}
//=================================================================================================//
} // namespace SPH
#endif // SPHINXSYS_VARIABLE_ARRAY_SYCL_HPP
