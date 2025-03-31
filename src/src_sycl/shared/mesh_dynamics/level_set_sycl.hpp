#ifndef LEVEL_SET_SYCL_HPP
#define LEVEL_SET_SYCL_HPP

#include "mesh_dynamics_sycl.hpp"
#include "execution_sycl.h"
#include "sphinxsys_variable_sycl.hpp"
#include "sphinxsys_constant_sycl.hpp"
#include "level_set.hpp"
namespace SPH
{
//=================================================================================================//
template <class KernelType>
void MultilevelLevelSet::configOperationExecutionPolicy(const ParallelDevicePolicy &par_device,
                                                        KernelType *kernel)
{
    device_clean_interface_ = makeUnique<CleanInterface<ParallelDevicePolicy, KernelType>>(
            *mesh_data_set_.back(), kernel, global_h_ratio_vec_.back());
    device_correct_topology_ = makeUnique<CorrectTopology<ParallelDevicePolicy, KernelType>>(
            *mesh_data_set_.back(), kernel, global_h_ratio_vec_.back());

    host_clean_interface_ = nullptr;
    host_correct_topology_ = nullptr;
    clean_interface_ = std::bind(
        &CleanInterface<ParallelDevicePolicy, KernelType>::exec, device_clean_interface_.get(), _1);
    correct_topology_ = std::bind(
        &CorrectTopology<ParallelDevicePolicy, KernelType>::exec, device_correct_topology_.get(), _1);
}
//=================================================================================================//
} // namespace SPH
#endif