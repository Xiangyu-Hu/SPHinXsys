#ifndef LEVEL_SET_SYCL_HPP
#define LEVEL_SET_SYCL_HPP

#include "level_set.h"

namespace SPH
{
//=================================================================================================//
void MultilevelLevelSet::configOperationExecutionPolicy(const ParallelDevicePolicy &par_device,
                                                        Kernel *kernel)
{
    device_clean_interface_ = makeUnique<CleanInterface<ParallelDevicePolicy>>(*mesh_data_set_.back(), kernel_wrapper.getKernel(par_device), global_h_ratio_vec_.back());
    device_correct_topology_ = makeUnique<CorrectTopology<ParallelDevicePolicy>>(*mesh_data_set_.back(), kernel_wrapper.getKernel(par_device), global_h_ratio_vec_.back());

    host_clean_interface_ = nullptr;
    host_correct_topology_ = nullptr;
    clean_interface_ = std::bind(&CleanInterface<ParallelPolicy>::exec, device_clean_interface_.get(), _1);
    correct_topology_ = std::bind(&CorrectTopology<ParallelPolicy>::exec, device_correct_topology_.get(), _1);
}
//=================================================================================================//
} // namespace SPH
#endif