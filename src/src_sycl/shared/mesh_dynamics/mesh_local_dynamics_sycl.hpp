#ifndef MESH_LOCAL_DYNAMICS_SYCL_H
#define MESH_LOCAL_DYNAMICS_SYCL_H

#include "mesh_local_dynamics.hpp"
#include "execution_sycl.h"

template <class KernelType>
KernelType* KernelWrapper<KernelType>::getKernel(const ParallelDevicePolicy &par_device)
{
    if(!device_kernel_ptr_){
        device_kernel_ptr_ = allocateDeviceOnly<KernelType>(1);
        copyToDevice(host_kernel_ptr_, device_kernel_ptr_, 1);
    }
    return device_kernel_ptr_;
}

#endif