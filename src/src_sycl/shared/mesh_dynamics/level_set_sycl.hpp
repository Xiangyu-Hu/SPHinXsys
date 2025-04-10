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
void MutilevelLevelSet::finishInitialization(const ParallelDevicePolicy &par_device)
{
    device_kernel_ = makeUnique<SingularVariable<KernelWendlandC2CK>>(
        "levelset_kernel", KernelWendlandC2CK(*host_kernel_));
    initializeMeshVariables(par_device, device_kernel_);
    configOperationExecutionPolicy(par_device, device_kernel_);
}
//=================================================================================================//
} // namespace SPH
#endif