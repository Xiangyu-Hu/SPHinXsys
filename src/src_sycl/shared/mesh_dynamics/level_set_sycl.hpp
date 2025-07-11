#ifndef LEVEL_SET_SYCL_HPP
#define LEVEL_SET_SYCL_HPP

#include "mesh_iterators_sycl.hpp"
#include "implementation_sycl.h"
#include "sphinxsys_variable_sycl.hpp"
#include "sphinxsys_constant_sycl.hpp"
#include "level_set.h"
namespace SPH
{
//=================================================================================================//
void MultilevelLevelSet::finishInitialization(const ParallelDevicePolicy &par_device)
{
    device_kernel_ = makeUnique<SingularVariable<KernelTabulatedCK>>(
        "levelset_kernel", KernelTabulatedCK(*host_kernel_));
    initializeMeshVariables(par_device, device_kernel_->DelegatedData(par_device));
    configOperationExecutionPolicy(par_device, device_kernel_->DelegatedData(par_device));
}
//=================================================================================================//
} // namespace SPH
#endif