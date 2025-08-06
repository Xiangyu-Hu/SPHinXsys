#ifndef LEVEL_SET_SYCL_HPP
#define LEVEL_SET_SYCL_HPP

#include "implementation_sycl.h"
#include "level_set.h"
#include "mesh_iterators_sycl.hpp"
#include "sphinxsys_constant_sycl.hpp"
#include "sphinxsys_variable_sycl.hpp"
namespace SPH
{
//=================================================================================================//
void MultilevelLevelSet::finishInitialization(const ParallelDevicePolicy &par_device)
{
    initializeMeshVariables(par_device, kernel_->DelegatedData(par_device));
    if (usage_type_ == UsageType::Volumetric)
    {
        initializeKernelIntegralVariables(par_device, kernel_->DelegatedData(par_device));
        configLevelSetPostProcesses(par_device, kernel_->DelegatedData(par_device));
    }
    sync_mesh_variable_data_ = [&]()
    { this->syncMeshVariableData(par_device); };
}
//=================================================================================================//
} // namespace SPH
#endif