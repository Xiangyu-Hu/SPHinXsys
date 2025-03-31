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
MultilevelLevelSet::MultilevelLevelSet(
  const ParallelDevicePolicy &par_device, BoundingBox tentative_bounds,
  MeshWithGridDataPackagesType* coarse_data, Shape &shape, SPHAdaptation &sph_adaptation)
  : BaseMeshField("LevelSet_" + shape.getName()), shape_(shape), total_levels_(1)
{
    Kernel *origin_kernel = sph_adaptation.getKernel();
    device_kernel_ = makeUnique<SingularVariable<KernelWendlandC2CK>>(
        "levelset_kernel", KernelWendlandC2CK(*origin_kernel));
    Real reference_data_spacing = coarse_data->DataSpacing() * 0.5;
    Real global_h_ratio = sph_adaptation.ReferenceSpacing() / reference_data_spacing;
    global_h_ratio_vec_.push_back(global_h_ratio);

    initializeLevel(par_device, 0, reference_data_spacing, global_h_ratio,
                    tentative_bounds, device_kernel_->DelegatedData(par_device),
                    coarse_data);
    
    configOperationExecutionPolicy(par_device, device_kernel_->DelegatedData(par_device));
}
//=================================================================================================//
MultilevelLevelSet::MultilevelLevelSet(
  const ParallelDevicePolicy &par_device, BoundingBox tentative_bounds, Real reference_data_spacing,
  size_t total_levels, Shape &shape, SPHAdaptation &sph_adaptation)
  : BaseMeshField("LevelSet_" + shape.getName()), shape_(shape), total_levels_(total_levels)
{
    Kernel *origin_kernel = sph_adaptation.getKernel();
    device_kernel_ = makeUnique<SingularVariable<KernelWendlandC2CK>>(
        "levelset_kernel", KernelWendlandC2CK(*origin_kernel));
    Real global_h_ratio = sph_adaptation.ReferenceSpacing() / reference_data_spacing;
    global_h_ratio_vec_.push_back(global_h_ratio);

    initializeLevel(par_device, 0, reference_data_spacing, global_h_ratio,
                    tentative_bounds, device_kernel_->DelegatedData(par_device));

    for (size_t level = 1; level < total_levels_; ++level) {
        reference_data_spacing *= 0.5;  // Halve the data spacing
        global_h_ratio *= 2;            // Double the ratio
        global_h_ratio_vec_.push_back(global_h_ratio);

        initializeLevel(par_device, level, reference_data_spacing, global_h_ratio,
                        tentative_bounds, device_kernel_->DelegatedData(par_device),
                        mesh_data_set_[level - 1]);
    }

    configOperationExecutionPolicy(par_device, device_kernel_->DelegatedData(par_device));
}
//=================================================================================================//
} // namespace SPH
#endif