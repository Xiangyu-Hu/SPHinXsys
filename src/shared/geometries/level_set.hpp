#ifndef LEVEL_SET_HPP
#define LEVEL_SET_HPP

#include "level_set.h"

namespace SPH
{
//=================================================================================================//
template <class ExecutionPolicy>
MultilevelLevelSet::MultilevelLevelSet(
    const ExecutionPolicy &ex_policy, BoundingBox tentative_bounds, MeshWithGridDataPackagesType* coarse_data,
    Shape &shape, SPHAdaptation &sph_adaptation)
    : BaseMeshField("LevelSet_" + shape.getName()), shape_(shape), total_levels_(1)
{
    Kernel *origin_kernel = sph_adaptation.getKernel();
    kernel_ = makeUnique<SingularVariable<KernelWendlandC2CK>>("levelset_kernel", KernelWendlandC2CK(*origin_kernel));
    Real reference_data_spacing = coarse_data->DataSpacing() * 0.5;
    Real global_h_ratio = sph_adaptation.ReferenceSpacing() / reference_data_spacing;
    global_h_ratio_vec_.push_back(global_h_ratio);

    initializeLevel(ex_policy, 0, reference_data_spacing, global_h_ratio,
                    tentative_bounds, kernel_->DelegatedData(ex_policy),
                    coarse_data);
    
    configOperationExecutionPolicy(ex_policy, kernel_->DelegatedData(ex_policy));
}
//=================================================================================================//
template <class ExecutionPolicy>
MultilevelLevelSet::MultilevelLevelSet(
    const ExecutionPolicy &ex_policy, BoundingBox tentative_bounds, Real reference_data_spacing,
    size_t total_levels, Shape &shape, SPHAdaptation &sph_adaptation)
    : BaseMeshField("LevelSet_" + shape.getName()), shape_(shape), total_levels_(total_levels)
{
    Kernel *origin_kernel = sph_adaptation.getKernel();
    kernel_ = makeUnique<SingularVariable<KernelWendlandC2CK>>("levelset_kernel", KernelWendlandC2CK(*origin_kernel));
    Real global_h_ratio = sph_adaptation.ReferenceSpacing() / reference_data_spacing;
    global_h_ratio_vec_.push_back(global_h_ratio);

    initializeLevel(ex_policy, 0, reference_data_spacing, global_h_ratio,
                    tentative_bounds, kernel_->DelegatedData(ex_policy));

    for (size_t level = 1; level < total_levels_; ++level) {
        reference_data_spacing *= 0.5;  // Halve the data spacing
        global_h_ratio *= 2;            // Double the ratio
        global_h_ratio_vec_.push_back(global_h_ratio);

        initializeLevel(ex_policy, level, reference_data_spacing, global_h_ratio,
                        tentative_bounds, kernel_->DelegatedData(ex_policy), mesh_data_set_[level - 1]);
    }

    configOperationExecutionPolicy(ex_policy, kernel_->DelegatedData(ex_policy));
}
//=================================================================================================//
template <class ExecutionPolicy, class KernelType>
void MultilevelLevelSet::configOperationExecutionPolicy(const ExecutionPolicy &ex_policy, KernelType *kernel)
{
    host_clean_interface_ = makeUnique<CleanInterface<ParallelPolicy, KernelType>>(*mesh_data_set_.back(), kernel, global_h_ratio_vec_.back());
    host_correct_topology_ = makeUnique<CorrectTopology<ParallelPolicy, KernelType>>(*mesh_data_set_.back(), kernel, global_h_ratio_vec_.back());

    device_clean_interface_ = nullptr;
    device_correct_topology_ = nullptr;
    clean_interface_ = std::bind(&CleanInterface<ParallelPolicy, KernelType>::exec, host_clean_interface_.get(), _1);
    correct_topology_ = std::bind(&CorrectTopology<ParallelPolicy, KernelType>::exec, host_correct_topology_.get(), _1);
}
//=================================================================================================//
template <class ExecutionPolicy, class KernelType>
void MultilevelLevelSet::initializeLevel(const ExecutionPolicy &ex_policy, size_t level,
                                         Real reference_data_spacing, Real global_h_ratio,
                                         BoundingBox tentative_bounds, KernelType *kernel, 
                                         MeshWithGridDataPackagesType* coarse_data)
{
    mesh_data_set_.push_back(
            mesh_data_ptr_vector_keeper_
                .template createPtr<MeshWithGridDataPackagesType>(tentative_bounds, reference_data_spacing, 4));

    RegisterMeshVariable().exec(mesh_data_set_[level]);

    if (coarse_data == nullptr) {
        MeshAllDynamicsCK<execution::ParallelPolicy, InitializeDataInACell> initialize_data_in_a_cell(*mesh_data_set_[level], shape_);
        initialize_data_in_a_cell.exec();
    } else {
        MeshAllDynamicsCK<execution::ParallelPolicy, InitializeDataInACellFromCoarse> initialize_data_in_a_cell_from_coarse(*mesh_data_set_[level], *coarse_data, shape_);
        initialize_data_in_a_cell_from_coarse.exec();
    }

    /* All initializations in `FinishDataPackages` are achieved on CPU. */
    FinishDataPackages finish_data_packages(*mesh_data_set_[level], shape_);
    finish_data_packages.exec();
    MeshInnerDynamicsCK<ExecutionPolicy, UpdateLevelSetGradient> update_level_set_gradient{*mesh_data_set_[level]};
    MeshInnerDynamicsCK<ExecutionPolicy, UpdateKernelIntegrals<KernelType>> update_kernel_integrals{*mesh_data_set_[level], kernel, global_h_ratio};
    update_level_set_gradient.exec();
    update_kernel_integrals.exec();

    registerProbes(ex_policy, level);
    cell_package_index_set_.push_back(mesh_data_set_[level]->cell_package_index_.DelegatedDataField(ex_policy));
    meta_data_cell_set_.push_back(mesh_data_set_[level]->meta_data_cell_.DelegatedDataField(ex_policy));
}
//=================================================================================================//
template <class ExecutionPolicy>
void MultilevelLevelSet::registerProbes(const ExecutionPolicy &ex_policy, size_t level)
{
    probe_signed_distance_set_.push_back(
        probe_signed_distance_vector_keeper_
            .template createPtr<ProbeSignedDistance>(ex_policy, mesh_data_set_[level]));
    probe_normal_direction_set_.push_back(
        probe_normal_direction_vector_keeper_
            .template createPtr<ProbeNormalDirection>(ex_policy, mesh_data_set_[level]));
    probe_level_set_gradient_set_.push_back(
        probe_level_set_gradient_vector_keeper_
            .template createPtr<ProbeLevelSetGradient>(ex_policy, mesh_data_set_[level]));
    probe_kernel_integral_set_.push_back(
        probe_kernel_integral_vector_keeper_
            .template createPtr<ProbeKernelIntegral>(ex_policy, mesh_data_set_[level]));
    probe_kernel_gradient_integral_set_.push_back(
        probe_kernel_gradient_integral_vector_keeper_
            .template createPtr<ProbeKernelGradientIntegral>(ex_policy, mesh_data_set_[level]));
}
//=================================================================================================//
}
  
#endif