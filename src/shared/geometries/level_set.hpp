#ifndef LEVEL_SET_HPP
#define LEVEL_SET_HPP

#include "level_set.h"

namespace SPH
{
//=================================================================================================//
template <class ExecutionPolicy, class KernelType>
void MultilevelLevelSet::configOperationExecutionPolicy(const ExecutionPolicy &ex_policy,
                                                        KernelType *kernel)
{
    sync_mesh_variable_data_ = [&]()
    { this->syncMeshVariableData(ex_policy); };
    clean_interface_keeper_ = makeUnique<CleanInterface<ExecutionPolicy, KernelType>>(
        *mesh_data_set_.back(), kernel, global_h_ratio_vec_.back());
    correct_topology_keeper_ = makeUnique<CorrectTopology<ExecutionPolicy, KernelType>>(
        *mesh_data_set_.back(), kernel, global_h_ratio_vec_.back());
}
//=================================================================================================//
template <class ExecutionPolicy, class KernelType>
void MultilevelLevelSet::initializeMeshVariables(const ExecutionPolicy &ex_policy, KernelType *kernel)
{
    for (size_t level = 0; level < total_levels_; level++)
    {
        MeshInnerDynamics<ExecutionPolicy, UpdateLevelSetGradient>
            update_level_set_gradient{*mesh_data_set_[level]};
        MeshInnerDynamics<ExecutionPolicy, UpdateKernelIntegrals<KernelType>>
            update_kernel_integrals{*mesh_data_set_[level], kernel, global_h_ratio_vec_[level]};
        update_level_set_gradient.exec();
        update_kernel_integrals.exec();

        registerProbes(ex_policy, level);
        cell_package_index_set_.push_back(
            mesh_data_set_[level]->cell_package_index_.DelegatedData(ex_policy));
        meta_data_cell_set_.push_back(
            mesh_data_set_[level]->meta_data_cell_.DelegatedData(ex_policy));
    }
}
//=================================================================================================//
template <class ExecutionPolicy>
ProbeSignedDistance MultilevelLevelSet::getProbeSignedDistance(const ExecutionPolicy &ex_policy)
{
    return ProbeSignedDistance(ex_policy, mesh_data_set_[total_levels_ - 1]);
}
//=================================================================================================//
template <class ExecutionPolicy>
ProbeNormalDirection MultilevelLevelSet::getProbeNormalDirection(const ExecutionPolicy &ex_policy)
{
    return ProbeNormalDirection(ex_policy, mesh_data_set_[total_levels_ - 1]);
}
//=================================================================================================//
template <class ExecutionPolicy>
ProbeKernelGradientIntegral MultilevelLevelSet::getProbeKernelGradientIntegral(const ExecutionPolicy &ex_policy)
{
    return ProbeKernelGradientIntegral(ex_policy, mesh_data_set_[total_levels_ - 1]);
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
    probe_kernel_second_gradient_integral_set_.push_back(
        probe_kernel_second_gradient_integral_vector_keeper_
            .template createPtr<ProbeKernelSecondGradientIntegral>(ex_policy, mesh_data_set_[level]));
}
//=================================================================================================//
} // namespace SPH

#endif