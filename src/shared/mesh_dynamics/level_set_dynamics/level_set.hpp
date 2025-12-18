#ifndef LEVEL_SET_HPP
#define LEVEL_SET_HPP

#include "level_set.h"

namespace SPH
{
//=================================================================================================//
template <class ExecutionPolicy>
void LevelSet::configLevelSetPostProcesses(const ExecutionPolicy &ex_policy)
{
    clean_interface_keeper_ = makeUnique<CleanInterface<ExecutionPolicy>>(
        *this, resolution_levels_ - 1, *neighbor_method_set_.back(), refinement_ratio_);
    correct_topology_keeper_ = makeUnique<CorrectTopology<ExecutionPolicy>>(
        *this, resolution_levels_ - 1, *neighbor_method_set_.back());
}
//=================================================================================================//
template <class ExecutionPolicy>
void LevelSet::initializePackageVariables(const ExecutionPolicy &ex_policy)
{
    for (size_t level = 0; level < resolution_levels_; level++)
    {
        MeshInnerDynamics<ExecutionPolicy, UpdateLevelSetGradient>
            update_level_set_gradient(*this, level);
        update_level_set_gradient.exec();
    }
}
//=================================================================================================//
template <class ExecutionPolicy>
void LevelSet::finishInitialization(const ExecutionPolicy &ex_policy, UsageType usage_type)
{
    initializePackageVariables(ex_policy);
    registerProbes(execution::par_host); // register probes on host

    if (usage_type == UsageType::Volumetric)
    {
        initializeKernelIntegralVariables(ex_policy);
        registerKernelIntegralProbes(execution::par_host); // register probes on host
        configLevelSetPostProcesses(ex_policy);
    }

    sync_mesh_variables_to_write_ = [&]() // for latter usage
    { this->syncPackageVariablesToWrite(ex_policy); 
    this->syncCellVariablesToWrite(ex_policy); };

    sync_mesh_variables_to_probe_ = [&]() // for latter usage
    { this->syncPackageVariablesToProbe(ex_policy); };

    this->syncPackageVariablesToProbe(ex_policy);
}
//=================================================================================================//
template <class ExecutionPolicy>
void LevelSet::initializeKernelIntegralVariables(const ExecutionPolicy &ex_policy)
{
    for (size_t level = 0; level < resolution_levels_; level++)
    {
        MeshInnerDynamics<ExecutionPolicy, UpdateKernelIntegrals>
            update_kernel_integrals(*this, level, *neighbor_method_set_[level]);
        update_kernel_integrals.exec();
    }
}
//=================================================================================================//
template <class ExecutionPolicy>
ProbeSignedDistance LevelSet::getProbeSignedDistance(const ExecutionPolicy &ex_policy)
{
    return ProbeSignedDistance(ex_policy, this, resolution_levels_ - 1);
}
//=================================================================================================//
template <class ExecutionPolicy>
ProbeNormalDirection LevelSet::getProbeNormalDirection(const ExecutionPolicy &ex_policy)
{
    return ProbeNormalDirection(ex_policy, this, resolution_levels_ - 1);
}
//=================================================================================================//
template <class ExecutionPolicy>
ProbeKernelGradientIntegral LevelSet::getProbeKernelGradientIntegral(const ExecutionPolicy &ex_policy)
{
    return ProbeKernelGradientIntegral(ex_policy, this, resolution_levels_ - 1);
}
//=================================================================================================//
template <class ExecutionPolicy>
void LevelSet::registerProbes(const ExecutionPolicy &ex_policy)
{
    for (size_t level = 0; level < resolution_levels_; level++)
    {
        probe_signed_distance_set_.push_back(
            probe_signed_distance_vector_keeper_
                .template createPtr<ProbeSignedDistance>(ex_policy, this, level));
        probe_normal_direction_set_.push_back(
            probe_normal_direction_vector_keeper_
                .template createPtr<ProbeNormalDirection>(ex_policy, this, level));
        probe_level_set_gradient_set_.push_back(
            probe_level_set_gradient_vector_keeper_
                .template createPtr<ProbeLevelSetGradient>(ex_policy, this, level));

        addMeshVariableToProbe<Real>("LevelSet");
        addMeshVariableToProbe<Vecd>("LevelSetGradient"); // shared with normal direction

        mesh_index_handler_set_.push_back(&getMeshLevel(level));
        cell_pkg_index_ = mcv_cell_pkg_index_->DelegatedData(ex_policy);
        pkg_type_ = dv_pkg_type_->DelegatedData(ex_policy);
    }
}
//=================================================================================================//
template <class ExecutionPolicy>
void LevelSet::registerKernelIntegralProbes(const ExecutionPolicy &ex_policy)
{
    for (size_t level = 0; level < resolution_levels_; level++)
    {
        probe_kernel_integral_set_.push_back(
            probe_kernel_integral_vector_keeper_
                .template createPtr<ProbeKernelIntegral>(ex_policy, this, level));
        probe_kernel_gradient_integral_set_.push_back(
            probe_kernel_gradient_integral_vector_keeper_
                .template createPtr<ProbeKernelGradientIntegral>(ex_policy, this, level));
        probe_kernel_second_gradient_integral_set_.push_back(
            probe_kernel_second_gradient_integral_vector_keeper_
                .template createPtr<ProbeKernelSecondGradientIntegral>(ex_policy, this, level));
        addMeshVariableToProbe<Real>("KernelWeight");
        addMeshVariableToProbe<Vecd>("KernelGradient");
        addMeshVariableToProbe<Matd>("KernelSecondGradient");
    }
}
//=================================================================================================//
} // namespace SPH
#endif
