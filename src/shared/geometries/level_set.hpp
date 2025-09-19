#ifndef LEVEL_SET_HPP
#define LEVEL_SET_HPP

#include "level_set.h"

namespace SPH
{
//=================================================================================================//
template <class ExecutionPolicy>
void MultilevelLevelSet::configLevelSetPostProcesses(const ExecutionPolicy &ex_policy)
{
    clean_interface_keeper_ = makeUnique<CleanInterface<ExecutionPolicy>>(
        *mesh_data_set_.back(), *neighbor_method_set_.back(), refinement_ratio_);
    correct_topology_keeper_ = makeUnique<CorrectTopology<ExecutionPolicy>>(
        *mesh_data_set_.back(), *neighbor_method_set_.back());
}
//=================================================================================================//
template <class ExecutionPolicy>
void MultilevelLevelSet::initializeMeshVariables(const ExecutionPolicy &ex_policy)
{
    for (size_t level = 0; level < total_levels_; level++)
    {
        MeshInnerDynamics<ExecutionPolicy, UpdateLevelSetGradient>
            update_level_set_gradient{*mesh_data_set_[level]};
        update_level_set_gradient.exec();
    }
}
//=================================================================================================//
template <class ExecutionPolicy>
void MultilevelLevelSet::syncMeshVariablesToWrite(ExecutionPolicy &ex_policy)
{
    for (size_t l = 0; l != total_levels_; l++)
        mesh_data_set_[l]->syncMeshVariablesToWrite(ex_policy);
}
//=================================================================================================//
template <class ExecutionPolicy>
void MultilevelLevelSet::syncBKGMeshVariablesToWrite(ExecutionPolicy &ex_policy)
{
    for (size_t l = 0; l != total_levels_; l++)
        mesh_data_set_[l]->syncBKGMeshVariablesToWrite(ex_policy);
}
//=================================================================================================//
template <class ExecutionPolicy>
void MultilevelLevelSet::syncMeshVariablesToProbe(ExecutionPolicy &ex_policy)
{
    for (size_t l = 0; l != total_levels_; l++)
        mesh_data_set_[l]->syncMeshVariablesToProbe(ex_policy);
}
//=================================================================================================//
template <class ExecutionPolicy>
void MultilevelLevelSet::finishInitialization(const ExecutionPolicy &ex_policy, UsageType usage_type)
{
    initializeMeshVariables(ex_policy);
    registerProbes(execution::par_host); // register probes on host

    if (usage_type == UsageType::Volumetric)
    {
        initializeKernelIntegralVariables(ex_policy);
        registerKernelIntegralProbes(execution::par_host); // register probes on host
        configLevelSetPostProcesses(ex_policy);
    }

    sync_mesh_variables_to_write_ = [&]() // for latter usage
    { this->syncMeshVariablesToWrite(ex_policy); };

    sync_bkg_mesh_variables_to_write_ = [&]() // for latter usage
    { this->syncBKGMeshVariablesToWrite(ex_policy); };

    sync_mesh_variables_to_probe_ = [&]() // for latter usage
    { this->syncMeshVariablesToProbe(ex_policy); };

    this->syncMeshVariablesToProbe(ex_policy);
}
//=================================================================================================//
template <class ExecutionPolicy>
void MultilevelLevelSet::initializeKernelIntegralVariables(const ExecutionPolicy &ex_policy)
{
    for (size_t level = 0; level < total_levels_; level++)
    {
        MeshInnerDynamics<ExecutionPolicy, UpdateKernelIntegrals>
            update_kernel_integrals{*mesh_data_set_[level], *neighbor_method_set_[level]};
        update_kernel_integrals.exec();
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
void MultilevelLevelSet::registerProbes(const ExecutionPolicy &ex_policy)
{
    for (size_t level = 0; level < total_levels_; level++)
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

        mesh_data_set_[level]->addMeshVariableToProbe<Real>("LevelSet");
        mesh_data_set_[level]->addMeshVariableToProbe<Vecd>("LevelSetGradient"); // shared with normal direction

        cell_pkg_index_set_.push_back(mesh_data_set_[level]->getCellPackageIndex().DelegatedData(ex_policy));
        pkg_cell_info_set_.push_back(mesh_data_set_[level]->dvPkgCellInfo().DelegatedData(ex_policy));
    }
}
//=================================================================================================//
template <class ExecutionPolicy>
void MultilevelLevelSet::registerKernelIntegralProbes(const ExecutionPolicy &ex_policy)
{
    for (size_t level = 0; level < total_levels_; level++)
    {
        probe_kernel_integral_set_.push_back(
            probe_kernel_integral_vector_keeper_
                .template createPtr<ProbeKernelIntegral>(ex_policy, mesh_data_set_[level]));
        probe_kernel_gradient_integral_set_.push_back(
            probe_kernel_gradient_integral_vector_keeper_
                .template createPtr<ProbeKernelGradientIntegral>(ex_policy, mesh_data_set_[level]));
        probe_kernel_second_gradient_integral_set_.push_back(
            probe_kernel_second_gradient_integral_vector_keeper_
                .template createPtr<ProbeKernelSecondGradientIntegral>(ex_policy, mesh_data_set_[level]));
        mesh_data_set_[level]->addMeshVariableToProbe<Real>("KernelWeight");
        mesh_data_set_[level]->addMeshVariableToProbe<Vecd>("KernelGradient");
        mesh_data_set_[level]->addMeshVariableToProbe<Matd>("KernelSecondGradient");
    }
}
//=================================================================================================//
template <typename DataType>
void MultilevelLevelSet::addMeshVariableToWrite(const std::string &variable_name)
{
    for (size_t level = 0; level < total_levels_; level++)
    {
        mesh_data_set_[level]->addMeshVariableToWrite<DataType>(variable_name);
    }
}
//=================================================================================================//
template <typename DataType>
void MultilevelLevelSet::addBKGMeshVariableToWrite(const std::string &variable_name)
{
    for (size_t level = 0; level < total_levels_; level++)
    {
        mesh_data_set_[level]->addBKGMeshVariableToWrite<DataType>(variable_name);
    }
}
//=================================================================================================//
} // namespace SPH
#endif
