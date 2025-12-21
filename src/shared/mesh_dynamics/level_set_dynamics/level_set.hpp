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
void LevelSet::registerProbes(const ExecutionPolicy &ex_policy)
{
    probe_signed_distance_ =
        probe_signed_distance_keeper_.template createPtr<
            ProbeLevelSet<Real>>(ex_policy, *this, "LevelSet");

    probe_level_set_gradient_ =
        probe_level_set_gradient_keeper_.template createPtr<
            ProbeLevelSet<Vecd>>(ex_policy, *this, "LevelSetGradient");

    addPackageVariableToProbe<Real>("LevelSet");
    addPackageVariableToProbe<Vecd>("LevelSetGradient"); // shared with normal direction
}
//=================================================================================================//
template <class ExecutionPolicy>
void LevelSet::registerKernelIntegralProbes(const ExecutionPolicy &ex_policy)
{
    probe_kernel_integral_ =
        probe_kernel_integral_keeper_.template createPtr<
            ProbeLevelSet<Real>>(ex_policy, *this, "KernelWeight");

    probe_kernel_gradient_integral_ =
        probe_kernel_gradient_integral_keeper_.template createPtr<
            ProbeLevelSet<Vecd>>(ex_policy, *this, "KernelGradient");

    probe_kernel_second_gradient_integral_ =
        probe_kernel_second_gradient_integral_keeper_.template createPtr<
            ProbeLevelSet<Matd>>(ex_policy, *this, "KernelSecondGradient");

    addPackageVariableToProbe<Real>("KernelWeight");
    addPackageVariableToProbe<Vecd>("KernelGradient");
    addPackageVariableToProbe<Matd>("KernelSecondGradient");
}
//=================================================================================================//
template <typename DataType>
template <class ExecutionPolicy>
LevelSet::ProbeLevelSet<DataType>::ProbeLevelSet(
    const ExecutionPolicy &ex_policy, LevelSet &encloser, const std::string &variable_name)
    : BaseProbe(ex_policy, encloser, variable_name),
      global_h_ratio_(encloser.ca_global_h_ratio_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <typename DataType>
DataType LevelSet::ProbeLevelSet<DataType>::operator()(const Vecd &position)
{
    if (this->resolution_levels_ == 1)
    {
        return BaseProbe::probeInResolutionLevel(0, position);
    }

    UnsignedInt proble_level = BaseProbe::locateResolutionLevelByPackageType(1, position);
    return BaseProbe::probeInResolutionLevel(proble_level, position);
}
//=================================================================================================//
template <typename DataType>
DataType LevelSet::ProbeLevelSet<DataType>::operator()(const Vecd &position, Real h_ratio)
{
    if (this->resolution_levels_ == 1)
    {
        return BaseProbe::probeInResolutionLevel(0, position);
    }

    UnsignedInt coarse_level = 0;
    for (size_t level = this->resolution_levels_ - 1; level != 0; --level)
    {
        if (h_ratio > global_h_ratio_[level])
        {
            coarse_level = level; // jump out the loop!
            break;
        }
    }

    Real alpha = (global_h_ratio_[coarse_level + 1] - h_ratio) /
                 (global_h_ratio_[coarse_level + 1] - global_h_ratio_[coarse_level]);
    return BaseProbe::probeBetweenResolutionLevels(coarse_level, alpha, position);
}
//=================================================================================================//
} // namespace SPH
#endif
