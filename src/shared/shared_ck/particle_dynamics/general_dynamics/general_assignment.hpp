#ifndef GENERAL_ASSIGNMENT_HPP
#define GENERAL_ASSIGNMENT_HPP

#include "general_assignment.h"

#include "base_body_part.h"
#include "complex_geometry.h"

namespace SPH
{
//=================================================================================================//
template <class DynamicsIdentifier, typename AssignmentFunctionType>
template <typename... Args>
VariableAssignment<DynamicsIdentifier, AssignmentFunctionType>::VariableAssignment(
    DynamicsIdentifier &identifier, Args &&...args)
    : BaseLocalDynamics<DynamicsIdentifier>(identifier),
      assignment_method_(this->particles_, std::forward<Args>(args)...) {}
//=================================================================================================//
template <class DynamicsIdentifier, typename AssignmentFunctionType>
template <class ExecutionPolicy, class EncloserType>
VariableAssignment<DynamicsIdentifier, AssignmentFunctionType>::UpdateKernel::UpdateKernel(
    const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : assign_(ex_policy, encloser.assignment_method_) {}
//=================================================================================================//
template <typename ConditionType>
template <typename... Args>
VariableAssignment<OrientedBoxByCell, ConditionType>::VariableAssignment(
    OrientedBoxByCell &oriented_box_part, Args &&...args)
    : BaseLocalDynamics<OrientedBoxByCell>(oriented_box_part),
      sv_oriented_box_(oriented_box_part.svOrientedBox()),
      assignment_method_(this->particles_, std::forward<Args>(args)...),
      dv_pos_(particles_->getVariableByName<Vecd>("Position")) {}
//=================================================================================================//
template <typename ConditionType>
template <class ExecutionPolicy, class EncloserType>
VariableAssignment<OrientedBoxByCell, ConditionType>::UpdateKernel::UpdateKernel(
    const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : oriented_box_(encloser.sv_oriented_box_->DelegatedData(ex_policy)),
      assign_(ex_policy, encloser.assignment_method_),
      pos_(encloser.dv_pos_->DelegatedDataView(ex_policy)) {}
//=================================================================================================//
template <typename ConditionType>
void VariableAssignment<OrientedBoxByCell, ConditionType>::UpdateKernel::update(
    size_t index_i, Real dt)
{
    if (oriented_box_->checkContain(pos_[index_i]))
    {
        assign_(index_i);
    }
}
//=================================================================================================//
template <typename DistributionType>
template <typename... Args>
SpatialDistribution<DistributionType>::SpatialDistribution(
    BaseParticles *particles, const std::string &variable_name, Args &&...args)
    : dv_variable_(particles->template getVariableByName<DataType>(variable_name)),
      dv_pos_(particles->template getVariableByName<Vecd>("Position")),
      distribution_(std::forward<Args>(args)...) {}
//=================================================================================================//
template <typename DistributionType>
template <class ExecutionPolicy, class EncloserType>
SpatialDistribution<DistributionType>::SpatialDistribution::ComputingKernel::
    ComputingKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : variable_(encloser.dv_variable_->DelegatedDataView(ex_policy)),
      pos_(encloser.dv_pos_->DelegatedDataView(ex_policy)),
      distribution_(encloser.distribution_) {}
//=================================================================================================//
template <typename DataType>
ConstantValue<DataType>::ConstantValue(
    BaseParticles *particles, const std::string &variable_name, DataType constant_value)
    : dv_variable_(particles->template getVariableByName<DataType>(variable_name)),
      constant_value_(constant_value) {}
//=================================================================================================//
template <typename DataType>
template <class ExecutionPolicy, class EncloserType>
ConstantValue<DataType>::ComputingKernel::
    ComputingKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : variable_(encloser.dv_variable_->DelegatedDataView(ex_policy)),
      constant_value_(encloser.constant_value_) {}
//=================================================================================================//
} // namespace SPH
#endif // GENERAL_ASSIGNMENT_HPP
