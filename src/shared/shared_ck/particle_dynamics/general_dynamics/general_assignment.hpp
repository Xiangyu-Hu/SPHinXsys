#ifndef GENERAL_ASSIGNMENT_HPP
#define GENERAL_ASSIGNMENT_HPP

#include "general_assignment.h"

namespace SPH
{
//=================================================================================================//
template <class DynamicsIdentifier, typename AssignmentFunctionType>
template <typename... Args>
VariableAssignment<DynamicsIdentifier, AssignmentFunctionType>::VariableAssignment(
    DynamicsIdentifier &identifier, const std::string &variable_name, Args &&...args)
    : BaseLocalDynamics<DynamicsIdentifier>(identifier),
      dv_variable_(this->particles_->template registerStateVariable<DataType>(variable_name)),
      assignment_method_(this->particles_, std::forward<Args>(args)...) {}
//=================================================================================================//
template <class DynamicsIdentifier, typename AssignmentFunctionType>
template <class ExecutionPolicy, class EncloserType>
VariableAssignment<DynamicsIdentifier, AssignmentFunctionType>::UpdateKernel::UpdateKernel(
    const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : variable_(encloser.dv_variable_->DelegatedData(ex_policy)),
      assign_(ex_policy, encloser.assignment_method_) {}
//=================================================================================================//
template <class DynamicsIdentifier, typename AssignmentFunctionType>
void VariableAssignment<DynamicsIdentifier, AssignmentFunctionType>::UpdateKernel::update(
    UnsignedInt index_i, Real dt)
{
    variable_[index_i] = assign_(index_i);
}
//=================================================================================================//
} // namespace SPH
#endif // GENERAL_ASSIGNMENT_HPP
