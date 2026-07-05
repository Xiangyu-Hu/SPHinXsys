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
template <class DynamicsIdentifier, typename AssignmentFunctionType>
template <typename... Args>
VariableEntryAssignment<DynamicsIdentifier, AssignmentFunctionType>::VariableEntryAssignment(
    DynamicsIdentifier &identifier, const std::string &variable_name,
    const std::string &entry_name, Args &&...args)
    : BaseLocalDynamics<DynamicsIdentifier>(identifier),
      dv_variable_(this->particles_->template registerStateVariable<DataType>(variable_name)),
      entry_name_(entry_name), assignment_method_(this->particles_, std::forward<Args>(args)...){};
//=================================================================================================//
template <class DynamicsIdentifier, typename AssignmentFunctionType>
template <class ExecutionPolicy, class EncloserType>
VariableEntryAssignment<DynamicsIdentifier, AssignmentFunctionType>::UpdateKernel::UpdateKernel(
    const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : entry_(encloser.dv_variable_->DelegatedEntryView(ex_policy, encloser.entry_name_)),
      assign_(ex_policy, encloser.assignment_method_) {}
//=================================================================================================//
template <class DynamicsIdentifier, typename AssignmentFunctionType>
void VariableEntryAssignment<DynamicsIdentifier, AssignmentFunctionType>::UpdateKernel::update(
    UnsignedInt index_i, Real dt)
{
    entry_[index_i] = assign_(index_i);
}
//=================================================================================================//
template <class DynamicsIdentifier, typename AssignmentFunctionType>
template <typename... Args>
MultiEntryVariableAssignment<DynamicsIdentifier, AssignmentFunctionType>::MultiEntryVariableAssignment(
    DynamicsIdentifier &identifier, const std::string &variable_name, Args &&...args)
    : BaseLocalDynamics<DynamicsIdentifier>(identifier),
      dv_variable_(this->particles_->template getVariableByName<DataType>(variable_name)),
      assignment_method_(this->particles_, std::forward<Args>(args)...)
{
    if (dv_variable_->getWidth() != assignment_method_.Width())
    {
        std::cout << "\n Error: the width of variable '" << variable_name
                  << "' does not match the assignment method!" << std::endl;
        exit(1);
    }
}
//=================================================================================================//
template <class DynamicsIdentifier, typename AssignmentFunctionType>
template <class ExecutionPolicy, class EncloserType>
MultiEntryVariableAssignment<DynamicsIdentifier, AssignmentFunctionType>::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : multi_entry_(encloser.dv_variable_->DelegatedMultiEntryView(ex_policy)),
      assign_(ex_policy, encloser.assignment_method_){};
//=================================================================================================//
template <class DynamicsIdentifier, typename AssignmentFunctionType>
void MultiEntryVariableAssignment<DynamicsIdentifier, AssignmentFunctionType>::UpdateKernel::update(
    UnsignedInt index_i, Real dt)
{
    for (UnsignedInt entry = 0; entry < multi_entry_.Width(); ++entry)
    {
        multi_entry_[index_i][entry] = assign_(index_i, entry);
    }
}
//=================================================================================================//
} // namespace SPH
#endif // GENERAL_ASSIGNMENT_HPP
