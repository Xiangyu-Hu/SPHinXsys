#ifndef GENERAL_REDUCE_CK_HPP
#define GENERAL_REDUCE_CK_HPP

#include "general_reduce_ck.h"

namespace SPH
{
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
TotalKineticEnergyCK::ReduceKernel::
    ReduceKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : mass_(encloser.dv_mass_->DelegatedData(ex_policy)),
      vel_(encloser.dv_vel_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
TotalMechanicalEnergyCK::ReduceKernel::
    ReduceKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : TotalKineticEnergyCK::ReduceKernel(ex_policy, encloser),
      gravity_(encloser.gravity_),
      pos_(encloser.dv_pos_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <typename DataType, class DynamicsIdentifier>
QuantitySum<DataType, DynamicsIdentifier>::
    QuantitySum(DynamicsIdentifier &identifier, const std::string &variable_name)
    : BaseLocalDynamicsReduce<ReduceSum<DataType>, DynamicsIdentifier>(identifier),
      dv_variable_(this->particles_->template getVariableByName<DataType>(variable_name))
{
    this->quantity_name_ = "Total" + variable_name;
}
//=================================================================================================//
template <typename DataType, class DynamicsIdentifier>
template <class ExecutionPolicy, class EncloserType>
QuantitySum<DataType, DynamicsIdentifier>::ReduceKernel::
    ReduceKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : variable_(encloser.dv_variable_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <typename DataType, class DynamicsIdentifier>
QuantityAverage<DataType, DynamicsIdentifier>::
    QuantityAverage(DynamicsIdentifier &identifier, const std::string &variable_name)
    : BaseDynamicsType(identifier),
      dv_variable_(this->particles_->template getVariableByName<DataType>(variable_name))
{
    this->quantity_name_ = "Average" + variable_name;
}
//=================================================================================================//
template <typename DataType, class DynamicsIdentifier>
template <class ExecutionPolicy, class EncloserType>
QuantityAverage<DataType, DynamicsIdentifier>::ReduceKernel::
    ReduceKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : variable_(encloser.dv_variable_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class DynamicsIdentifier>
UpperFrontInAxisDirectionCK<DynamicsIdentifier>::UpperFrontInAxisDirectionCK(
    DynamicsIdentifier &identifier, const std::string &name, int axis)
    : BaseLocalDynamicsReduce<ReduceMax, DynamicsIdentifier>(identifier),
      axis_(axis), dv_pos_(this->particles_->template getVariableByName<Vecd>("Position"))
{
    this->quantity_name_ = name;
}
//=================================================================================================//
template <class DynamicsIdentifier>
template <class ExecutionPolicy, class EncloserType>
UpperFrontInAxisDirectionCK<DynamicsIdentifier>::ReduceKernel::
    ReduceKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : axis_(encloser.axis_), pos_(encloser.dv_pos_->DelegatedData(ex_policy)) {}
//=================================================================================================//
} // namespace SPH
#endif // GENERAL_REDUCE_CK_HPP
