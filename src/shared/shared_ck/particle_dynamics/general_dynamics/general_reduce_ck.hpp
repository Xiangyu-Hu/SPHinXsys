#ifndef GENERAL_REDUCE_CK_HPP
#define GENERAL_REDUCE_CK_HPP

#include "general_reduce_ck.h"

namespace SPH
{
//=================================================================================================//
template <typename DataType, typename NormType, class DynamicsIdentifier>
VariableNormCK<DataType, NormType, DynamicsIdentifier>::
    VariableNormCK(DynamicsIdentifier &identifier, const std::string &variable_name)
    : BaseLocalDynamicsReduce<NormType, DynamicsIdentifier>(identifier),
      dv_variable_(this->particles_->template getVariableByName<DataType>(variable_name)) {}
//=================================================================================================//
template <typename DataType, typename NormType, class DynamicsIdentifier>
template <class ExecutionPolicy, class EncloserType>
VariableNormCK<DataType, NormType, DynamicsIdentifier>::ReduceKernel::
    ReduceKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : variable_(encloser.dv_variable_->DelegatedData(ex_policy)) {}
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
template <class ReduceSumType>
template <class DynamicsIdentifier>
AverageCK<ReduceSumType>::AverageCK(DynamicsIdentifier &identifier, const std::string &variable_name)
    : ReduceSumType(identifier, variable_name),
      sv_total_sample_size_(identifier.getTotalSize())
{
    this->quantity_name_ = "Average" + variable_name;
}
//=================================================================================================//
} // namespace SPH
#endif // GENERAL_REDUCE_CK_HPP
