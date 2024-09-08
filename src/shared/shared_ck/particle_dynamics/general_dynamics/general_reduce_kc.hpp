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
    : variable_(encloser.dv_variable_->template DelegatedDataField(ex_policy)) {}
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
TotalKineticEnergyCK::ReduceKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : mass_(encloser.dv_mass_->DelegatedDataField(ex_policy)),
      vel_(encloser.dv_vel_->DelegatedDataField(ex_policy)) {}
//=================================================================================================//
template <class ExecutionPolicy, class EncloserType>
TotalMechanicalEnergyCK::ReduceKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : TotalKineticEnergyCK::ReduceKernel(ex_policy, encloser),
      gravity_(encloser.gravity_),
      pos_(encloser.dv_pos_->DelegatedDataField(ex_policy)) {}
//=================================================================================================//
} // namespace SPH
#endif // GENERAL_REDUCE_CK_HPP
