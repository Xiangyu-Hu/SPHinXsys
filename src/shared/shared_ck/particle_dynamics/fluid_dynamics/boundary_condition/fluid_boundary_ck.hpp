#ifndef FLUID_BOUNDARY_CK_HPP
#define FLUID_BOUNDARY_CK_HPP

#include "fluid_boundary_ck.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
template <class AlignedBoxPartType, class ConditionFunction, class TransformFunction>
InflowConditionCK<AlignedBoxPartType, ConditionFunction, TransformFunction>::
    InflowConditionCK(AlignedBoxPartType &aligned_box_part)
    : BaseLocalDynamics<AlignedBoxPartType>(aligned_box_part),
      sv_aligned_box_(this->particles_->addUniqueSingularVariable<AlignedBox>(
          "AlignedBox", aligned_box_part.getAlignedBox())),
      condition_function_(this->particles_),
      dv_pos_(this->particles_->getStateVariableByName<Vecd>("Position")) {}
//=================================================================================================//
template <class AlignedBoxPartType, class ConditionFunction, class TransformFunction>
template <class ExecutionPolicy, class EncloserType>
InflowConditionCK<AlignedBoxPartType, ConditionFunction, TransformFunction>::
    UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser),
    aligned_box_(encloser.sv_aligned_box_->DelegatedData(ex_policy)),
    condition_(ex_policy, encloser.condition_function_){}
//=================================================================================================//
template <class AlignedBoxPartType, class ConditionFunction, class TransformFunction>

void InflowConditionCK<AlignedBoxPartType, ConditionFunction, TransformFunction>::
    UpdateKernel::update(size_t index_i, Real dt)
{
    condition_(aligned_box_, index_i);
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
#endif // FLUID_BOUNDARY_CK_HPP
