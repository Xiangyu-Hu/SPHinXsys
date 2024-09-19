#pragma once

#include "kernel_summation.h"

namespace SPH
{
//=================================================================================================//
template <class DataDelegationType>
template <class BaseRelationType>
NablaWV<DataDelegationType>::NablaWV(BaseRelationType &base_relation)
    : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation),
      kernel_sum_(this->particles_->template registerStateVariable<Vecd>("KernelSummation")) {}
//=================================================================================================//
} // namespace SPH