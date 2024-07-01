#pragma once

#include "surface_indication.h"

namespace SPH
{
//=================================================================================================//
template <class DataDelegationType>
template <class BaseRelationType>
FreeSurfaceIndication<DataDelegationType>::FreeSurfaceIndication(BaseRelationType &base_relation)
    : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation),
      indicator_(*this->particles_->template registerSharedVariable<int>("Indicator")),
      pos_div_(*this->particles_->template registerSharedVariable<Real>("PositionDivergence")),
      Vol_(*this->particles_->template getVariableDataByName<Real>("VolumetricMeasure")),
      threshold_by_dimensions_(0.75 * Dimensions) {}
//=================================================================================================//
} // namespace SPH
