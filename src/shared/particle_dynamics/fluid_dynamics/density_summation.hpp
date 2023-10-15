#pragma once

#include "density_summation.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
template <class DataDelegationType>
template <class BaseRelationType>
DensitySummation<Base, DataDelegationType>::DensitySummation(BaseRelationType &base_relation)
    : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation),
      rho_(this->particles_->rho_), mass_(this->particles_->mass_),
      rho_sum_(*this->particles_->template registerSharedVariable<Real>("DensitySummation")),
      rho0_(this->sph_body_.base_material_->ReferenceDensity()),
      inv_sigma0_(1.0 / this->sph_body_.sph_adaptation_->LatticeNumberDensity()),
      W0_(this->sph_body_.sph_adaptation_->getKernel()->W0(ZeroVecd)) {}
//=================================================================================================//
template <class InteractionType>
template <typename... Args>
DensitySummation<FreeStream<InteractionType>>::
    DensitySummation(Args &&...args)
    : DensitySummation<InteractionType>(std::forward<Args>(args)...),
      indicator_(*this->particles_->template getVariableByName<int>("Indicator")){};
//=================================================================================================//
template <class InteractionType>
void DensitySummation<FreeStream<InteractionType>>::update(size_t index_i, Real dt)
{
    if (this->rho_sum_[index_i] < this->rho0_ && isNearFreeSurface(index_i))
    {
        this->reinitializeDensity(index_i);
    }
    else
    {
        this->assignDensity(index_i);
    }
}
//=================================================================================================//
template <class InteractionType>
bool DensitySummation<FreeStream<InteractionType>>::isNearFreeSurface(size_t index_i)
{
    bool is_near_surface = false;
    const Neighborhood &inner_neighborhood = this->inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        if (indicator_[inner_neighborhood.j_[n]] == 1)
        {
            is_near_surface = true;
            break;
        }
    }
    return is_near_surface;
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
