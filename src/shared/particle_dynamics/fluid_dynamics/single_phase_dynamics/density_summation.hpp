#pragma once

#include "density_summation.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
template <class DataDelegationType>
template <class BaseRelationType>
BaseDensitySummation<DataDelegationType>::BaseDensitySummation(BaseRelationType &base_relation)
    : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation),
      rho_(this->particles_->rho_), mass_(this->particles_->mass_),
      rho0_(this->sph_body_.base_material_->ReferenceDensity()),
      inv_sigma0_(1.0 / this->sph_body_.sph_adaptation_->LatticeNumberDensity())
{
    this->particles_->registerVariable(rho_sum_, "DensitySummation");
}
//=================================================================================================//
template <class DensitySummationType>
void DensitySummationFreeSurface<DensitySummationType>::update(size_t index_i, Real dt)
{
    this->rho_[index_i] = ReinitializedDensity(this->rho_sum_[index_i], this->rho0_);
}
//=================================================================================================//
template <class DensitySummationFreeSurfaceType>
template <typename... ConstructorArgs>
DensitySummationFreeStream<DensitySummationFreeSurfaceType>::
    DensitySummationFreeStream(ConstructorArgs &&...args)
    : DensitySummationFreeSurfaceType(std::forward<ConstructorArgs>(args)...),
      indicator_(*this->particles_->template getVariableByName<int>("Indicator")){};
//=================================================================================================//
template <class DensitySummationFreeSurfaceType>
void DensitySummationFreeStream<DensitySummationFreeSurfaceType>::update(size_t index_i, Real dt)
{
    if (this->rho_sum_[index_i] < this->rho0_ && isNearFreeSurface(index_i))
    {
        this->rho_[index_i] = this->ReinitializedDensity(this->rho_sum_[index_i], this->rho0_);
    }
    else
    {
        this->rho_[index_i] = this->rho_sum_[index_i];
    }
}
//=================================================================================================//
template <class DensitySummationFreeSurfaceType>
bool DensitySummationFreeStream<DensitySummationFreeSurfaceType>::isNearFreeSurface(size_t index_i)
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
