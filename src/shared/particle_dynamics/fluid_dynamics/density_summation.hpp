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
      rho_(this->particles_->template getVariableDataByName<Real>("Density")),
      mass_(this->particles_->template getVariableDataByName<Real>("Mass")),
      rho_sum_(this->particles_->template registerStateVariable<Real>("DensitySummation")),
      Vol_(this->particles_->template getVariableDataByName<Real>("VolumetricMeasure")),
      rho0_(this->sph_body_.getBaseMaterial().ReferenceDensity()),
      inv_sigma0_(1.0 / this->sph_body_.getSPHAdaptation().LatticeNumberDensity()),
      W0_(this->sph_body_.getSPHAdaptation().getKernel()->W0(ZeroVecd)) {}
//=================================================================================================//
template <typename... SummationType>
template <typename... Args>
DensitySummation<Inner<FreeSurface, SummationType...>>::DensitySummation(Args &&...args)
    : DensitySummation<Inner<SummationType...>>(std::forward<Args>(args)...) {}
//=================================================================================================//
template <typename... SummationType>
void DensitySummation<Inner<FreeSurface, SummationType...>>::update(size_t index_i, Real dt)
{
    this->rho_[index_i] = SMAX(this->rho_sum_[index_i], this->rho0_);
}
//=================================================================================================//
template <typename NearSurfaceType, typename... SummationType>
template <typename... Args>
DensitySummation<Inner<NearSurfaceType, SummationType...>>::DensitySummation(Args &&...args)
    : DensitySummation<Inner<SummationType...>>(std::forward<Args>(args)...),
      indicator_(this->particles_->template getVariableDataByName<int>("Indicator")){};
//=================================================================================================//
template <typename NearSurfaceType, typename... SummationType>
void DensitySummation<Inner<NearSurfaceType, SummationType...>>::update(size_t index_i, Real dt)
{
    this->rho_[index_i] =
        isNearFreeSurface(index_i)
            ? near_surface_rho_(this->rho_sum_[index_i], this->rho0_, this->rho_[index_i])
            : this->rho_sum_[index_i];
}
//=================================================================================================//
template <typename NearSurfaceType, typename... SummationType>
bool DensitySummation<Inner<NearSurfaceType, SummationType...>>::isNearFreeSurface(size_t index_i)
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
