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
template <typename... SummationType>
DensitySummation<Inner<FreeStream, SummationType...>>::
    DensitySummation(BaseInnerRelation &inner_relation, Real free_stream_pressure)
    : DensitySummation<Inner<SummationType...>>(inner_relation),
      indicator_(*this->particles_->template getVariableByName<int>("Indicator"))
{
    free_stream_rho_ = DynamicCast<Fluid>(this, this->particles_->getBaseMaterial())
                           .DensityFromPressure(free_stream_pressure);
};
//=================================================================================================//
template <typename... SummationType>
void DensitySummation<Inner<FreeStream, SummationType...>>::update(size_t index_i, Real dt)
{
    if (this->rho_sum_[index_i] < free_stream_rho_ && isNearFreeSurface(index_i))
    {
        this->rho_[index_i] = reinitializeDensity(this->rho_sum_[index_i], free_stream_rho_, this->rho_[index_i]);
    }
    else
    {
        this->rho_[index_i] = this->rho_sum_[index_i];
    }
}
//=================================================================================================//
template <typename... SummationType>
Real DensitySummation<Inner<FreeStream, SummationType...>>::reinitializeDensity(Real rho_sum, Real rho0, Real rho)
{
    return rho_sum + SMAX(0.0, (rho - rho_sum)) * rho0 / rho;
}
//=================================================================================================//
template <typename... SummationType>
bool DensitySummation<Inner<FreeStream, SummationType...>>::isNearFreeSurface(size_t index_i)
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
