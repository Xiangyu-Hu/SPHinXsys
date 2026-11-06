#ifndef SUBSCALE_BUOYANCY_HPP
#define SUBSCALE_BUOYANCY_HPP

#include "subscale_buoyancy.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
template <typename... Parameters>
SubscaleBuoyancy<Contact<WithUpdate, Parameters...>>::SubscaleBuoyancy(
    Contact<Parameters...> &water_air_relation, const Gravity &gravity, Real threshold)
    : BaseDynamicsType(water_air_relation), ForcePriorCK(this->particles_, "SubscaleBuoyancy"),
      gravity_(gravity), threshold_by_dimensions_(threshold * Dimensions),
      dv_mass(this->particles_->template getVariableByName<Real>("Mass")),
      dv_isolation_factor_(nullptr), // from contact body
      dv_subscale_buoyancy_(this->getCurrentForce())
{
    dv_isolation_factor_ = this->contact_particles_->template getVariableByName<Real>("IsolationFactor");
}
//=================================================================================================//
template <typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
SubscaleBuoyancy<Contact<WithUpdate, Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : BaseDynamicsType::InteractKernel(ex_policy, encloser),
      gravity_(encloser.gravity_),
      threshold_by_dimensions_(encloser.threshold_by_dimensions_),
      Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
      mass_(encloser.dv_mass->DelegatedData(ex_policy)),
      isolation_factor_(encloser.dv_isolation_factor_->DelegatedData(ex_policy)),
      subscale_buoyancy_(encloser.dv_subscale_buoyancy_->DelegatedData(ex_policy)) {};
//=================================================================================================//
template <typename... Parameters>
void SubscaleBuoyancy<Contact<WithUpdate, Parameters...>>::InteractKernel::interact(size_t index_i, Real dt)
{
    Vecd buoyancy_i = Vecd::Zero();
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        Real W_ij_Vi = this->W_ij(index_i, index_j) * Vol_[index_i];

        // Reverse gravity when air particle isolation factor less than threshold
        // Additional mass is a assume for only fraction of exclused volume
        // so that the summation is only contributed on water particles
        // Also, the original gravity should be canceled
        Vecd buoyancy_j = -SMAX(1.0 - isolation_factor_[index_j] / threshold_by_dimensions_, 0.0) *
                          2.0 * gravity_.InducedAcceleration(Vecd::Zero(), 0.0);
        buoyancy_i += W_ij_Vi * buoyancy_j;
    }

    subscale_buoyancy_[index_i] = buoyancy_i * mass_[index_i];
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
#endif // SUBSCALE_BUOYANCY_HPP
