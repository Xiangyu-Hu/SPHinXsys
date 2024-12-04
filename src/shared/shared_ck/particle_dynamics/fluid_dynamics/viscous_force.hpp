#pragma once

#include "viscous_force.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
template <template <typename...> class RelationType, typename... Parameters>
template <class BaseRelationType>
ViscousForceCK<Base, RelationType<Parameters...>>::
    ViscousForceCK(BaseRelationType &base_relation)
    : Interaction<RelationType<Parameters...>>(base_relation),
      ForcePriorCK(base_relation.getSPHBody(), "ViscousForce"),
      dv_rho_(this->particles_->template getVariableByName<Real>("Density")),
      dv_mass_(this->particles_->template getVariableByName<Real>("Mass")),
      dv_Vol_(this->particles_->template getVariableByName<Real>("VolumetricMeasure")),
      dv_vel_(this->particles_->template getVariableByName<Vecd>("Velocity")),
      dv_viscous_force_(this->particles_->template registerStateVariableOnly<Vecd>("ViscousForce")),
      smoothing_length_(this->sph_body_.sph_adaptation_->ReferenceSmoothingLength()) {}
//=================================================================================================//
template <template <typename...> class RelationType, typename... Parameters>
template <class ExecutionPolicy, typename... Args>
ViscousForceCK<Base, RelationType<Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy,
                   ViscousForceCK<Base, RelationType<Parameters...>> &encloser,
                   Args &&...args)
    : Interaction<RelationType<Parameters...>>::
          InteractKernel(ex_policy, encloser, std::forward<Args>(args)...),
      ForcePriorCK::UpdateKernel(ex_policy, encloser),
      rho_(encloser.dv_rho_->DelegatedDataField(ex_policy)),
      mass_(encloser.dv_mass_->DelegatedDataField(ex_policy)),
      Vol_(encloser.dv_Vol_->DelegatedDataField(ex_policy)),
      vel_(encloser.dv_vel_->DelegatedDataField(ex_policy)),
      viscous_force_(encloser.dv_viscous_force_->DelegatedDataField(ex_policy)),
      smoothing_length_(encloser.smoothing_length_) {}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
