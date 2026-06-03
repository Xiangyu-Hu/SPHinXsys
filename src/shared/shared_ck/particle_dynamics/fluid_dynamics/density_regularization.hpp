#ifndef DENSITY_REGULARIZATION_HPP
#define DENSITY_REGULARIZATION_HPP

#include "density_regularization.h"

#include "base_particles.hpp"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
template <template <typename...> class RelationType, typename... Parameters>
template <class DynamicsIdentifier>
DensitySummationCK<Base, RelationType<Parameters...>>::
    DensitySummationCK(DynamicsIdentifier &identifier)
    : Interaction<RelationType<Parameters...>>(identifier),
      dv_mass_(this->particles_->template getVariableByName<Real>("Mass")),
      dv_rho_sum_(this->particles_->template registerStateVariable<Real>("DensitySummation")),
      rho0_(this->sph_body_->getMatterMaterial().ReferenceDensity()) {}
//=================================================================================================//
template <typename... Parameters>
template <class DynamicsIdentifier>
DensitySummationCK<Inner<Parameters...>>::DensitySummationCK(DynamicsIdentifier &identifier)
    : BaseInteraction(identifier) {}
//=================================================================================================//
template <typename... Parameters>
template <class ExecutionPolicy, class Encloser>
DensitySummationCK<Inner<Parameters...>>::InteractKernel::InteractKernel(
    const ExecutionPolicy &ex_policy, Encloser &encloser)
    : BaseInteraction::InteractKernel(ex_policy, encloser),
      mass_(encloser.dv_mass_->DelegatedDataView(ex_policy)),
      rho_sum_(encloser.dv_rho_sum_->DelegatedDataView(ex_policy)),
      zero_(Vecd::Zero()) {}
//=================================================================================================//
template <typename... Parameters>
void DensitySummationCK<Inner<Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    Real sigma = this->W0(index_i, zero_) * mass_[index_i];
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        sigma += this->W_ij(index_i, index_j) * mass_[index_j];
    }
    rho_sum_[index_i] = sigma;
}
//=================================================================================================//
template <typename... Parameters>
template <class DynamicsIdentifier>
DensitySummationCK<Contact<Parameters...>>::DensitySummationCK(DynamicsIdentifier &identifier)
    : BaseInteraction(identifier)
{
    for (size_t k = 0; k != this->contact_particles_.size(); ++k)
    {
        Real rho0_k = this->contact_bodies_[k]->getMatterMaterial().ReferenceDensity();
        contact_inv_rho0_.push_back(1.0 / rho0_k);
        dv_contact_mass_.push_back(this->contact_particles_[k]->template getVariableByName<Real>("Mass"));
    }
}
//=================================================================================================//
template <typename... Parameters>
template <class ExecutionPolicy, class Encloser>
DensitySummationCK<Contact<Parameters...>>::InteractKernel::InteractKernel(
    const ExecutionPolicy &ex_policy, Encloser &encloser, size_t contact_index)
    : BaseInteraction::InteractKernel(ex_policy, encloser, contact_index),
      rho0_(encloser.rho0_), contact_inv_rho0_(encloser.contact_inv_rho0_[contact_index]),
      rho_sum_(encloser.dv_rho_sum_->DelegatedDataView(ex_policy)),
      contact_mass_(encloser.dv_contact_mass_[contact_index]->DelegatedDataView(ex_policy)) {}
//=================================================================================================//
template <typename... Parameters>
void DensitySummationCK<Contact<Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    Real sigma(0);
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        sigma += this->W_ij(index_i, index_j) * contact_inv_rho0_ * contact_mass_[index_j];
    }
    rho_sum_[index_i] += sigma * rho0_;
}
//=================================================================================================//
template <class DynamicsIdentifier, class FlowType, typename... ParticleScopes>
DensityRegularization<DynamicsIdentifier, FlowType, ParticleScopes...>::
    DensityRegularization(DynamicsIdentifier &identifier)
    : BaseLocalDynamics<DynamicsIdentifier>(identifier),
      rho0_(this->sph_body_->getMatterMaterial().ReferenceDensity()),
      dv_rho_(this->particles_->template getVariableByName<Real>("Density")),
      dv_rho_sum_(this->particles_->template getVariableByName<Real>("DensitySummation")),
      regularization_method_(*this), within_scope_method_(this->particles_)
{
    static_assert(std::is_base_of<WithinScope, ParticleScopeTypeCK<ParticleScopes...>>::value,
                  "WithinScope is not the base of ParticleScope!");
}
//=================================================================================================//
template <class DynamicsIdentifier, class FlowType, typename... ParticleScopes>
template <class ExecutionPolicy, class Encloser>
DensityRegularization<DynamicsIdentifier, FlowType, ParticleScopes...>::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, Encloser &encloser)
    : rho0_(encloser.rho0_),
      rho_(encloser.dv_rho_->DelegatedData(ex_policy)),
      rho_sum_(encloser.dv_rho_sum_->DelegatedData(ex_policy)),
      regularization_(ex_policy, encloser.regularization_method_),
      particle_scope_(ex_policy, encloser.within_scope_method_) {}
//=================================================================================================//
template <class DynamicsIdentifier, class FlowType, typename... ParticleScopes>
void DensityRegularization<DynamicsIdentifier, FlowType, ParticleScopes...>::
    UpdateKernel::update(size_t index_i, Real dt)
{
    if (this->particle_scope_(index_i))
        this->rho_[index_i] = regularization_(index_i, this->rho_sum_[index_i]);
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
#endif // DENSITY_REGULARIZATION_HPP
