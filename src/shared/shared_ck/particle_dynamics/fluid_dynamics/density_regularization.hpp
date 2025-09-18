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
DensityRegularization<Base, RelationType<Parameters...>>::
    DensityRegularization(DynamicsIdentifier &identifier)
    : Interaction<RelationType<Parameters...>>(identifier),
      dv_rho_(this->particles_->template getVariableByName<Real>("Density")),
      dv_mass_(this->particles_->template getVariableByName<Real>("Mass")),
      dv_rho_sum_(this->particles_->template registerStateVariable<Real>("DensitySummation")),
      rho0_(this->sph_body_->getBaseMaterial().ReferenceDensity()),
      inv_sigma0_(1.0 / this->sph_adaptation_->LatticeNumberDensity()) {}
//=================================================================================================//
template <template <typename...> class RelationType, typename... Parameters>
template <class ExecutionPolicy, typename... Args>
DensityRegularization<Base, RelationType<Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy,
                   DensityRegularization<Base, RelationType<Parameters...>> &encloser,
                   Args &&...args)
    : Interaction<RelationType<Parameters...>>::
          InteractKernel(ex_policy, encloser, std::forward<Args>(args)...),
      rho_(encloser.dv_rho_->DelegatedData(ex_policy)),
      mass_(encloser.dv_mass_->DelegatedData(ex_policy)),
      rho_sum_(encloser.dv_rho_sum_->DelegatedData(ex_policy)),
      Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
      rho0_(encloser.rho0_), inv_sigma0_(encloser.inv_sigma0_) {}
//=================================================================================================//
template <typename RegularizationType, typename ParticleScopeType, typename... Parameters>
DensityRegularization<Inner<WithUpdate, RegularizationType, ParticleScopeType, Parameters...>>::
    DensityRegularization(Inner<Parameters...> &inner_relation)
    : DensityRegularization<Base, Inner<Parameters...>>(inner_relation),
      regularization_method_(this->particles_),
      within_scope_method_(this->particles_)
{
}
//=================================================================================================//
template <typename RegularizationType, typename ParticleScopeType, typename... Parameters>
template <class ExecutionPolicy>
DensityRegularization<Inner<WithUpdate, RegularizationType, ParticleScopeType, Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy,
                   DensityRegularization<Inner<WithUpdate, RegularizationType, ParticleScopeType, Parameters...>> &encloser)
    : DensityRegularization<Base, Inner<Parameters...>>::InteractKernel(ex_policy, encloser),
      W0_(this->W(ZeroData<Vecd>::value)) {}
//=================================================================================================//
template <typename RegularizationType, typename ParticleScopeType, typename... Parameters>
void DensityRegularization<Inner<WithUpdate, RegularizationType, ParticleScopeType, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    Real sigma = W0_;
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
        sigma += this->W_ij(index_i, this->neighbor_index_[n]);

    this->rho_sum_[index_i] = sigma * this->rho0_ * this->inv_sigma0_;
}
//=================================================================================================//
template <typename RegularizationType, typename ParticleScopeType, typename... Parameters>
template <class ExecutionPolicy>
DensityRegularization<Inner<WithUpdate, RegularizationType, ParticleScopeType, Parameters...>>::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy,
                 DensityRegularization<Inner<WithUpdate, RegularizationType, ParticleScopeType, Parameters...>> &encloser)
    : DensityRegularization<Base, Inner<Parameters...>>::InteractKernel(ex_policy, encloser),
      regularization_(ex_policy, encloser.regularization_method_, *this),
      particle_scope_(ex_policy, encloser.within_scope_method_, *this)
{
}
//=================================================================================================//
template <typename RegularizationType, typename ParticleScopeType, typename... Parameters>
void DensityRegularization<Inner<WithUpdate, RegularizationType, ParticleScopeType, Parameters...>>::
    UpdateKernel::update(size_t index_i, Real dt)
{
    if (this->particle_scope_(index_i))
        this->rho_[index_i] = regularization_(index_i, this->rho_sum_[index_i]);
}
//=================================================================================================//
template <typename... Parameters>
DensityRegularization<Contact<Parameters...>>::
    DensityRegularization(Contact<Parameters...> &contact_relation)
    : DensityRegularization<Base, Contact<Parameters...>>(contact_relation)
{
    for (size_t k = 0; k != this->contact_particles_.size(); ++k)
    {
        Real rho0_k = this->contact_bodies_[k]->getBaseMaterial().ReferenceDensity();
        contact_inv_rho0_.push_back(1.0 / rho0_k);
        dv_contact_mass_.push_back(this->contact_particles_[k]->template getVariableByName<Real>("Mass"));
    }
}
//=================================================================================================//
template <typename... Parameters>
template <class ExecutionPolicy>
DensityRegularization<Contact<Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy,
                   DensityRegularization<Contact<Parameters...>> &encloser,
                   size_t contact_index)
    : DensityRegularization<Base, Contact<Parameters...>>::
          InteractKernel(ex_policy, encloser, contact_index),
      contact_inv_rho0_k_(encloser.contact_inv_rho0_[contact_index]),
      contact_mass_k_(encloser.dv_contact_mass_[contact_index]->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <typename... Parameters>
void DensityRegularization<Contact<Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    Real sigma(0);
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        sigma += this->W_ij(index_i, index_j) * contact_inv_rho0_k_ * contact_mass_k_[index_j];
    }
    this->rho_sum_[index_i] += sigma * this->rho0_ * this->rho0_ *
                               this->inv_sigma0_ / this->mass_[index_i];
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
#endif // DENSITY_REGULARIZATION_HPP
