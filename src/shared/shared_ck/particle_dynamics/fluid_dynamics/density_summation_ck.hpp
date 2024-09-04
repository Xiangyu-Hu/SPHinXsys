#ifndef DENSITY_SUMMATION_CK_HPP
#define DENSITY_SUMMATION_CK_HPP

#include "density_summation_ck.h"

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
      dv_rho_(this->particles_->template getVariableByName<Real>("Density")),
      dv_mass_(this->particles_->template getVariableByName<Real>("Mass")),
      dv_rho_sum_(this->particles_->template getVariableByName<Real>("DensitySummation")),
      dv_Vol_(this->particles_->template getVariableByName<Real>("VolumetricMeasure")),
      rho0_(this->sph_body_.getBaseMaterial().ReferenceDensity()),
      inv_sigma0_(1.0 / this->sph_adaptation_->LatticeNumberDensity()) {}
//=================================================================================================//
template <template <typename...> class RelationType, typename... Parameters>
template <class ExecutionPolicy, typename... Args>
DensitySummationCK<Base, RelationType<Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy,
                   DensitySummationCK<Base, RelationType<Parameters...>> &encloser,
                   Args &&...args)
    : Interaction<RelationType<Parameters...>>::
          InteractKernel(ex_policy, encloser, std::forward<Args>(args)...),
      rho_(encloser.dv_rho_->DelegatedDataField(ex_policy)),
      mass_(encloser.dv_mass_->DelegatedDataField(ex_policy)),
      rho_sum_(encloser.dv_rho_sum_->DelegatedDataField(ex_policy)),
      Vol_(encloser.dv_Vol_->DelegatedDataField(ex_policy)),
      rho0_(encloser.rho0_), inv_sigma0_(encloser.inv_sigma0_) {}
//=================================================================================================//
template <typename RegularizationType, typename... Parameters>
DensitySummationCK<Inner<RegularizationType, Parameters...>>::
    DensitySummationCK(InnerRelation &inner_relation)
    : DensitySummationCK<Base, Inner<Parameters...>>(inner_relation),
      regularization_method_(this->particles_) {}
//=================================================================================================//
template <typename RegularizationType, typename... Parameters>
template <class ExecutionPolicy>
DensitySummationCK<Inner<RegularizationType, Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy,
                   DensitySummationCK<Inner<RegularizationType, Parameters...>> &encloser)
    : DensitySummationCK<Base, Inner<Parameters...>>::InteractKernel(ex_policy, encloser),
      W0_(this->kernel_.W(ZeroData<Vecd>::value)) {}
//=================================================================================================//
template <typename RegularizationType, typename... Parameters>
void DensitySummationCK<Inner<RegularizationType, Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    Real sigma = W0_;
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
        sigma += this->W_ij(index_i, this->neighbor_index_[n]);

    this->rho_sum_[index_i] = sigma * this->rho0_ * this->inv_sigma0_;
}
//=================================================================================================//
template <typename RegularizationType, typename... Parameters>
template <class ExecutionPolicy>
DensitySummationCK<Inner<RegularizationType, Parameters...>>::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy,
                 DensitySummationCK<Inner<RegularizationType, Parameters...>> &encloser)
    : DensitySummationCK<Base, Inner<Parameters...>>::InteractKernel(ex_policy, encloser),
      regularization_(ex_policy, encloser.regularization_method_, *this) {}
//=================================================================================================//
template <typename RegularizationType, typename... Parameters>
void DensitySummationCK<Inner<RegularizationType, Parameters...>>::
    UpdateKernel::update(size_t index_i, Real dt)
{
    this->rho_[index_i] = regularization_(this->rho_sum_[index_i]);
    this->Vol_[index_i] = this->mass_[index_i] / this->rho_[index_i];
}
//=================================================================================================//
template <typename... Parameters>
DensitySummationCK<Contact<Parameters...>>::DensitySummationCK(ContactRelation &contact_relation)
    : DensitySummationCK<Base, Contact<Parameters...>>(contact_relation)
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
DensitySummationCK<Contact<Parameters...>>::InteractKernel::
    InteractKernel(const ExecutionPolicy &ex_policy,
                   DensitySummationCK<Contact<Parameters...>> &encloser,
                   size_t contact_index)
    : DensitySummationCK<Base, Contact<Parameters...>>::
          InteractKernel(ex_policy, encloser, contact_index),
      contact_inv_rho0_k_(encloser.contact_inv_rho0_[contact_index]),
      contact_mass_k_(encloser.dv_contact_mass_[contact_index]->DelegatedDataField(ex_policy)) {}
//=================================================================================================//
template <typename... Parameters>
void DensitySummationCK<Contact<Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    Real sigma(0);
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        sigma += this->W_ij(index_i, index_j) * contact_inv_rho0_k_ * contact_mass_k_[index_j];
    }
    this->rho_sum_[index_i] = sigma * this->rho0_ * this->rho0_ *
                              this->inv_sigma0_ / this->mass_[index_i];
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
#endif // DENSITY_SUMMATION_CK_HPP
