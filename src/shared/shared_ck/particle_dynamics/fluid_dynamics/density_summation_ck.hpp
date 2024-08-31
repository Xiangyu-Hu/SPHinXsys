#ifndef DENSITY_SUMMATION_CK_HPP
#define DENSITY_SUMMATION_CK_HPP

#include "density_summation_ck.h"

#include "base_particles.hpp"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
template <class LocalDynamicsType>
template <class DynamicsIdentifier>
DensitySummationCK<LocalDynamicsType>::DensitySummationCK(DynamicsIdentifier &identifier)
    : LocalDynamicsType(identifier),
      dv_rho_(this->particles_->template getVariableByName<Real>("Density")),
      dv_mass_(this->particles_->template getVariableByName<Real>("Mass")),
      dv_rho_sum_(this->particles_->template getVariableByName<Real>("DensitySummation")),
      dv_Vol_(this->particles_->template getVariableByName<Real>("VolumetricMeasure")),
      rho0_(this->sph_body_.getBaseMaterial().ReferenceDensity()),
      inv_sigma0_(1.0 / this->sph_adaptation_->LatticeNumberDensity()) {}
//=================================================================================================//
template <class LocalDynamicsType>
template <class ExecutionPolicy>
DensitySummationCK<LocalDynamicsType>::ComputingKernel::
    ComputingKernel(const ExecutionPolicy &ex_policy,
                    DensitySummationCK<LocalDynamicsType> &encloser)
    : LocalDynamicsType::ComputingKernel(ex_policy, encloser),
      rho_(encloser.dv_rho_->DelegatedDataField(ex_policy)),
      mass_(encloser.dv_mass_->DelegatedDataField(ex_policy)),
      rho_sum_(encloser.dv_rho_sum_->DelegatedDataField(ex_policy)),
      Vol_(encloser.dv_Vol_->DelegatedDataField(ex_policy)),
      rho0_(encloser.rho0_), inv_sigma0_(encloser.inv_sigma0_) {}
//=================================================================================================//
template <typename... RegularizationType, typename... Parameters>
DensitySummationCK<Inner<Regularization<RegularizationType...>, Parameters...>>::
    DensitySummationCK(InnerRelation &inner_relation)
    : DensitySummationCK<LocalInteractionDynamics<Inner<Parameters...>>>(inner_relation) {}
//=================================================================================================//
template <typename... RegularizationType, typename... Parameters>
template <class ExecutionPolicy>
DensitySummationCK<Inner<Regularization<RegularizationType...>, Parameters...>>::
    ComputingKernel::ComputingKernel(const ExecutionPolicy &ex_policy,
                                     DensitySummationCK<Inner<Parameters...>> &encloser)
    : DensitySummationCK<LocalInteractionDynamics<Inner<Parameters...>>>::
          ComputingKernel(ex_policy, encloser),
      W0_(this->kernel_.W(ZeroData<Vecd>::value)), regularization_(*this) {}
//=================================================================================================//
template <typename... RegularizationType, typename... Parameters>
void DensitySummationCK<Inner<Regularization<RegularizationType...>, Parameters...>>::
    ComputingKernel::interaction(size_t index_i, Real dt)
{
    Real sigma = W0_;
    for (UnsignedInt n = this->FistNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
        sigma += this->W_ij(index_i, this->NeighborIndex(n));

    this->rho_sum_[index_i] = sigma * this->rho0_ * this->inv_sigma0_;
}
//=================================================================================================//
template <typename... RegularizationType, typename... Parameters>
void DensitySummationCK<Inner<Regularization<RegularizationType...>, Parameters...>>::
    ComputingKernel::update(size_t index_i, Real dt)
{
    this->rho_[index_i] = regularization_(this->rho_sum_[index_i]);
    this->Vol_[index_i] = this->mass_[index_i] / this->rho_[index_i];
}
//=================================================================================================//
template <typename... Parameters>
DensitySummationCK<Contact<Parameters...>>::DensitySummationCK(ContactRelation &contact_relation)
    : DensitySummationCK<LocalInteractionDynamics<Contact<Parameters...>>>(contact_relation)
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
DensitySummationCK<Contact<Parameters...>>::ComputingKernel::
    ComputingKernel(const ExecutionPolicy &ex_policy,
                    DensitySummationCK<Contact<Parameters...>> &encloser,
                    size_t contact_index)
    : DensitySummationCK<LocalInteractionDynamics<Contact<Parameters...>>>::
          ComputingKernel(ex_policy, encloser, contact_index),
      contact_inv_rho0_k_(encloser.contact_inv_rho0_),
      contact_mass_k_(encloser.contact_mass_[contact_index]->DelegatedDataField(ex_policy)) {}
//=================================================================================================//
template <typename... Parameters>
void DensitySummationCK<Contact<Parameters...>>::ComputingKernel::
    interaction(size_t index_i, Real dt)
{
    Real sigma(0);
    for (UnsignedInt n = this->FistNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->NeighborIndex(n);
        sigma += this->W_ij(index_i, index_j) * contact_inv_rho0_k_ * contact_mass_k_[index_j];
    }
    this->rho_sum_[index_i] = sigma * this->rho0_ * this->rho0_ *
                              this->inv_sigma0_ / this->mass_[index_i];
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
#endif // DENSITY_SUMMATION_CK_HPP
