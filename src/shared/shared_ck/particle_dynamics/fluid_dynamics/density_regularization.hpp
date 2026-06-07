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
CompressionSummation<Base, RelationType<Parameters...>>::CompressionSummation(
    DynamicsIdentifier &identifier)
    : Interaction<RelationType<Parameters...>>(identifier),
      dv_Vol_ref_(this->particles_->template registerStateVariableFrom<Real>(
          "VolumetricMeasureRef", "VolumetricMeasure")),
      dv_compression_sum_(this->particles_->template registerStateVariable<Real>(
          "CompressionSummation"))
{
    this->particles_->template addEvolvingVariable<Real>(dv_Vol_ref_);
}
//=================================================================================================//
template <typename... Parameters>
template <class DynamicsIdentifier>
CompressionSummation<Inner<Parameters...>>::CompressionSummation(DynamicsIdentifier &identifier)
    : BaseInteraction(identifier) {}
//=================================================================================================//
template <typename... Parameters>
template <class ExecutionPolicy, class Encloser>
CompressionSummation<Inner<Parameters...>>::InteractKernel::InteractKernel(
    const ExecutionPolicy &ex_policy, Encloser &encloser)
    : BaseInteraction::InteractKernel(ex_policy, encloser),
      zero_(Vecd::Zero()), Vol_ref_(encloser.dv_Vol_ref_->DelegatedDataView(ex_policy)),
      compression_sum_(encloser.dv_compression_sum_->DelegatedDataView(ex_policy)) {}
//=================================================================================================//
template <typename... Parameters>
void CompressionSummation<Inner<Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    Real sigma = this->W0(index_i, zero_) * Vol_ref_[index_i];
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        sigma += this->W_ij(index_i, index_j) * Vol_ref_[index_j];
    }
    compression_sum_[index_i] = sigma;
}
//=================================================================================================//
template <typename... Parameters>
template <class DynamicsIdentifier>
CompressionSummation<Contact<Parameters...>>::CompressionSummation(DynamicsIdentifier &identifier)
    : BaseInteraction(identifier)
{
    for (size_t k = 0; k != this->contact_particles_.size(); ++k)
    {
        BaseParticles *contact_particles = this->contact_particles_[k];
        dv_contact_Vol_ref_.push_back(
            contact_particles->template registerStateVariableFrom<Real>(
                "VolumetricMeasureRef", "VolumetricMeasure"));
    }
}
//=================================================================================================//
template <typename... Parameters>
template <class ExecutionPolicy, class Encloser>
CompressionSummation<Contact<Parameters...>>::InteractKernel::InteractKernel(
    const ExecutionPolicy &ex_policy, Encloser &encloser, size_t contact_index)
    : BaseInteraction::InteractKernel(ex_policy, encloser, contact_index),
      compression_sum_(encloser.dv_compression_sum_->DelegatedDataView(ex_policy)),
      contact_Vol_ref_(encloser.dv_contact_Vol_ref_[contact_index]->DelegatedDataView(ex_policy)) {}
//=================================================================================================//
template <typename... Parameters>
void CompressionSummation<Contact<Parameters...>>::
    InteractKernel::interact(size_t index_i, Real dt)
{
    Real sigma(0);
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        UnsignedInt index_j = this->neighbor_index_[n];
        sigma += this->W_ij(index_i, index_j) * contact_Vol_ref_[index_j];
    }
    compression_sum_[index_i] = sigma;
}
//=================================================================================================//
template <class DynamicsIdentifier, class FluidType, class FlowType, typename... ParticleScopes>
DensityRegularization<DynamicsIdentifier, FluidType, FlowType, ParticleScopes...>::
    DensityRegularization(DynamicsIdentifier &identifier)
    : BaseLocalDynamics<DynamicsIdentifier>(identifier),
      fluid_(DynamicCast<FluidType>(this, this->sph_body_->getMatterMaterial())),
      dv_rho_(this->particles_->template getVariableByName<Real>("Density")),
      dv_compression_sum_(this->particles_->template getVariableByName<Real>("CompressionSummation")),
      regularization_method_(this->particles_), within_scope_method_(this->particles_)
{
    static_assert(std::is_base_of<WithinScope, ParticleScopeTypeCK<ParticleScopes...>>::value,
                  "WithinScope is not the base of ParticleScope!");
}
//=================================================================================================//
template <class DynamicsIdentifier, class FluidType, class FlowType, typename... ParticleScopes>
template <class ExecutionPolicy, class Encloser>
DensityRegularization<DynamicsIdentifier, FluidType, FlowType, ParticleScopes...>::UpdateKernel::
    UpdateKernel(const ExecutionPolicy &ex_policy, Encloser &encloser)
    : eos_(ex_policy, encloser.fluid_), rho_(encloser.dv_rho_->DelegatedDataView(ex_policy)),
      compression_sum_(encloser.dv_compression_sum_->DelegatedDataView(ex_policy)),
      regularization_(ex_policy, encloser.regularization_method_),
      particle_scope_(ex_policy, encloser.within_scope_method_) {}
//=================================================================================================//
template <class DynamicsIdentifier, class FluidType, class FlowType, typename... ParticleScopes>
void DensityRegularization<DynamicsIdentifier, FluidType, FlowType, ParticleScopes...>::
    UpdateKernel::update(size_t index_i, Real dt)
{
    if (this->particle_scope_(index_i))
    {
        rho_[index_i] = regularization_(
            index_i, compression_sum_[index_i], eos_.getReferenceDensity(index_i));
    }
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
#endif // DENSITY_REGULARIZATION_HPP
