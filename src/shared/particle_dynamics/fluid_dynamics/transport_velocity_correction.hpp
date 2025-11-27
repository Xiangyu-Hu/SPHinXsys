#pragma once

#include "transport_velocity_correction.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
template <class DataDelegationType, class KernelCorrectionType, class ParticleScope>
template <class BaseRelationType>
TransportVelocityCorrection<Base, DataDelegationType, KernelCorrectionType, ParticleScope>::
    TransportVelocityCorrection(BaseRelationType &base_relation)
    : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation),
      kernel_gradient_integral_(this->particles_->template registerStateVariableData<Vecd>("KernelGradientIntegral")),
      kernel_correction_(this->particles_), within_scope_(this->particles_)
{
    static_assert(std::is_base_of<WithinScope, ParticleScope>::value,
                  "WithinScope is not the base of ParticleScope!");
}
//=================================================================================================//
template <class AdaptationType, class LimiterType, typename... CommonControlTypes>
TransportVelocityCorrection<Inner<AdaptationType, LimiterType>, CommonControlTypes...>::
    TransportVelocityCorrection(BaseInnerRelation &inner_relation, Real coefficient)
    : TransportVelocityCorrection<Base, DataDelegateInner, CommonControlTypes...>(inner_relation),
      h_ref_(this->getSPHAdaptation().ReferenceSmoothingLength()),
      correction_scaling_(coefficient * h_ref_ * h_ref_),
      Vol_(this->particles_->template getVariableDataByName<Real>("VolumetricMeasure")),
      pos_(this->particles_->template getVariableDataByName<Vecd>("Position")),
      h_ratio_(DynamicCast<AdaptationType>(this, this->getSPHAdaptation())),
      limiter_(h_ref_ * h_ref_)
{
    static_assert(std::is_base_of<Limiter, LimiterType>::value,
                  "Limiter is not the base of LimiterType!");
}
//=================================================================================================//
template <class AdaptationType, class LimiterType, typename... CommonControlTypes>
void TransportVelocityCorrection<Inner<AdaptationType, LimiterType>, CommonControlTypes...>::
    interaction(size_t index_i, Real dt)
{
    if (this->within_scope_(index_i))
    {
        Vecd inconsistency = Vecd::Zero();
        const Neighborhood &inner_neighborhood = this->inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];
            // acceleration for transport velocity
            inconsistency -= (this->kernel_correction_(index_i) + this->kernel_correction_(index_j, index_i)) *
                             inner_neighborhood.dW_ij_[n] * this->Vol_[index_j] * inner_neighborhood.e_ij_[n];
        }
        this->kernel_gradient_integral_[index_i] = inconsistency;
    }
}
//=================================================================================================//
template <class AdaptationType, class LimiterType, typename... CommonControlTypes>
void TransportVelocityCorrection<Inner<AdaptationType, LimiterType>, CommonControlTypes...>::
    update(size_t index_i, Real dt)
{
    if (this->within_scope_(index_i))
    {
        Real inv_h_ratio = 1.0 / h_ratio_(index_i);
        Real squared_norm = this->kernel_gradient_integral_[index_i].squaredNorm();
        pos_[index_i] += correction_scaling_ * limiter_(squared_norm) *
                         this->kernel_gradient_integral_[index_i] * inv_h_ratio * inv_h_ratio;
    }
}
//=================================================================================================//
template <typename... CommonControlTypes>
TransportVelocityCorrection<Contact<Boundary>, CommonControlTypes...>::
    TransportVelocityCorrection(BaseContactRelation &contact_relation)
    : TransportVelocityCorrection<Base, DataDelegateContact, CommonControlTypes...>(contact_relation)
{
    for (size_t k = 0; k != this->contact_particles_.size(); ++k)
    {
        wall_Vol_.push_back(this->contact_particles_[k]->template getVariableDataByName<Real>("VolumetricMeasure"));
    }
};
//=================================================================================================//
template <typename... CommonControlTypes>
void TransportVelocityCorrection<Contact<Boundary>, CommonControlTypes...>::
    interaction(size_t index_i, Real dt)
{
    if (this->within_scope_(index_i))
    {
        Vecd inconsistency = Vecd::Zero();
        for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
        {
            Real *wall_Vol_k = wall_Vol_[k];
            Neighborhood &contact_neighborhood = (*this->contact_configuration_[k])[index_i];
            for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
            {
                size_t index_j = contact_neighborhood.j_[n];
                // acceleration for transport velocity
                inconsistency -= 2.0 * this->kernel_correction_(index_i) * contact_neighborhood.dW_ij_[n] *
                                 wall_Vol_k[index_j] * contact_neighborhood.e_ij_[n];
            }
        }
        this->kernel_gradient_integral_[index_i] += inconsistency;
    }
}
//=================================================================================================//
template <class KernelCorrectionType, typename... CommonControlTypes>
TransportVelocityCorrection<Contact<>, KernelCorrectionType, CommonControlTypes...>::
    TransportVelocityCorrection(BaseContactRelation &contact_relation)
    : TransportVelocityCorrection<Base, DataDelegateContact, KernelCorrectionType, CommonControlTypes...>(
          contact_relation)
{
    for (size_t k = 0; k != this->contact_particles_.size(); ++k)
    {
        contact_kernel_corrections_.push_back(KernelCorrectionType(this->contact_particles_[k]));
        contact_Vol_.push_back(this->contact_particles_[k]->template getVariableDataByName<Real>("VolumetricMeasure"));
    }
}
//=================================================================================================//
template <class KernelCorrectionType, typename... CommonControlTypes>
void TransportVelocityCorrection<Contact<>, KernelCorrectionType, CommonControlTypes...>::
    interaction(size_t index_i, Real dt)
{
    if (this->within_scope_(index_i))
    {
        Vecd inconsistency = Vecd::Zero();
        for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
        {
            Real *Vol_k = this->contact_Vol_[k];
            Neighborhood &contact_neighborhood = (*this->contact_configuration_[k])[index_i];
            KernelCorrectionType &kernel_correction_k = this->contact_kernel_corrections_[k];
            for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
            {
                size_t index_j = contact_neighborhood.j_[n];
                // acceleration for transport velocity
                inconsistency -= (this->kernel_correction_(index_i) + kernel_correction_k(index_j)) *
                                 contact_neighborhood.dW_ij_[n] * Vol_k[index_j] * contact_neighborhood.e_ij_[n];
            }
        }
        this->kernel_gradient_integral_[index_i] += inconsistency;
    }
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
