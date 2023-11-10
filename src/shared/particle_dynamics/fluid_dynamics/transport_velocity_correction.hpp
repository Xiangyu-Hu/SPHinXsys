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
      transport_acc_(*this->particles_->template registerSharedVariable<Vecd>("TransportAcceleration")),
      kernel_correction_(this->particles_), checkWithinScope(this->particles_) {}
//=================================================================================================//
template <class ResolutionType, typename... CommonControlTypes>
TransportVelocityCorrection<Inner<ResolutionType>, CommonControlTypes...>::
    TransportVelocityCorrection(BaseInnerRelation &inner_relation, Real coefficient)
    : TransportVelocityCorrection<Base, FluidDataInner, CommonControlTypes...>(inner_relation),
      correction_scaling_(coefficient * pow(this->sph_body_.sph_adaptation_->ReferenceSmoothingLength(), 2)),
      pos_(this->particles_->pos_), h_ratio_(this->particles_) {}
//=================================================================================================//
template <class ResolutionType, typename... CommonControlTypes>
void TransportVelocityCorrection<Inner<ResolutionType>, CommonControlTypes...>::
    interaction(size_t index_i, Real dt)
{
    if (this->checkWithinScope(index_i))
    {
        Vecd acceleration_trans = Vecd::Zero();
        const Neighborhood &inner_neighborhood = this->inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];
            // acceleration for transport velocity
            acceleration_trans -= (this->kernel_correction_(index_i) + this->kernel_correction_(index_j)) *
                                  inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
        }
        this->transport_acc_[index_i] = acceleration_trans;
    }
}
//=================================================================================================//
template <class ResolutionType, typename... CommonControlTypes>
void TransportVelocityCorrection<Inner<ResolutionType>, CommonControlTypes...>::
    update(size_t index_i, Real dt)
{
    if (this->checkWithinScope(index_i))
    {
        Real inv_h_ratio = 1.0 / h_ratio_(index_i);
        pos_[index_i] += correction_scaling_ * this->transport_acc_[index_i] * inv_h_ratio * inv_h_ratio;
    }
}
//=================================================================================================//
template <typename... CommonControlTypes>
TransportVelocityCorrection<ContactBoundary<>, CommonControlTypes...>::
    TransportVelocityCorrection(BaseContactRelation &contact_relation)
    : TransportVelocityCorrection<Base, FluidContactData, CommonControlTypes...>(contact_relation) {}
//=================================================================================================//
template <typename... CommonControlTypes>
void TransportVelocityCorrection<ContactBoundary<>, CommonControlTypes...>::
    interaction(size_t index_i, Real dt)
{
    if (this->checkWithinScope(index_i))
    {
        Vecd acceleration_trans = Vecd::Zero();
        for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
        {
            Neighborhood &contact_neighborhood = (*this->contact_configuration_[k])[index_i];
            for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
            {
                // acceleration for transport velocity
                acceleration_trans -= 2.0 * this->kernel_correction_(index_i) *
                                      contact_neighborhood.dW_ijV_j_[n] * contact_neighborhood.e_ij_[n];
            }
        }
        this->transport_acc_[index_i] += acceleration_trans;
    }
}
//=================================================================================================//
template <class KernelCorrectionType, typename... CommonControlTypes>
TransportVelocityCorrection<Contact<>, KernelCorrectionType, CommonControlTypes...>::
    TransportVelocityCorrection(BaseContactRelation &contact_relation)
    : TransportVelocityCorrection<Base, FluidContactData, KernelCorrectionType, CommonControlTypes...>(
          contact_relation)
{
    for (size_t k = 0; k != this->contact_particles_.size(); ++k)
    {
        contact_kernel_corrections_.push_back(KernelCorrectionType(this->contact_particles_[k]));
    }
}
//=================================================================================================//
template <class KernelCorrectionType, typename... CommonControlTypes>
void TransportVelocityCorrection<Contact<>, KernelCorrectionType, CommonControlTypes...>::
    interaction(size_t index_i, Real dt)
{
    if (this->checkWithinScope(index_i))
    {
        Vecd acceleration_trans = Vecd::Zero();
        for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
        {
            Neighborhood &contact_neighborhood = (*this->contact_configuration_[k])[index_i];
            KernelCorrectionType &kernel_correction_k = this->contact_kernel_corrections_[k];
            for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
            {
                size_t index_j = contact_neighborhood.j_[n];
                // acceleration for transport velocity
                acceleration_trans -= (this->kernel_correction_(index_i) + kernel_correction_k(index_j)) *
                                      contact_neighborhood.dW_ijV_j_[n] * contact_neighborhood.e_ij_[n];
            }
        }
        this->transport_acc_[index_i] += acceleration_trans;
    }
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
