#pragma once

#include "transport_velocity_correction.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
template <class DataDelegationType, class KernelCorrectionType, class ResolutionType, class ParticleScope>
template <class BaseRelationType>
TransportVelocityCorrection<Base, DataDelegationType, KernelCorrectionType, ResolutionType, ParticleScope>::
    TransportVelocityCorrection(BaseRelationType &base_relation, Real coefficient)
    : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation),
      correction_scaling_(coefficient * pow(sph_body_.sph_adaptation_->ReferenceSmoothingLength(), 2)),
      pos_(this->particles_->pos_), kernel_correction_(this->particles_), h_ratio_(this->particles_),
      checkWithinScope(this->particles_) {}
//=================================================================================================//
template <typename... CommonControlTypes>
TransportVelocityCorrection<Inner<>, CommonControlTypes...>::
    TransportVelocityCorrection(BaseInnerRelation &inner_relation, Real coefficient)
    : TransportVelocityCorrection<Base, FluidDataInner, CommonControlTypes...>(
          inner_relation, coefficient) {}
//=================================================================================================//
template <typename... CommonControlTypes>
void TransportVelocityCorrection<Inner<>, CommonControlTypes...>::
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

        Real inv_h_ratio = 1.0 / this->h_ratio_(index_i);
        this->pos_[index_i] += this->correction_scaling_ * inv_h_ratio * inv_h_ratio * acceleration_trans;
    }
}
//=================================================================================================//
template <typename... CommonControlTypes>
TransportVelocityCorrection<ContactBoundary<>, CommonControlTypes...>::
    TransportVelocityCorrection(BaseContactRelation &contact_relation, Real coefficient)
    : TransportVelocityCorrection<Base, FluidContactData, CommonControlTypes...>(
          contact_relation, coefficient) {}
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
        Real inv_h_ratio = 1.0 / this->h_ratio_(index_i);
        this->pos_[index_i] += this->correction_scaling_ * inv_h_ratio * inv_h_ratio * acceleration_trans;
    }
}
//=================================================================================================//
template <class KernelCorrectionType, typename... CommonControlTypes>
TransportVelocityCorrection<Contact<>, KernelCorrectionType, CommonControlTypes...>::
    TransportVelocityCorrection(BaseContactRelation &contact_relation, Real coefficient)
    : TransportVelocityCorrection<Base, FluidContactData, KernelCorrectionType, CommonControlTypes...>(
          contact_relation, coefficient)
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
        Real inv_h_ratio = 1.0 / this->h_ratio_(index_i);
        this->pos_[index_i] += this->correction_scaling_ * inv_h_ratio * inv_h_ratio * acceleration_trans;
    }
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH