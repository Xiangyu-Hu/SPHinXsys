#pragma once

#include "viscous_dynamics.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
template <class DataDelegationType>
template <class BaseRelationType>
ViscousForce<DataDelegationType>::ViscousForce(BaseRelationType &base_relation)
    : ForcePrior(base_relation.getSPHBody(), "ViscousForce"),
      DataDelegationType(base_relation),
      rho_(this->particles_->template getVariableDataByName<Real>("Density")),
      mass_(this->particles_->template getVariableDataByName<Real>("Mass")),
      Vol_(this->particles_->template getVariableDataByName<Real>("VolumetricMeasure")),
      vel_(this->particles_->template getVariableDataByName<Vecd>("Velocity")),
      viscous_force_(this->particles_->template registerStateVariable<Vecd>("ViscousForce")),
      smoothing_length_(this->sph_body_.getSPHAdaptation().ReferenceSmoothingLength()) {}
//=================================================================================================//
template <typename ViscosityType, class KernelCorrectionType>
ViscousForce<Inner<>, ViscosityType, KernelCorrectionType>::
    ViscousForce(BaseInnerRelation &inner_relation)
    : ViscousForce<DataDelegateInner>(inner_relation),
      mu_(particles_), kernel_correction_(particles_)
{
    static_assert(std::is_base_of<ParticleAverage, ViscosityType>::value,
                  "ParticleAverage is not the base of ViscosityType!");
}
//=================================================================================================//
template <typename ViscosityType, class KernelCorrectionType>
void ViscousForce<Inner<>, ViscosityType, KernelCorrectionType>::interaction(size_t index_i, Real dt)
{
    Vecd force = Vecd::Zero();
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        const Vecd &e_ij = inner_neighborhood.e_ij_[n];

        // viscous force
        Vecd vel_derivative = (vel_[index_i] - vel_[index_j]) /
                              (inner_neighborhood.r_ij_[n] + 0.01 * smoothing_length_);
        force += e_ij.dot((kernel_correction_(index_i) + kernel_correction_(index_j)) * e_ij) *
                 mu_(index_i, index_j) * vel_derivative *
                 inner_neighborhood.dW_ij_[n] * Vol_[index_j];
    }

    viscous_force_[index_i] = force * Vol_[index_i];
}
//=================================================================================================//
template <typename ViscosityType, class KernelCorrectionType>
ViscousForce<Inner<AngularConservative>, ViscosityType, KernelCorrectionType>::
    ViscousForce(BaseInnerRelation &inner_relation)
    : ViscousForce<DataDelegateInner>(inner_relation),
      mu_(particles_), kernel_correction_(particles_) {}
//=================================================================================================//
template <typename ViscosityType, class KernelCorrectionType>
void ViscousForce<Inner<AngularConservative>, ViscosityType, KernelCorrectionType>::
    interaction(size_t index_i, Real dt)
{
    Vecd force = Vecd::Zero();
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        const Vecd &e_ij = inner_neighborhood.e_ij_[n];
        Real r_ij = inner_neighborhood.r_ij_[n];

        /** The following viscous force is given in Monaghan 2005 (Rep. Prog. Phys.), it seems that
         * this formulation is more accurate than the previous one for Taylor-Green-Vortex flow. */
        Real v_r_ij = (vel_[index_i] - vel_[index_j]).dot(e_ij);
        Real eta_ij = Real(Dimensions + 2) * mu_(index_i, index_j) * v_r_ij /
                      (r_ij + 0.01 * smoothing_length_);
        force += e_ij.dot((kernel_correction_(index_i) + kernel_correction_(index_j)) * e_ij) *
                 eta_ij * inner_neighborhood.dW_ij_[n] * Vol_[index_j] * e_ij;
    }

    viscous_force_[index_i] = force * Vol_[index_i];
}
//=================================================================================================//
template <typename ViscosityType, class KernelCorrectionType>
ViscousForce<Contact<Wall>, ViscosityType, KernelCorrectionType>::
    ViscousForce(BaseContactRelation &wall_contact_relation)
    : BaseViscousForceWithWall(wall_contact_relation),
      mu_(particles_), kernel_correction_(particles_) {}
//=================================================================================================//
template <typename ViscosityType, class KernelCorrectionType>
void ViscousForce<Contact<Wall>, ViscosityType, KernelCorrectionType>::
    interaction(size_t index_i, Real dt)
{
    Vecd force = Vecd::Zero();
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        Vecd *vel_ave_k = wall_vel_ave_[k];
        Real *wall_Vol_k = wall_Vol_[k];
        const Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            Real r_ij = contact_neighborhood.r_ij_[n];
            const Vecd &e_ij = contact_neighborhood.e_ij_[n];

            Vecd vel_derivative = 2.0 * (vel_[index_i] - vel_ave_k[index_j]) /
                                  (r_ij + 0.01 * smoothing_length_);
            force += 2.0 * e_ij.dot(kernel_correction_(index_i) * e_ij) * mu_(index_i, index_i) *
                     vel_derivative * contact_neighborhood.dW_ij_[n] * wall_Vol_k[index_j];
        }
    }

    viscous_force_[index_i] += force * Vol_[index_i];
}
//=================================================================================================//
template <typename ViscosityType, class KernelCorrectionType>
ViscousForce<Contact<Wall, AngularConservative>, ViscosityType, KernelCorrectionType>::
    ViscousForce(BaseContactRelation &wall_contact_relation)
    : BaseViscousForceWithWall(wall_contact_relation),
      mu_(particles_), kernel_correction_(particles_){};
//=================================================================================================//
template <typename ViscosityType, class KernelCorrectionType>
void ViscousForce<Contact<Wall, AngularConservative>, ViscosityType, KernelCorrectionType>::
    interaction(size_t index_i, Real dt)
{
    Vecd force = Vecd::Zero();
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        Vecd *vel_ave_k = wall_vel_ave_[k];
        Real *wall_Vol_k = wall_Vol_[k];
        const Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            const Vecd &e_ij = contact_neighborhood.e_ij_[n];
            Real r_ij = contact_neighborhood.r_ij_[n];

            Real v_r_ij = 2.0 * (vel_[index_i] - vel_ave_k[index_j]).dot(e_ij);
            Real eta_ij = Real(Dimensions + 2) * mu_(index_i, index_i) * v_r_ij /
                          (r_ij + 0.01 * smoothing_length_);
            force += 2.0 * e_ij.dot(kernel_correction_(index_i) * e_ij) * mu_(index_i, index_i) *
                     eta_ij * contact_neighborhood.dW_ij_[n] * wall_Vol_k[index_j] * e_ij;
        }
    }

    viscous_force_[index_i] += force * Vol_[index_i];
}
//=================================================================================================//
template <typename ViscosityType, class KernelCorrectionType>
ViscousForce<Contact<>, ViscosityType, KernelCorrectionType>::
    ViscousForce(BaseContactRelation &contact_relation)
    : ViscousForce<DataDelegateContact>(contact_relation),
      kernel_correction_(particles_)
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        contact_mu_.emplace_back(ViscosityType(particles_, contact_particles_[k]));
        contact_kernel_corrections_.emplace_back(KernelCorrectionType(contact_particles_[k]));
        contact_vel_.push_back(contact_particles_[k]->template getVariableDataByName<Vecd>("Velocity"));
        wall_Vol_.push_back(contact_particles_[k]->template getVariableDataByName<Real>("VolumetricMeasure"));
    }
}
//=================================================================================================//
template <typename ViscosityType, class KernelCorrectionType>
void ViscousForce<Contact<>, ViscosityType, KernelCorrectionType>::
    interaction(size_t index_i, Real dt)
{
    Vecd force = Vecd::Zero();
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        auto &contact_mu_k = contact_mu_[k];
        Vecd *vel_k = contact_vel_[k];
        Real *wall_Vol_k = wall_Vol_[k];
        KernelCorrectionType &kernel_correction_k = contact_kernel_corrections_[k];
        const Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            const Vecd &e_ij = contact_neighborhood.e_ij_[n];

            Vecd vel_derivative = (vel_[index_i] - vel_k[index_j]) /
                                  (contact_neighborhood.r_ij_[n] + 0.01 * smoothing_length_);
            force += e_ij.dot((kernel_correction_(index_i) + kernel_correction_k(index_j)) * e_ij) *
                     contact_mu_k(index_i, index_j) * vel_derivative *
                     contact_neighborhood.dW_ij_[n] * wall_Vol_k[index_j];
        }
    }
    viscous_force_[index_i] += force * Vol_[index_i];
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
