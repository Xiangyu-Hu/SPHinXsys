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
    : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation),
      rho_(*this->particles_->template getVariableByName<Real>("Density")),
      mass_(*this->particles_->template getVariableByName<Real>("Mass")),
      Vol_(*this->particles_->template getVariableByName<Real>("VolumetricMeasure")),
      vel_(*this->particles_->template getVariableByName<Vecd>("Velocity")),
      viscous_force_(*this->particles_->template registerSharedVariable<Vecd>("ViscousForce")),
      smoothing_length_(this->sph_body_.sph_adaptation_->ReferenceSmoothingLength()) {}
//=================================================================================================//
template <class ViscosityType>
ViscousForce<Inner<>, ViscosityType>::ViscousForce(BaseInnerRelation &inner_relation)
    : ViscousForce<FluidDataInner>(inner_relation),
      ForcePrior(&base_particles_, "ViscousForce"), mu_(&base_particles_)
{
    static_assert(std::is_base_of<ParticleAverage, ViscosityType>::value,
                  "ParticleAverage is not the base of ViscosityType!");
}
//=================================================================================================//
template <class ViscosityType>
void ViscousForce<Inner<>, ViscosityType>::interaction(size_t index_i, Real dt)
{
    Vecd force = Vecd::Zero();
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];

        // viscous force
        Vecd vel_derivative = (vel_[index_i] - vel_[index_j]) /
                              (inner_neighborhood.r_ij_[n] + 0.01 * smoothing_length_);
        force += 2.0 * mu_(index_i, index_j) * vel_derivative *
                 inner_neighborhood.dW_ij_[n] * Vol_[index_j];
    }

    viscous_force_[index_i] = force * Vol_[index_i];
}
//=================================================================================================//
template <typename ViscosityType>
void ViscousForce<Inner<AngularConservative>, ViscosityType>::interaction(size_t index_i, Real dt)
{
    Vecd force = Vecd::Zero();
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Vecd &e_ij = inner_neighborhood.e_ij_[n];
        Real r_ij = inner_neighborhood.r_ij_[n];

        /** The following viscous force is given in Monaghan 2005 (Rep. Prog. Phys.), it seems that
         * this formulation is more accurate than the previous one for Taylor-Green-Vortex flow. */
        Real v_r_ij = (vel_[index_i] - vel_[index_j]).dot(e_ij);
        Real eta_ij = 2.0 * Real(Dimensions + 2) * mu_(index_i, index_j) * v_r_ij /
                      (r_ij + 0.01 * smoothing_length_);
        force += eta_ij * inner_neighborhood.dW_ij_[n] * Vol_[index_j] * e_ij;
    }

    viscous_force_[index_i] = force * Vol_[index_i];
}
//=================================================================================================//
template <typename ViscosityType>
void ViscousForce<Contact<Wall>, ViscosityType>::interaction(size_t index_i, Real dt)
{
    Vecd force = Vecd::Zero();
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        StdLargeVec<Vecd> &vel_ave_k = *(wall_vel_ave_[k]);
        StdLargeVec<Real> &wall_Vol_k = *(wall_Vol_[k]);
        Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            Real r_ij = contact_neighborhood.r_ij_[n];

            Vecd vel_derivative = 2.0 * (vel_[index_i] - vel_ave_k[index_j]) /
                                  (r_ij + 0.01 * smoothing_length_);
            force += 2.0 * mu_(index_i, index_i) * vel_derivative *
                     contact_neighborhood.dW_ij_[n] * wall_Vol_k[index_j];
        }
    }

    viscous_force_[index_i] += force * Vol_[index_i];
}
//=================================================================================================//
template <typename ViscosityType>
ViscousForce<Contact<Wall, AngularConservative>, ViscosityType>::
    ViscousForce(BaseContactRelation &wall_contact_relation)
    : BaseViscousForceWithWall(wall_contact_relation),
      mu_(&base_particles_),
      distance_from_wall_(*this->particles_->template registerSharedVariable<Vecd>("DistanceFromWall")){};
//=================================================================================================//
template <typename ViscosityType>
void ViscousForce<Contact<Wall, AngularConservative>, ViscosityType>::interaction(size_t index_i, Real dt)
{
    Vecd force = Vecd::Zero();
    const Vecd &distance_from_wall = distance_from_wall_[index_i];
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        StdLargeVec<Vecd> &vel_ave_k = *(wall_vel_ave_[k]);
        StdLargeVec<Real> &wall_Vol_k = *(wall_Vol_[k]);
        Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            Vecd &e_ij = contact_neighborhood.e_ij_[n];
            Real r_ij = contact_neighborhood.r_ij_[n];

            Vecd distance_diff = distance_from_wall - r_ij * e_ij;
            Real factor = 1.0 - distance_from_wall.dot(distance_diff) / distance_from_wall.squaredNorm();
            Real v_r_ij = factor * (vel_[index_i] - vel_ave_k[index_j]).dot(e_ij);
            Real eta_ij = 2.0 * Real(Dimensions + 2) * mu_(index_i, index_i) * v_r_ij /
                          (r_ij + 0.01 * smoothing_length_);
            force += eta_ij * contact_neighborhood.dW_ij_[n] * wall_Vol_k[index_j] * e_ij;
        }
    }

    viscous_force_[index_i] += force * Vol_[index_i];
}
//=================================================================================================//
template <typename ViscosityType>
ViscousForce<Contact<>, ViscosityType>::ViscousForce(BaseContactRelation &contact_relation)
    : ViscousForce<FluidContactData>(contact_relation)
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        contact_mu_.emplace_back(ViscosityType(particles_, contact_particles_[k]));
        contact_vel_.push_back(contact_particles_[k]->template getVariableByName<Vecd>("Velocity"));
        wall_Vol_.push_back(contact_particles_[k]->template getVariableByName<Real>("VolumetricMeasure"));
    }
}
//=================================================================================================//
template <typename ViscosityType>
void ViscousForce<Contact<>, ViscosityType>::interaction(size_t index_i, Real dt)
{
    Vecd force = Vecd::Zero();
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        auto &contact_mu_k = contact_mu_[k];
        StdLargeVec<Vecd> &vel_k = *(contact_vel_[k]);
        StdLargeVec<Real> &wall_Vol_k = *(wall_Vol_[k]);
        Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            Vecd vel_derivative = (vel_[index_i] - vel_k[index_j]) /
                                  (contact_neighborhood.r_ij_[n] + 0.01 * smoothing_length_);
            force += 2.0 * contact_mu_k(index_i, index_j) * vel_derivative *
                     contact_neighborhood.dW_ij_[n] * wall_Vol_k[index_j];
        }
    }
    viscous_force_[index_i] += force * Vol_[index_i];
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
