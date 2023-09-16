/**
 * @file 	fluid_dynamics_complex.hpp
 * @brief 	Here, we define the algorithm classes for complex fluid dynamics,
 * 			which is involving with either solid walls (with suffix WithWall)
 * 			or/and other bodies treated as wall for the fluid (with suffix Complex).
 * @author	Chi Zhang and Xiangyu Hu
 */
#pragma once

#include "fluid_dynamics_complex.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
template <template <class DataDelegationType> class BaseInteractionType>
InteractionWithWall<BaseInteractionType>::InteractionWithWall(BaseContactRelation &wall_contact_relation)
    : BaseInteractionType<FSIContactData>(wall_contact_relation)
{
    for (size_t k = 0; k != this->contact_particles_.size(); ++k)
    {
        Real rho0_k = this->contact_bodies_[k]->base_material_->ReferenceDensity();
        wall_inv_rho0_.push_back(1.0 / rho0_k);
        wall_mass_.push_back(&(this->contact_particles_[k]->mass_));
        wall_vel_ave_.push_back(this->contact_particles_[k]->AverageVelocity());
        wall_acc_ave_.push_back(this->contact_particles_[k]->AverageAcceleration());
        wall_n_.push_back(&(this->contact_particles_[k]->n_));
    }
}
//=================================================================================================//
template <class RiemannSolverType>
BaseIntegration1stHalfWithWall<RiemannSolverType>::
    BaseIntegration1stHalfWithWall(BaseContactRelation &wall_contact_relation)
    : InteractionWithWall<BaseIntegration>(wall_contact_relation),
      riemann_solver_(this->fluid_, this->fluid_) {}
//=================================================================================================//
template <class RiemannSolverType>
void BaseIntegration1stHalfWithWall<RiemannSolverType>::interaction(size_t index_i, Real dt)
{
    Vecd acceleration = Vecd::Zero();
    Real rho_dissipation(0);
    for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
    {
        StdLargeVec<Vecd> &acc_ave_k = *(this->wall_acc_ave_[k]);
        Neighborhood &wall_neighborhood = (*this->contact_configuration_[k])[index_i];
        for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
        {
            size_t index_j = wall_neighborhood.j_[n];
            Vecd &e_ij = wall_neighborhood.e_ij_[n];
            Real dW_ijV_j = wall_neighborhood.dW_ijV_j_[n];
            Real r_ij = wall_neighborhood.r_ij_[n];

            Real face_wall_external_acceleration = (this->acc_prior_[index_i] - acc_ave_k[index_j]).dot(-e_ij);
            Real p_in_wall = this->p_[index_i] + this->rho_[index_i] * r_ij * SMAX(Real(0), face_wall_external_acceleration);
            acceleration -= (this->p_[index_i] + p_in_wall) * dW_ijV_j * e_ij;
            rho_dissipation += this->riemann_solver_.DissipativeUJump(this->p_[index_i] - p_in_wall) * dW_ijV_j;
        }
    }
    this->acc_[index_i] += acceleration / this->rho_[index_i];
    this->drho_dt_[index_i] += rho_dissipation * this->rho_[index_i];
}
//=================================================================================================//
template <class RiemannSolverType>
BaseExtendIntegration1stHalfWithWall<RiemannSolverType>::
    BaseExtendIntegration1stHalfWithWall(BaseContactRelation &wall_contact_relation, Real penalty_strength)
    : BaseIntegration1stHalfWithWall<RiemannSolverType>(wall_contact_relation),
      penalty_strength_(penalty_strength) {}
//=================================================================================================//
template <class RiemannSolverType>
void BaseExtendIntegration1stHalfWithWall<RiemannSolverType>::interaction(size_t index_i, Real dt)
{
    BaseIntegration1stHalfWithWall<RiemannSolverType>::interaction(index_i, dt);

    Real rho_i = this->rho_[index_i];
    Real penalty_pressure = this->p_[index_i];
    Vecd acceleration = Vecd::Zero();
    for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
    {
        Real particle_spacing_j1 = 1.0 / this->contact_bodies_[k]->sph_adaptation_->ReferenceSpacing();
        Real particle_spacing_ratio2 = 1.0 / (this->sph_body_.sph_adaptation_->ReferenceSpacing() * particle_spacing_j1);
        particle_spacing_ratio2 *= 0.1 * particle_spacing_ratio2;

        StdLargeVec<Vecd> &n_k = *(this->wall_n_[k]);
        Neighborhood &wall_neighborhood = (*this->contact_configuration_[k])[index_i];
        for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
        {
            size_t index_j = wall_neighborhood.j_[n];
            Vecd &e_ij = wall_neighborhood.e_ij_[n];
            Real dW_ijV_j = wall_neighborhood.dW_ijV_j_[n];
            Real r_ij = wall_neighborhood.r_ij_[n];
            Vecd &n_j = n_k[index_j];

            /** penalty method to prevent particle running into boundary */
            Real projection = e_ij.dot(n_j);
            Real delta = 2.0 * projection * r_ij * particle_spacing_j1;
            Real beta = delta < 1.0 ? (1.0 - delta) * (1.0 - delta) * particle_spacing_ratio2 : 0.0;
            // penalty must be positive so that the penalty force is pointed to fluid inner domain
            Real penalty = penalty_strength_ * beta * fabs(projection * penalty_pressure);

            // penalty force induced acceleration
            acceleration -= 2.0 * penalty * n_j * dW_ijV_j / rho_i;
        }
    }
    this->acc_[index_i] += acceleration;
}
//=================================================================================================//
template <class RiemannSolverType>
BaseIntegration2ndHalfWithWall<RiemannSolverType>::
    BaseIntegration2ndHalfWithWall(BaseContactRelation &wall_contact_relation)
    : InteractionWithWall<BaseIntegration>(wall_contact_relation),
      riemann_solver_(this->fluid_, this->fluid_) {}
//=================================================================================================//
template <class RiemannSolverType>
void BaseIntegration2ndHalfWithWall<RiemannSolverType>::interaction(size_t index_i, Real dt)
{
    Real density_change_rate = 0.0;
    Vecd p_dissipation = Vecd::Zero();
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        StdLargeVec<Vecd> &vel_ave_k = *(wall_vel_ave_[k]);
        StdLargeVec<Vecd> &n_k = *(wall_n_[k]);
        Neighborhood &wall_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
        {
            size_t index_j = wall_neighborhood.j_[n];
            Vecd &e_ij = wall_neighborhood.e_ij_[n];
            Real dW_ijV_j = wall_neighborhood.dW_ijV_j_[n];

            Vecd vel_in_wall = 2.0 * vel_ave_k[index_j] - vel_[index_i];
            density_change_rate += (vel_[index_i] - vel_in_wall).dot(e_ij) * dW_ijV_j;
            Real u_jump = 2.0 * (vel_[index_i] - vel_ave_k[index_j]).dot(n_k[index_j]);
            p_dissipation += riemann_solver_.DissipativePJump(u_jump) * dW_ijV_j * n_k[index_j];
        }
    }
    drho_dt_[index_i] += density_change_rate * rho_[index_i];
    acc_[index_i] += p_dissipation / rho_[index_i];
}
//=================================================================================================//
} // namespace fluid_dynamics
  //=================================================================================================//
} // namespace SPH
  //=================================================================================================//