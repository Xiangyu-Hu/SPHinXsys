/* -------------------------------------------------------------------------*
 *								SPHinXsys									*
 * -------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle*
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
 * physical accurate simulation and aims to model coupled industrial dynamic*
 * systems including fluid, solid, multi-body dynamics and beyond with SPH	*
 * (smoothed particle hydrodynamics), a meshless computational method using	*
 * particle discretization.													*                                                                         *
 * ------------------------------------------------------------------------*/
#pragma once

#include "common_weakly_compressible_eulerian_classes.h"

namespace SPH
{
//=================================================================================================//
template <class BaseIntegrationType>
template <class BaseBodyRelationType>
InteractionWithWall<BaseIntegrationType>::
    InteractionWithWall(BaseBodyRelationType &base_body_relation, BaseContactRelation &wall_contact_relation)
    : BaseIntegrationType(base_body_relation), fluid_dynamics::FluidWallData(wall_contact_relation)
{
    if (&base_body_relation.getSPHBody() != &wall_contact_relation.getSPHBody())
    {
        std::cout << "\n Error: the two body_relations do not have the same source body!" << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
    }

    for (size_t k = 0; k != fluid_dynamics::FluidWallData::contact_particles_.size(); ++k)
    {
        Real rho_0_k = fluid_dynamics::FluidWallData::contact_bodies_[k]->base_material_->ReferenceDensity();
        wall_inv_rho0_.push_back(1.0 / rho_0_k);
        wall_vel_ave_.push_back(fluid_dynamics::FluidWallData::contact_particles_[k]->AverageVelocity());
        wall_acc_ave_.push_back(fluid_dynamics::FluidWallData::contact_particles_[k]->AverageAcceleration());
        wall_n_.push_back(&(fluid_dynamics::FluidWallData::contact_particles_[k]->n_));
    }
}
//=================================================================================================//
template <class BaseViscousAccelerationType>
template <class BaseBodyRelationType>
ViscousWithWall<BaseViscousAccelerationType>::
    ViscousWithWall(BaseBodyRelationType &base_body_relation, BaseContactRelation &wall_contact_relation)
    : InteractionWithWall<BaseViscousAccelerationType>(base_body_relation, wall_contact_relation){};
//=================================================================================================//
template <class BaseViscousAccelerationType>
void ViscousWithWall<BaseViscousAccelerationType>::interaction(size_t index_i, Real dt)
{
    BaseViscousAccelerationType::interaction(index_i, dt);

    Real rho_i = this->rho_[index_i];
    const Vecd &vel_i = this->vel_[index_i];

    Vecd acceleration = Vecd::Zero();
    Vecd vel_derivative = Vecd::Zero();
    for (size_t k = 0; k < fluid_dynamics::FluidWallData::contact_configuration_.size(); ++k)
    {
        StdLargeVec<Vecd> &vel_ave_k = *(this->wall_vel_ave_[k]);
        Neighborhood &contact_neighborhood = (*fluid_dynamics::FluidWallData::contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            Real r_ij = contact_neighborhood.r_ij_[n];

            vel_derivative = 2.0 * (vel_i - vel_ave_k[index_j]) / (r_ij + 0.01 * this->smoothing_length_);
            acceleration += 2.0 * this->mu_ * vel_derivative * contact_neighborhood.dW_ijV_j_[n] / rho_i;
        }
    }

    this->dmom_dt_prior_[index_i] += acceleration * rho_i;
}
//=================================================================================================//
template <class RiemannSolverType>
BaseIntegration1stHalf<RiemannSolverType>::BaseIntegration1stHalf(BaseInnerRelation &inner_relation, Real limiter_parameter)
    : EulerianBaseIntegration(inner_relation), limiter_input_(limiter_parameter), riemann_solver_(this->fluid_, this->fluid_, limiter_input_) {}
//=================================================================================================//
template <class RiemannSolverType>
void BaseIntegration1stHalf<RiemannSolverType>::interaction(size_t index_i, Real dt)
{
    FluidState state_i(rho_[index_i], vel_[index_i], p_[index_i]);
    Vecd momentum_change_rate = dmom_dt_prior_[index_i];
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Real dW_ijV_j = inner_neighborhood.dW_ijV_j_[n];
        Vecd &e_ij = inner_neighborhood.e_ij_[n];

        FluidState state_j(rho_[index_j], vel_[index_j], p_[index_j]);
        FluidStarState interface_state = riemann_solver_.getInterfaceState(state_i, state_j, e_ij);
        Real rho_star = this->fluid_.DensityFromPressure(interface_state.p_);

        momentum_change_rate -= 2.0 *
                                ((rho_star * interface_state.vel_) * interface_state.vel_.transpose() + interface_state.p_ * Matd::Identity()) * e_ij * dW_ijV_j;
    }
    dmom_dt_[index_i] = momentum_change_rate;
}
//=================================================================================================//
template <class RiemannSolverType>
void BaseIntegration1stHalf<RiemannSolverType>::update(size_t index_i, Real dt)
{
    mom_[index_i] += dmom_dt_[index_i] * dt;
    vel_[index_i] = mom_[index_i] / rho_[index_i];
}
//=================================================================================================//
template <class BaseIntegration1stHalfType>
void BaseIntegration1stHalfWithWall<BaseIntegration1stHalfType>::interaction(size_t index_i, Real dt)
{
    BaseIntegration1stHalfType::interaction(index_i, dt);

    FluidState state_i(this->rho_[index_i], this->vel_[index_i], this->p_[index_i]);

    Vecd momentum_change_rate = Vecd::Zero();
    for (size_t k = 0; k < fluid_dynamics::FluidWallData::contact_configuration_.size(); ++k)
    {
        StdLargeVec<Vecd> &n_k = *(this->wall_n_[k]);
        Neighborhood &wall_neighborhood = (*fluid_dynamics::FluidWallData::contact_configuration_[k])[index_i];
        for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
        {
            size_t index_j = wall_neighborhood.j_[n];
            Vecd &e_ij = wall_neighborhood.e_ij_[n];
            Real dW_ijV_j = wall_neighborhood.dW_ijV_j_[n];

            Vecd vel_in_wall = -state_i.vel_;
            Real p_in_wall = state_i.p_;
            Real rho_in_wall = state_i.rho_;
            FluidState state_j(rho_in_wall, vel_in_wall, p_in_wall);
            FluidStarState interface_state = this->riemann_solver_.getInterfaceState(state_i, state_j, n_k[index_j]);
            Real rho_star = this->fluid_.DensityFromPressure(interface_state.p_);

            momentum_change_rate -= 2.0 * ((rho_star * interface_state.vel_) * interface_state.vel_.transpose() + interface_state.p_ * Matd::Identity()) * e_ij * dW_ijV_j;
        }
    }
    this->dmom_dt_[index_i] += momentum_change_rate;
}
//=================================================================================================//
template <class RiemannSolverType>
BaseIntegration2ndHalf<RiemannSolverType>::BaseIntegration2ndHalf(BaseInnerRelation &inner_relation, Real limiter_parameter)
    : EulerianBaseIntegration(inner_relation), limiter_input_(limiter_parameter), riemann_solver_(this->fluid_, this->fluid_, limiter_input_) {}
//=================================================================================================//
template <class RiemannSolverType>
void BaseIntegration2ndHalf<RiemannSolverType>::interaction(size_t index_i, Real dt)
{
    FluidState state_i(rho_[index_i], vel_[index_i], p_[index_i]);
    Real density_change_rate = 0.0;
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Vecd &e_ij = inner_neighborhood.e_ij_[n];
        Real dW_ijV_j = inner_neighborhood.dW_ijV_j_[n];

        FluidState state_j(rho_[index_j], vel_[index_j], p_[index_j]);
        FluidStarState interface_state = riemann_solver_.getInterfaceState(state_i, state_j, e_ij);

        Real rho_star = this->fluid_.DensityFromPressure(interface_state.p_);
        density_change_rate -= 2.0 * (rho_star * interface_state.vel_).dot(e_ij) * dW_ijV_j;
    }
    drho_dt_[index_i] = density_change_rate;
}
//=================================================================================================//
template <class RiemannSolverType>
void BaseIntegration2ndHalf<RiemannSolverType>::update(size_t index_i, Real dt)
{
    rho_[index_i] += drho_dt_[index_i] * dt;
    p_[index_i] = fluid_.getPressure(rho_[index_i]);
}
//=================================================================================================//
template <class BaseIntegration2ndHalfType>
void BaseIntegration2ndHalfWithWall<BaseIntegration2ndHalfType>::interaction(size_t index_i, Real dt)
{
    BaseIntegration2ndHalfType::interaction(index_i, dt);

    FluidState state_i(this->rho_[index_i], this->vel_[index_i], this->p_[index_i]);
    Real density_change_rate = 0.0;
    for (size_t k = 0; k < fluid_dynamics::FluidWallData::contact_configuration_.size(); ++k)
    {
        StdLargeVec<Vecd> &n_k = *(this->wall_n_[k]);
        Neighborhood &wall_neighborhood = (*fluid_dynamics::FluidWallData::contact_configuration_[k])[index_i];
        for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
        {
            size_t index_j = wall_neighborhood.j_[n];
            Vecd &e_ij = wall_neighborhood.e_ij_[n];
            Real dW_ijV_j = wall_neighborhood.dW_ijV_j_[n];

            Vecd vel_in_wall = -state_i.vel_;
            Real p_in_wall = state_i.p_;
            Real rho_in_wall = state_i.rho_;

            FluidState state_j(rho_in_wall, vel_in_wall, p_in_wall);
            FluidStarState interface_state = this->riemann_solver_.getInterfaceState(state_i, state_j, n_k[index_j]);
            Real rho_star = this->fluid_.DensityFromPressure(interface_state.p_);

            density_change_rate -= 2.0 * (rho_star * interface_state.vel_).dot(e_ij) * dW_ijV_j;
        }
    }
    this->drho_dt_[index_i] += density_change_rate;
}
//=================================================================================================//
} // namespace SPH
  //=================================================================================================//