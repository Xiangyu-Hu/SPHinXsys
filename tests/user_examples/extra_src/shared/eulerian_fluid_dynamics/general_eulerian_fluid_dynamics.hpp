#pragma once

#include "general_eulerian_fluid_dynamics.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
template <class RiemannSolverType>
ICEIntegration1stHalf<RiemannSolverType>::ICEIntegration1stHalf(BaseInnerRelation &inner_relation)
    : BaseIntegration(inner_relation), riemann_solver_(this->fluid_, this->fluid_),
      acc_prior_(particles_->acc_prior_)
{
    particles_->registerVariable(mom_, "Momentum");
    particles_->registerVariable(dmom_dt_, "MomentumChangeRate");
};
//=================================================================================================//
template <class RiemannSolverType>
void ICEIntegration1stHalf<RiemannSolverType>::interaction(size_t index_i, Real dt)
{
    FluidState state_i(rho_[index_i], vel_[index_i], p_[index_i]);
    Vecd momentum_change_rate = Vecd::Zero();
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Real dW_ijV_j = inner_neighborhood.dW_ijV_j_[n];
        Vecd &e_ij = inner_neighborhood.e_ij_[n];
        FluidState state_j(rho_[index_j], vel_[index_j], p_[index_j]);

        Vec2d wave_speeds = riemann_solver_.getBoundingWaveSpeeds(state_i, state_j, e_ij);
        Real s_l = wave_speeds[0];
        Real s_r = wave_speeds[1];

        if (s_l < 0 && s_r > 0)
        {
            Matd flux_l = rho_[index_i] * vel_[index_i] * vel_[index_i].transpose() + p_[index_i] * Matd::Identity();
            Matd flux_r = rho_[index_j] * vel_[index_j] * vel_[index_j].transpose() + p_[index_j] * Matd::Identity();
            DissipationState dissipation_state = riemann_solver_.getDissipationState(state_i, state_j, e_ij);
            Matd dissipation_term = s_l * s_r * dissipation_state.momentum_dissipation_ / (s_r - s_l);
            momentum_change_rate -= 2 * ((s_r * flux_l - s_l * flux_r) / (s_r - s_l) + dissipation_term) * e_ij * dW_ijV_j;
        }
        else if(s_l > 0)
        {
            std::cout << "s_l is greater than 0" << std::endl;
        }
        else if (s_r < 0)
        {
            std::cout << "s_r is smaller than 0" << std::endl;
        }
    }
    dmom_dt_[index_i] = momentum_change_rate;
};
//=================================================================================================//
template <class RiemannSolverType>
void ICEIntegration1stHalf<RiemannSolverType>::update(size_t index_i, Real dt)
{
    mom_[index_i] += (dmom_dt_[index_i] + rho_[index_i] * acc_prior_[index_i]) * dt;
    vel_[index_i] = mom_[index_i] / rho_[index_i];
};
//=================================================================================================//
template <class ICEIntegration1stHalfType>
void ICEIntegration1stHalfWithWall<ICEIntegration1stHalfType>::interaction(size_t index_i, Real dt)
{
    ICEIntegration1stHalfType::interaction(index_i, dt);

    FluidState state_i(this->rho_[index_i], this->vel_[index_i], this->p_[index_i]);
    Vecd momentum_change_rate = Vecd::Zero();
    for (size_t k = 0; k < FluidWallData::contact_configuration_.size(); ++k)
    {
        StdLargeVec<Vecd> &n_k = *(this->wall_n_[k]);
        Neighborhood &wall_neighborhood = (*FluidWallData::contact_configuration_[k])[index_i];
        for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
        {
            size_t index_j = wall_neighborhood.j_[n];
            Vecd &e_ij = wall_neighborhood.e_ij_[n];
            Real dW_ijV_j = wall_neighborhood.dW_ijV_j_[n];

            Vecd vel_in_wall = -this->vel_[index_i];
            Real p_in_wall = this->p_[index_i];
            Real rho_in_wall = this->rho_[index_i];
            FluidState state_j(rho_in_wall, vel_in_wall, p_in_wall);

            Vec2d wave_speeds = this->riemann_solver_.getBoundingWaveSpeeds(state_i, state_j, e_ij);
            Real s_l = wave_speeds[0];
            Real s_r = wave_speeds[1];
            
            if (s_l < 0 && s_r > 0)
            {
                Matd flux_l = this->rho_[index_i] * this->vel_[index_i] * this->vel_[index_i].transpose() + this->p_[index_i] * Matd::Identity();
                Matd flux_r = rho_in_wall * vel_in_wall * vel_in_wall.transpose() + p_in_wall * Matd::Identity();
                DissipationState dissipation_state = this->riemann_solver_.getDissipationState(state_i, state_j, e_ij);
                Matd dissipation_term = s_l * s_r * dissipation_state.momentum_dissipation_ / (s_r - s_l);
                momentum_change_rate -= 2 * ((s_r * flux_l - s_l * flux_r) / (s_r - s_l) + dissipation_term) * e_ij * dW_ijV_j;
            }
        }
    };
    this->dmom_dt_[index_i] += momentum_change_rate;
}
//=================================================================================================//
template <class RiemannSolverType>
ICEIntegration2ndHalf<RiemannSolverType>::
    ICEIntegration2ndHalf(BaseInnerRelation &inner_relation)
    : BaseIntegration(inner_relation), riemann_solver_(this->fluid_, this->fluid_) {}
//=================================================================================================//
template<class RiemannSolverType>
void ICEIntegration2ndHalf<RiemannSolverType>::interaction(size_t index_i, Real dt)
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

        Vec2d wave_speeds = riemann_solver_.getBoundingWaveSpeeds(state_i, state_j, e_ij);
        Real s_l = wave_speeds[0];
        Real s_r = wave_speeds[1];

        if (s_l < 0 && s_r > 0)
        {
            Vecd flux_l = rho_[index_i] * vel_[index_i];
            Vecd flux_r = rho_[index_j] * vel_[index_j];
            DissipationState dissipation_state = riemann_solver_.getDissipationState(state_i, state_j, e_ij);
            Vecd dissipation_term = s_l * s_r * dissipation_state.density_dissipation_ / (s_r - s_l);
            density_change_rate -= 2 * ((s_r * flux_l - s_l * flux_r) / (s_r - s_l) + dissipation_term).dot(e_ij) * dW_ijV_j;
        }
    }
    drho_dt_[index_i] = density_change_rate;
};
//=================================================================================================//
template <class RiemannSolverType>
void ICEIntegration2ndHalf<RiemannSolverType>::update(size_t index_i, Real dt)
{
    rho_[index_i] += drho_dt_[index_i] * dt;
    p_[index_i] = fluid_.getPressure(rho_[index_i]);
}
//=================================================================================================//
template <class ICEIntegration2ndHalfType>
void ICEIntegration2ndHalfWithWall<ICEIntegration2ndHalfType>::interaction(size_t index_i, Real dt)
{
    ICEIntegration2ndHalfType::interaction(index_i, dt);

    FluidState state_i(this->rho_[index_i], this->vel_[index_i], this->p_[index_i]);
    Real density_change_rate = 0.0;
    for (size_t k = 0; k < FluidWallData::contact_configuration_.size(); ++k)
    {
        StdLargeVec<Vecd> &n_k = *(this->wall_n_[k]);
        Neighborhood &wall_neighborhood = (*FluidWallData::contact_configuration_[k])[index_i];
        for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
        {
            size_t index_j = wall_neighborhood.j_[n];
            Vecd &e_ij = wall_neighborhood.e_ij_[n];
            Real dW_ijV_j = wall_neighborhood.dW_ijV_j_[n];

            Vecd vel_in_wall = -state_i.vel_;
            Real p_in_wall = state_i.p_;
            Real rho_in_wall = state_i.rho_;
            FluidState state_j(rho_in_wall, vel_in_wall, p_in_wall);

           Vec2d wave_speeds = this->riemann_solver_.getBoundingWaveSpeeds(state_i, state_j, e_ij);
           Real s_l = wave_speeds[0];
           Real s_r = wave_speeds[1];

            if (s_l < 0 && s_r > 0)
            {
                Vecd flux_l = this->rho_[index_i] * this->vel_[index_i];
                Vecd flux_r = rho_in_wall * vel_in_wall;
                DissipationState dissipation_state = this->riemann_solver_.getDissipationState(state_i, state_j, e_ij);
                Vecd dissipation_term = s_l * s_r * dissipation_state.density_dissipation_ / (s_r - s_l);
                density_change_rate -= 2 * ((s_r * flux_l - s_l * flux_r) / (s_r - s_l) + dissipation_term).dot(e_ij) * dW_ijV_j;
            }
        }
    }
    this->drho_dt_[index_i] += density_change_rate;
};
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
