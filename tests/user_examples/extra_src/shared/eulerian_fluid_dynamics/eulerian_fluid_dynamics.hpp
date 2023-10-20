#pragma once

#include "eulerian_fluid_dynamics.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
template <class RiemannSolverType>
EulerianIntegration1stHalf<RiemannSolverType>::
    EulerianIntegration1stHalf(BaseInnerRelation &inner_relation)
    : BaseIntegration(inner_relation), riemann_solver_(this->fluid_, this->fluid_),
      acc_prior_(particles_->acc_prior_)
{
    particles_->registerVariable(mom_, "Momentum");
    particles_->registerVariable(dmom_dt_, "MomentumChangeRate");
}
//=================================================================================================//
template <class RiemannSolverType>
void EulerianIntegration1stHalf<RiemannSolverType>::interaction(size_t index_i, Real dt)
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
        FluidStarState interface_state = riemann_solver_.getInterfaceState(state_i, state_j, e_ij);
        Real rho_star = this->fluid_.DensityFromPressure(interface_state.p_);

        momentum_change_rate -= 2.0 * ((rho_star * interface_state.vel_) * interface_state.vel_.transpose() + interface_state.p_ * Matd::Identity()) * e_ij * dW_ijV_j;
    }
    dmom_dt_[index_i] = momentum_change_rate;
}
//=================================================================================================//
template <class RiemannSolverType>
void EulerianIntegration1stHalf<RiemannSolverType>::update(size_t index_i, Real dt)
{
    mom_[index_i] += (dmom_dt_[index_i] + rho_[index_i] * acc_prior_[index_i]) * dt;
    vel_[index_i] = mom_[index_i] / rho_[index_i];
}
//=================================================================================================//
template <class EulerianIntegration1stHalfType>
void EulerianIntegration1stHalfWithWall<EulerianIntegration1stHalfType>::interaction(size_t index_i, Real dt)
{
    EulerianIntegration1stHalfType::interaction(index_i, dt);

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
EulerianIntegration2ndHalf<RiemannSolverType>::
    EulerianIntegration2ndHalf(BaseInnerRelation &inner_relation)
    : BaseIntegration(inner_relation), riemann_solver_(this->fluid_, this->fluid_) {}
//=================================================================================================//
template <class RiemannSolverType>
void EulerianIntegration2ndHalf<RiemannSolverType>::interaction(size_t index_i, Real dt)
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
void EulerianIntegration2ndHalf<RiemannSolverType>::update(size_t index_i, Real dt)
{
    rho_[index_i] += drho_dt_[index_i] * dt;
    p_[index_i] = fluid_.getPressure(rho_[index_i]);
}
//=================================================================================================//
template <class EulerianIntegration2ndHalfType>
void EulerianIntegration2ndHalfWithWall<EulerianIntegration2ndHalfType>::interaction(size_t index_i, Real dt)
{
    EulerianIntegration2ndHalfType::interaction(index_i, dt);

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
            FluidStarState interface_state = this->riemann_solver_.getInterfaceState(state_i, state_j, n_k[index_j]);
            Real rho_star = this->fluid_.DensityFromPressure(interface_state.p_);

            density_change_rate -= 2.0 * (rho_star * interface_state.vel_).dot(e_ij) * dW_ijV_j;
        }
    }
    this->drho_dt_[index_i] += density_change_rate;
}
//=================================================================================================//
template <class RiemannSolverType>
EulerianIntegration1stHalfConsistency<RiemannSolverType>::
    EulerianIntegration1stHalfConsistency(BaseInnerRelation &inner_relation)
    : BaseIntegration(inner_relation), riemann_solver_(this->fluid_, this->fluid_),
      acc_prior_(particles_->acc_prior_)
{
    particles_->registerVariable(mom_, "Momentum");
    particles_->registerVariable(dmom_dt_, "MomentumChangeRate");
}
//=================================================================================================//
template <class RiemannSolverType>
void EulerianIntegration1stHalfConsistency<RiemannSolverType>::interaction(size_t index_i, Real dt)
{
    Vecd momentum_change_rate = Vecd::Zero();
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    Real p_i = this->p_[index_i];
    Real rho_i = this->rho_[index_i];
    Vecd vel_i = this->vel_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Real dW_ijV_j = inner_neighborhood.dW_ijV_j_[n];
        Vecd &e_ij = inner_neighborhood.e_ij_[n];

        Real p_j = this->p_[index_j];
        Real rho_j = this->rho_[index_j];
        Vecd vel_j = this->vel_[index_j];

        Real p_jump = (p_i - p_j);
        Real u_jump = (vel_i - vel_j).dot(-e_ij);

        Real rho_ave = 2.0 * rho_i * rho_j / (rho_i + rho_j);
        
        Matd rho_vel_sqr_ave = rho_ave * (rho_i * vel_i * vel_i.transpose() + rho_j * vel_j * vel_j.transpose()) / (rho_i + rho_j);
        Matd p_ave = (rho_i * p_i + rho_j * p_j) / (rho_i + rho_j) * Matd::Identity();

        //Matd rho_vel_sqr_ave = rho_ave * (rho_i * vel_i * vel_i.transpose() * B_[index_j] + rho_j * vel_j * vel_j.transpose() * B_[index_i]) / (rho_i + rho_j);
        //Matd p_ave = (rho_i * p_i * B_[index_j] + rho_j * p_j * B_[index_i]) / (rho_i + rho_j) * Matd::Identity();

        Matd rho_dissipation = rho_ave * (2 * (rho_i * vel_i + rho_j * vel_j) / (rho_i + rho_j) * riemann_solver_.DissipativeUJump(p_jump, u_jump) * (-e_ij).transpose() + 
                                         riemann_solver_.DissipativeUJump(p_jump, u_jump) * riemann_solver_.DissipativeUJump(p_jump, u_jump) * e_ij * e_ij.transpose());
        Matd p_dissipation = riemann_solver_.DissipativePJump(u_jump) * Matd::Identity();

        momentum_change_rate -= 2.0 * (rho_vel_sqr_ave + p_ave + rho_dissipation + p_dissipation) * e_ij * dW_ijV_j;
    }
    dmom_dt_[index_i] = momentum_change_rate;
};
//=================================================================================================//
template <class RiemannSolverType>
void EulerianIntegration1stHalfConsistency<RiemannSolverType>::update(size_t index_i, Real dt)
{
    mom_[index_i] += (dmom_dt_[index_i] + rho_[index_i] * acc_prior_[index_i]) * dt;
    vel_[index_i] = mom_[index_i] / rho_[index_i];
};
//=================================================================================================//
template <class EulerianIntegration1stHalfType>
void EulerianIntegration1stHalfWithWallConsistency<EulerianIntegration1stHalfType>::interaction(size_t index_i, Real dt)
{
    EulerianIntegration1stHalfType::interaction(index_i, dt);

    Vecd momentum_change_rate = Vecd::Zero();
    Real p_i = this->p_[index_i];
    Real rho_i = this->rho_[index_i];
    Vecd vel_i = this->vel_[index_i];
    for (size_t k = 0; k < FluidWallData::contact_configuration_.size(); ++k)
    {
        StdLargeVec<Vecd> &n_k = *(this->wall_n_[k]);
        StdLargeVec<Matd> &B_k = *(this->wall_B_[k]);
        Neighborhood &wall_neighborhood = (*FluidWallData::contact_configuration_[k])[index_i];
        for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
        {
            size_t index_j = wall_neighborhood.j_[n];
            Vecd &e_ij = wall_neighborhood.e_ij_[n];
            Real dW_ijV_j = wall_neighborhood.dW_ijV_j_[n];

            Vecd vel_in_wall = -vel_i;
            Real p_in_wall = p_i;
            Real rho_in_wall = rho_i;

            Real p_jump = (p_i - p_in_wall);
            Real u_jump = (vel_i - vel_in_wall).dot(-e_ij);

            Real rho_ave = 2 * rho_i * rho_in_wall / (rho_i + rho_in_wall);

            Matd rho_vel_sqr_ave = rho_ave * (rho_i * vel_i * vel_i.transpose() + rho_in_wall * vel_in_wall * vel_in_wall.transpose()) / (rho_i + rho_in_wall);
            Matd p_ave = (rho_i * p_i + rho_in_wall * p_in_wall) / (rho_i + rho_in_wall) * Matd::Identity();

            //Matd rho_vel_sqr_ave = rho_ave * (rho_i * vel_i * vel_i.transpose() * B_k[index_j] + rho_in_wall * vel_in_wall * vel_in_wall.transpose() * this->B_[index_i]) / (rho_i + rho_in_wall);
            //Matd p_ave = (rho_i * p_i * B_k[index_j] + rho_in_wall * p_in_wall * this->B_[index_i]) / (rho_i + rho_in_wall) * Matd::Identity();

            Matd rho_dissipation = rho_ave * (2 * (rho_i * vel_i + rho_in_wall * vel_in_wall) / (rho_i + rho_in_wall) * (this->riemann_solver_.DissipativeUJump(p_jump, u_jump)) * (-e_ij).transpose() +
                                              (this->riemann_solver_.DissipativeUJump(p_jump, u_jump)) * (this->riemann_solver_.DissipativeUJump(p_jump, u_jump)) * e_ij * e_ij.transpose());
            Matd p_dissipation = this->riemann_solver_.DissipativePJump(u_jump) * Matd::Identity();

            momentum_change_rate -= 2 * (rho_vel_sqr_ave + p_ave + rho_dissipation + p_dissipation) * e_ij * dW_ijV_j;
        }
    }
    this->dmom_dt_[index_i] += momentum_change_rate;
};
//=================================================================================================//
template <class RiemannSolverType>
EulerianIntegration2ndHalfConsistency<RiemannSolverType>::
    EulerianIntegration2ndHalfConsistency(BaseInnerRelation &inner_relation)
    : BaseIntegration(inner_relation), riemann_solver_(this->fluid_, this->fluid_){};
//=================================================================================================//
template <class RiemannSolverType>
void EulerianIntegration2ndHalfConsistency<RiemannSolverType>::interaction(size_t index_i, Real dt)
{
    Real density_change_rate = 0.0;
    Neighborhood& inner_neighborhood = inner_configuration_[index_i];
    Real p_i = this->p_[index_i];
    Real rho_i = this->rho_[index_i];
    Vecd vel_i = this->vel_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Vecd& e_ij = inner_neighborhood.e_ij_[n];
        Real dW_ijV_j = inner_neighborhood.dW_ijV_j_[n];

        Real p_j = this->p_[index_j];
        Real rho_j = this->rho_[index_j];
        Vecd vel_j = this->vel_[index_j];

        Real p_jump = p_i - p_j;
        Real u_jump = (vel_i - vel_j).dot(-e_ij);

        Real rho_ave = 2 * rho_i * rho_j / (rho_i + rho_j);

        Vecd rho_vel_ave = 0.5 * (rho_i * vel_i + rho_j * vel_j);
        //Vecd rho_vel_ave = 0.5 * (rho_i * this->B_[index_j] * vel_i + rho_j * this->B_[index_i] * vel_j);

        Vecd rho_dissipation = rho_ave * riemann_solver_.DissipativeUJump(p_jump, u_jump) * (-e_ij);
        density_change_rate -= 2 * (rho_vel_ave + rho_dissipation).dot(e_ij) * dW_ijV_j;
    }
    drho_dt_[index_i] = density_change_rate;
};
//=================================================================================================//
template <class RiemannSolverType>
void EulerianIntegration2ndHalfConsistency<RiemannSolverType>::update(size_t index_i, Real dt)
{
    rho_[index_i] += drho_dt_[index_i] * dt;
    p_[index_i] = fluid_.getPressure(rho_[index_i]);
};
//=================================================================================================//
template <class EulerianIntegration2ndHalfType>
void EulerianIntegration2ndHalfWithWallConsistency<EulerianIntegration2ndHalfType>::interaction(size_t index_i, Real dt)
{
    EulerianIntegration2ndHalfType::interaction(index_i, dt);

    Real density_change_rate = 0.0;
    Real p_i = this->p_[index_i];
    Real rho_i = this->rho_[index_i];
    Vecd vel_i = this->vel_[index_i];
    for (size_t k = 0; k < FluidWallData::contact_configuration_.size(); ++k)
    {
        StdLargeVec<Vecd> &n_k = *(this->wall_n_[k]);
        StdLargeVec<Matd> &B_k = *(this->wall_B_[k]);
        Neighborhood &wall_neighborhood = (*FluidWallData::contact_configuration_[k])[index_i];
        for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
        {
            size_t index_j = wall_neighborhood.j_[n];
            Vecd &e_ij = wall_neighborhood.e_ij_[n];
            Real dW_ijV_j = wall_neighborhood.dW_ijV_j_[n];

            Vecd vel_in_wall = -this->vel_[index_i];
            Real p_in_wall = this->p_[index_i];
            Real rho_in_wall = this->rho_[index_i];

            Real p_jump = p_i - p_in_wall;
            Real u_jump = (vel_i - vel_in_wall).dot(-e_ij);

            Real rho_ave = 2 * rho_i * rho_in_wall / (rho_i + rho_in_wall);

            Vecd rho_vel_ave = 0.5 * (rho_i * vel_i + rho_in_wall * vel_in_wall);
            //Vecd rho_vel_ave = 0.5 * (rho_i * B_k[index_j] * vel_i  + rho_in_wall * this->B_[index_i] * vel_in_wall);

            Vecd rho_dissipation = rho_ave * this->riemann_solver_.DissipativeUJump(p_jump, u_jump) * (-e_ij);
            density_change_rate -= 2 * (rho_vel_ave + rho_dissipation).dot(e_ij) * dW_ijV_j;
        }
    }
};
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
