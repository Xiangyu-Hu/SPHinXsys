#include "general_eulerian_fluid_dynamics.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
Vec2d NoRiemannSolverGE::getBoundingWaveSpeeds(const FluidState &state_i, const FluidState &state_j, const Vecd e_ij)
{
    Real s_l = -1;
    Real s_r = 1;
    return Vec2d(s_l, s_r);
};
//=================================================================================================//
DissipationState NoRiemannSolverGE::getDissipationState(const FluidState &state_i, const FluidState &state_j, const Vecd e_ij)
{
    return DissipationState(); //default return zero.
}
//=================================================================================================//
Vec2d HLLERiemannSolver::getBoundingWaveSpeeds(const FluidState &state_i, const FluidState &state_j, const Vecd e_ij)
{
    Real ul = -e_ij.dot(state_i.vel_);
    Real ur = -e_ij.dot(state_j.vel_);
    Real cl = this->fluid_i_.getSoundSpeed(state_i.p_, state_i.rho_);
    Real cr = this->fluid_j_.getSoundSpeed(state_j.p_, state_j.rho_);
    Real s_l = SMIN(ul - cl, ur - cr);
    Real s_r = SMAX(ul + cl, ur + cr);
    return Vec2d(s_l, s_r);
};
//=================================================================================================//
DissipationState HLLERiemannSolver::getDissipationState(const FluidState &state_i, const FluidState &state_j, const Vecd e_ij)
{
    Real ul = -e_ij.dot(state_i.vel_);
    Real ur = -e_ij.dot(state_j.vel_);
    Real cl = this->fluid_i_.getSoundSpeed(state_i.p_, state_i.rho_);
    Real cr = this->fluid_j_.getSoundSpeed(state_j.p_, state_j.rho_);
    Real clr = (cl * state_i.rho_ + cr * state_j.rho_) / (state_i.rho_ + state_j.rho_);
    
    Real s_l = SMIN(ul - cl, ur - cr);
    Real s_r = SMAX(ul + cl, ur + cr);
    Matd flux_lm = (state_i.rho_ * state_i.vel_ * state_i.vel_.transpose() + state_i.p_ * Matd::Identity());
    Matd flux_rm = (state_j.rho_ * state_j.vel_ * state_j.vel_.transpose() + state_j.p_ * Matd::Identity());
    Vecd state_lm = state_i.rho_ * state_i.vel_;
    Vecd state_rm = state_j.rho_ * state_j.vel_;
    Vecd flux_ld = (state_i.rho_ * state_i.vel_);
    Vecd flux_rd = (state_j.rho_ * state_j.vel_);
    Real state_ld = state_i.rho_;
    Real state_rd = state_j.rho_;

    Matd momentum_dissipation = SMIN<Real>(1.0 * SMAX<Real>((ul - ur) / clr, Real(0)), Real(1)) *
        (0.5 * (s_r + s_l) / (s_r - s_l) * (flux_lm - flux_rm) + s_r * s_l * (state_rm - state_lm) * (-e_ij).transpose() / (s_r - s_l));
    Vecd density_dissipation = SMIN<Real>(1.0 * SMAX<Real>((ul - ur) / clr, Real(0)), Real(1)) *
        ((0.5 * (s_r + s_l) / (s_r - s_l) * (flux_ld - flux_rd)) + s_r * s_l * (state_rd - state_ld) / (s_r - s_l) * (-e_ij));

    //Matd momentum_dissipation = SMIN<Real>(15 * SMAX<Real>((ul - ur) / clr, Real(0)), Real(1)) * (state_j.rho_ * state_j.vel_ - state_i.rho_ * state_i.vel_) * (-e_ij).transpose();
    //Vecd density_dissipation = SMIN<Real>(15 * SMAX<Real>((ul - ur) / clr, Real(0)), Real(1)) * (state_j.rho_ - state_i.rho_) * (-e_ij);
    return DissipationState(momentum_dissipation, density_dissipation);
};
//=================================================================================================//
SmearedSurfaceIndication::SmearedSurfaceIndication(BaseInnerRelation &inner_relation)
    : LocalDynamics(inner_relation.getSPHBody()), FluidDataInner(inner_relation),
      indicator_(*particles_->getVariableByName<int>("Indicator")),
      smeared_surface_(*particles_->getVariableByName<int>("SmearedSurface")) {}
//=================================================================================================//
void SmearedSurfaceIndication::interaction(size_t index_i, Real dt)
{
    bool is_near_surface = false;
    const Neighborhood &inner_neighborhood = this->inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        if (indicator_[inner_neighborhood.j_[n]] == 1)
        {
            is_near_surface = true;
            break;
        }
    }
    smeared_surface_[index_i] = is_near_surface;
}
//=================================================================================================//
NonReflectiveBoundaryCorrection::NonReflectiveBoundaryCorrection(BaseInnerRelation &inner_relation)
    : LocalDynamics(inner_relation.getSPHBody()), DataDelegateInner<BaseParticles>(inner_relation),
      fluid_(DynamicCast<WeaklyCompressibleFluid>(this, particles_->getBaseMaterial())),
      rho_farfield_(0.0), sound_speed_(0.0), vel_farfield_(Vecd::Zero()),
      rho_(particles_->rho_), p_(*particles_->getVariableByName<Real>("Pressure")),
      Vol_(particles_->Vol_), vel_(particles_->vel_),
      mom_(*particles_->getVariableByName<Vecd>("Momentum")), pos_(particles_->pos_),
      indicator_(*particles_->getVariableByName<int>("Indicator")),
      n_(*particles_->getVariableByName<Vecd>("NormalDirection"))
{
    particles_->registerVariable(inner_weight_summation_, "InnerWeightSummation");
    particles_->registerVariable(rho_average_, "DensityAverage");
    particles_->registerVariable(vel_normal_average_, "VelocityNormalAverage");
    particles_->registerVariable(vel_tangential_average_, "VelocityTangentialAverage");
    particles_->registerVariable(vel_average_, "VelocityAverage");
    particles_->registerVariable(smeared_surface_, "SmearedSurface");
};
//=================================================================================================//
void NonReflectiveBoundaryCorrection::interaction(size_t index_i, Real dt)
{
    if (indicator_[index_i] == 1 || smeared_surface_[index_i] == 1)
    {
        Real velocity_boundary_normal = vel_[index_i].dot(n_[index_i]);
        bool is_inflow = (n_[index_i][0] <= 0.0 || fabs(n_[index_i][1]) > fabs(n_[index_i][0]));
        bool is_subsonic = (fabs(velocity_boundary_normal) < sound_speed_);

        inner_weight_summation_[index_i] = 0.0;
        Real rho_summation = 0.0;
        Real vel_normal_summation = 0.0;
        Vecd vel_tangential_summation = Vecd::Zero();
        Vecd vel_summation = Vecd::Zero();
        size_t total_inner_neighbor_particles = 0;
        const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];
            if (indicator_[index_j] != 1)
            {
                Real W_ij = inner_neighborhood.W_ij_[n];
                inner_weight_summation_[index_i] += W_ij * Vol_[index_j];
                rho_summation += rho_[index_j];
                vel_normal_summation += vel_[index_j].dot(n_[index_i]);
                vel_tangential_summation += vel_[index_j] - vel_[index_j].dot(n_[index_i]) * n_[index_i];
                vel_summation += vel_[index_j];
                total_inner_neighbor_particles += 1;
            }
        }
        rho_average_[index_i] = rho_summation / (total_inner_neighbor_particles + TinyReal);
        vel_normal_average_[index_i] = vel_normal_summation / (total_inner_neighbor_particles + TinyReal);
        vel_tangential_average_[index_i] = vel_tangential_summation / (total_inner_neighbor_particles + TinyReal);

        if (!is_inflow)
        {
            if (!is_subsonic)
            {
                for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
                {
                    size_t index_j = inner_neighborhood.j_[n];
                    if (indicator_[index_j] != 1)
                    {
                        rho_summation += rho_[index_j];
                        vel_summation += vel_[index_j];
                        total_inner_neighbor_particles += 1;
                    }
                }
                rho_average_[index_i] = rho_summation / (total_inner_neighbor_particles + TinyReal);
                vel_average_[index_i] = vel_summation / (total_inner_neighbor_particles + TinyReal);
            }
        }
    }
}
//=================================================================================================//
void NonReflectiveBoundaryCorrection::update(size_t index_i, Real dt)
{
    if (indicator_[index_i] == 1 || smeared_surface_[index_i] == 1)
    {
        Real velocity_farfield_normal = vel_farfield_.dot(n_[index_i]);
        Real velocity_boundary_normal = vel_[index_i].dot(n_[index_i]);
        bool is_inflow = (n_[index_i][0] <= 0.0 || fabs(n_[index_i][1]) > fabs(n_[index_i][0]));
        bool is_subsonic = (fabs(velocity_boundary_normal) < sound_speed_);

        if (is_inflow)
        {
            if (is_subsonic)
            {
                rho_[index_i] = rho_average_[index_i] * inner_weight_summation_[index_i] + rho_farfield_ * (1.0 - inner_weight_summation_[index_i]);
                p_[index_i] = fluid_.getPressure(rho_[index_i]);
                Real vel_normal = vel_normal_average_[index_i] * inner_weight_summation_[index_i] + velocity_farfield_normal * (1.0 - inner_weight_summation_[index_i]);
                vel_[index_i] = vel_normal * n_[index_i] + (vel_farfield_ - velocity_farfield_normal * n_[index_i]);
                mom_[index_i] = rho_[index_i] * vel_[index_i];
            }
            else
            {
                vel_[index_i] = vel_farfield_;
                rho_[index_i] = rho_farfield_;
                mom_[index_i] = rho_[index_i] * vel_[index_i];
            }
        }
        else
        {
            if (is_subsonic)
            {
                rho_[index_i] = rho_average_[index_i] * inner_weight_summation_[index_i] + rho_farfield_ * (1.0 - inner_weight_summation_[index_i]);
                p_[index_i] = fluid_.getPressure(rho_[index_i]);
                Real vel_normal = vel_normal_average_[index_i] * inner_weight_summation_[index_i] + velocity_farfield_normal * (1.0 - inner_weight_summation_[index_i]);
                vel_[index_i] = vel_normal * n_[index_i] + vel_tangential_average_[index_i];
                mom_[index_i] = rho_[index_i] * vel_[index_i];
            }
            else
            {
                rho_[index_i] = rho_average_[index_i];
                vel_[index_i] = vel_average_[index_i];
                mom_[index_i] = rho_[index_i] * vel_[index_i];
            }
        }
    }
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
