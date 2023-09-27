
#include "common_compressible_eulerian_classes.h"

namespace SPH
{
//=================================================================================================//
EulerianCompressibleTimeStepInitialization::EulerianCompressibleTimeStepInitialization(SPHBody &sph_body, SharedPtr<Gravity> gravity_ptr)
    : TimeStepInitialization(sph_body, gravity_ptr), rho_(particles_->rho_), pos_(particles_->pos_), vel_(particles_->vel_),
      dmom_dt_prior_(*particles_->getVariableByName<Vecd>("OtherMomentumChangeRate")),
      dE_dt_prior_(*particles_->getVariableByName<Real>("OtherEnergyChangeRate")){};
//=================================================================================================//
void EulerianCompressibleTimeStepInitialization::update(size_t index_i, Real dt)
{
    dmom_dt_prior_[index_i] = rho_[index_i] * gravity_->InducedAcceleration(pos_[index_i]);
    dE_dt_prior_[index_i] = rho_[index_i] * (gravity_->InducedAcceleration(pos_[index_i])).dot(vel_[index_i]);
}
//=================================================================================================//
EulerianCompressibleAcousticTimeStepSize::EulerianCompressibleAcousticTimeStepSize(SPHBody &sph_body)
    : AcousticTimeStepSize(sph_body), rho_(particles_->rho_), p_(*particles_->getVariableByName<Real>("Pressure")), vel_(particles_->vel_),
      smoothing_length_(sph_body.sph_adaptation_->ReferenceSmoothingLength()), compressible_fluid_(CompressibleFluid(1.0, 1.4)){};
//=================================================================================================//
Real EulerianCompressibleAcousticTimeStepSize::reduce(size_t index_i, Real dt)
{
    return compressible_fluid_.getSoundSpeed(p_[index_i], rho_[index_i]) + vel_[index_i].norm();
}
//=================================================================================================//
Real EulerianCompressibleAcousticTimeStepSize::outputResult(Real reduced_value)
{
    return 0.6 / Dimensions * smoothing_length_ / (reduced_value + TinyReal);
}
//=================================================================================================//
NoRiemannSolverInCEM ::NoRiemannSolverInCEM(CompressibleFluid &compressible_fluid_i, CompressibleFluid &compressible_fluid_j)
    : compressible_fluid_i_(compressible_fluid_i), compressible_fluid_j_(compressible_fluid_j){};
//=================================================================================================//
CompressibleFluidStarState NoRiemannSolverInCEM::
    getInterfaceState(const CompressibleFluidState &state_i, const CompressibleFluidState &state_j, const Vecd &e_ij)
{
    Real p_star = 0.5 * (state_i.p_ + state_j.p_);
    Vecd v_star = 0.5 * (state_i.vel_ + state_j.vel_);
    Real rho_star = 0.5 * (state_i.rho_ + state_j.rho_);
    Real energy_star = 0.5 * (state_i.E_ + state_j.E_);

    return CompressibleFluidStarState(rho_star, v_star, p_star, energy_star);
}
//=================================================================================================//
Vec2d HLLCRiemannSolver::getSmallestAndLargestWaveSpeeds(const CompressibleFluidState &state_i, const CompressibleFluidState &state_j, const Vecd &e_ij)
{
    Real ul = -e_ij.dot(state_i.vel_);
    Real ur = -e_ij.dot(state_j.vel_);
    Real rhol_cl = compressible_fluid_i_.getSoundSpeed(state_i.p_, state_i.rho_) * state_i.rho_;
    Real rhor_cr = compressible_fluid_j_.getSoundSpeed(state_j.p_, state_j.rho_) * state_j.rho_;
    Real R_lf = state_j.rho_ / state_i.rho_;
    Real u_tlide = (ul + ur * R_lf) / (1.0 + R_lf);
    Real v_tlide = ((state_i.vel_ - ul * (-e_ij)).norm() + (state_j.vel_ - ur * (-e_ij)).norm() * R_lf) / (1.0 + R_lf);
    Real hl = (state_i.E_ + state_i.p_) / state_i.rho_;
    Real hr = (state_j.E_ + state_j.p_) / state_j.rho_;
    Real h_tlide = (hl + hr * R_lf) / (1.0 + R_lf);
    Real sound_tlide = sqrt((1.4 - 1.0) * (h_tlide - 0.5 * (u_tlide * u_tlide + v_tlide * v_tlide)));
    Real s_l = SMIN(ul - compressible_fluid_i_.getSoundSpeed(state_i.p_, state_i.rho_), u_tlide - sound_tlide);
    Real s_r = SMAX(ur + compressible_fluid_j_.getSoundSpeed(state_j.p_, state_j.rho_), u_tlide + sound_tlide);
    return Vec2d(s_l, s_r);
}
//=================================================================================================//
Real HLLCRiemannSolver::getContactWaveSpeed(const CompressibleFluidState &state_i, const CompressibleFluidState &state_j, const Vecd &e_ij)
{
    Real ul = -e_ij.dot(state_i.vel_);
    Real ur = -e_ij.dot(state_j.vel_);
    Real s_l = getSmallestAndLargestWaveSpeeds(state_i, state_j, e_ij)[0];
    Real s_r = getSmallestAndLargestWaveSpeeds(state_i, state_j, e_ij)[1];
    return (state_j.p_ - state_i.p_) / (state_i.rho_ * (s_l - ul) - state_j.rho_ * (s_r - ur)) +
           (state_i.rho_ * (s_l - ul) * ul - state_j.rho_ * (s_r - ur) * ur) / (state_i.rho_ * (s_l - ul) - state_j.rho_ * (s_r - ur));
}
//=================================================================================================//
CompressibleFluidStarState HLLCRiemannSolver::
    getInterfaceState(const CompressibleFluidState &state_i, const CompressibleFluidState &state_j, const Vecd &e_ij)
{
    Real ul = -e_ij.dot(state_i.vel_);
    Real ur = -e_ij.dot(state_j.vel_);
    Real s_l = getSmallestAndLargestWaveSpeeds(state_i, state_j, e_ij)[0];
    Real s_star = getContactWaveSpeed(state_i, state_j, e_ij);
    Real s_r = getSmallestAndLargestWaveSpeeds(state_i, state_j, e_ij)[1];
    Real p_star, rho_star, energy_star;
    Vecd vel_star;
    if (0.0 < s_l){ return CompressibleFluidStarState(state_i.rho_, state_i.vel_, state_i.p_, state_i.E_);}
    if ((s_l <= 0.0 && 0.0 <= s_star) || (s_star <= 0.0 && 0.0 <= s_r))
    {
        p_star = state_i.p_ + state_i.rho_ * (s_l - ul) * (s_star - ul);
        vel_star = (0.0 <= s_star) ? state_i.vel_ - e_ij * (s_star - ul) : state_j.vel_ - e_ij * (s_star - ur);
        rho_star = (0.0 <= s_star) ? state_i.rho_ * (s_l - ul) / (s_l - s_star) : state_j.rho_ * (s_r - ur) / (s_r - s_star);
        energy_star = (0.0 <= s_star) ? state_i.rho_ * (s_l - ul) / (s_l - s_star) * (state_i.E_ / state_i.rho_ + (s_star - ul) * (s_star + state_i.p_ / state_i.rho_ / (s_l - ul)))
                                      : state_j.rho_ * (s_r - ur) / (s_r - s_star) * (state_j.E_ / state_j.rho_ + (s_star - ur) * (s_star + state_j.p_ / state_j.rho_ / (s_r - ur)));
        return CompressibleFluidStarState(rho_star, vel_star, p_star, energy_star);
    }
    if (s_r < 0.0){ return CompressibleFluidStarState(state_j.rho_, state_j.vel_, state_j.p_, state_j.E_);}
}
//=================================================================================================//
Real HLLCWithLimiterRiemannSolver::getContactWaveSpeed(const CompressibleFluidState &state_i, const CompressibleFluidState &state_j, const Vecd &e_ij)
{
    Real ul = -e_ij.dot(state_i.vel_);
    Real ur = -e_ij.dot(state_j.vel_);
    Real s_l = getSmallestAndLargestWaveSpeeds(state_i, state_j, e_ij)[0];
    Real s_r = getSmallestAndLargestWaveSpeeds(state_i, state_j, e_ij)[1];
    Real clr = (compressible_fluid_i_.getSoundSpeed(state_i.p_, state_i.rho_) * state_i.rho_ 
        + compressible_fluid_j_.getSoundSpeed(state_j.p_, state_j.rho_) * state_j.rho_) / (state_i.rho_ + state_j.rho_);
    return (state_j.p_ - state_i.p_) * pow(SMIN(5.0 * SMAX((ul - ur) / clr, Real(0)), Real(1)), 2) 
        / (state_i.rho_ * (s_l - ul) - state_j.rho_ * (s_r - ur)) + (state_i.rho_ * (s_l - ul) * ul - state_j.rho_ * (s_r - ur) * ur) 
        / (state_i.rho_ * (s_l - ul) - state_j.rho_ * (s_r - ur));
}
//=================================================================================================//
CompressibleFluidStarState HLLCWithLimiterRiemannSolver::
    getInterfaceState(const CompressibleFluidState &state_i, const CompressibleFluidState &state_j, const Vecd &e_ij)
{
    Real ul = -e_ij.dot(state_i.vel_);
    Real ur = -e_ij.dot(state_j.vel_);
    Real clr = (compressible_fluid_i_.getSoundSpeed(state_i.p_, state_i.rho_) * state_i.rho_ 
        + compressible_fluid_j_.getSoundSpeed(state_j.p_, state_j.rho_) * state_j.rho_) / (state_i.rho_ + state_j.rho_);
    Real s_l = getSmallestAndLargestWaveSpeeds(state_i, state_j, e_ij)[0];
    Real s_star = getContactWaveSpeed(state_i, state_j, e_ij);
    Real s_r = getSmallestAndLargestWaveSpeeds(state_i, state_j, e_ij)[1];
    Real p_star, rho_star, energy_star;
    Vecd vel_star;
    if (0.0 < s_l){ return CompressibleFluidStarState(state_i.rho_, state_i.vel_, state_i.p_, state_i.E_);}
    if ((s_l <= 0.0 && 0.0 <= s_star) || (s_star <= 0.0 && 0.0 <= s_r))
    {
        p_star = 0.5 * (state_i.p_ + state_j.p_) + 0.5 * (state_i.rho_ * (s_l - ul) * (s_star - ul) 
            + state_j.rho_ * (s_r - ur) * (s_star - ur)) * SMIN(5.0 * SMAX((ul - ur) / clr, Real(0)), Real(1));
        vel_star = (0.0 <= s_star) ? state_i.vel_ - e_ij * (s_star - ul) : state_j.vel_ - e_ij * (s_star - ur);
        rho_star = (0.0 <= s_star) ? state_i.rho_ * (s_l - ul) / (s_l - s_star) : state_j.rho_ * (s_r - ur) / (s_r - s_star);
        energy_star = (0.0 <= s_star) 
            ? state_i.rho_ * (s_l - ul) / (s_l - s_star) * (state_i.E_ / state_i.rho_ + (s_star - ul) * (s_star + state_i.p_ / state_i.rho_ / (s_l - ul)))
            : state_j.rho_ * (s_r - ur) / (s_r - s_star) * (state_j.E_ / state_j.rho_ + (s_star - ur) * (s_star + state_j.p_ / state_j.rho_ / (s_r - ur)));
        return CompressibleFluidStarState(rho_star, vel_star, p_star, energy_star);
    }
    if (s_r < 0.0){ return CompressibleFluidStarState(state_j.rho_, state_j.vel_, state_j.p_, state_j.E_);}
}
//=================================================================================================//
EulerianCompressibleViscousAccelerationInner::EulerianCompressibleViscousAccelerationInner(BaseInnerRelation &inner_relation)
    : ViscousAccelerationInner(inner_relation),
      dE_dt_prior_(*particles_->getVariableByName<Real>("OtherEnergyChangeRate")),
      dmom_dt_prior_(*particles_->getVariableByName<Vecd>("OtherMomentumChangeRate")){};
//=================================================================================================//
void EulerianCompressibleViscousAccelerationInner::interaction(size_t index_i, Real dt)
{
    Real rho_i = rho_[index_i];
    const Vecd &vel_i = vel_[index_i];

    Vecd acceleration = Vecd::Zero();
    Vecd vel_derivative = Vecd::Zero();
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];

        // viscous force
        vel_derivative = (vel_i - vel_[index_j]) / (inner_neighborhood.r_ij_[n] + 0.01 * smoothing_length_);
        acceleration += 2.0 * mu_ * vel_derivative * inner_neighborhood.dW_ijV_j_[n] / rho_i;
    }
    dmom_dt_prior_[index_i] += rho_[index_i] * acceleration;
    dE_dt_prior_[index_i] += rho_[index_i] * acceleration.dot(vel_[index_i]);
}
//=================================================================================================//
BaseIntegrationInCompressible::BaseIntegrationInCompressible(BaseInnerRelation &inner_relation) 
    : BaseIntegration(inner_relation), compressible_fluid_(CompressibleFluid(1.0, 1.4)),
        Vol_(particles_->Vol_), E_(*particles_->getVariableByName<Real>("TotalEnergy")),
        dE_dt_(*particles_->getVariableByName<Real>("TotalEnergyChangeRate")),
        dE_dt_prior_(*particles_->getVariableByName<Real>("OtherEnergyChangeRate")),
        mom_(*particles_->getVariableByName<Vecd>("Momentum")),
        dmom_dt_(*particles_->getVariableByName<Vecd>("MomentumChangeRate")),
        dmom_dt_prior_(*particles_->getVariableByName<Vecd>("OtherMomentumChangeRate")){};
//=================================================================================================//
} // namespace SPH
  //=================================================================================================//