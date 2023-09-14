
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
NoRiemannSolverInCompressobleEulerianMethod ::NoRiemannSolverInCompressobleEulerianMethod(CompressibleFluid &compressible_fluid_i, CompressibleFluid &compressible_fluid_j)
    : compressible_fluid_i_(compressible_fluid_i), compressible_fluid_j_(compressible_fluid_j){};
//=================================================================================================//
CompressibleFluidStarState NoRiemannSolverInCompressobleEulerianMethod::
    getInterfaceState(const CompressibleFluidState &state_i, const CompressibleFluidState &state_j, const Vecd &e_ij)
{
    Real p_star = 0.5 * (state_i.p_ + state_j.p_);
    Vecd v_star = 0.5 * (state_i.vel_ + state_j.vel_);
    Real rho_star = 0.5 * (state_i.rho_ + state_j.rho_);
    Real energy_star = 0.5 * (state_i.E_ + state_j.E_);

    return CompressibleFluidStarState(rho_star, v_star, p_star, energy_star);
}
//=================================================================================================//
HLLCRiemannSolver::HLLCRiemannSolver(CompressibleFluid &compressible_fluid_i, CompressibleFluid &compressible_fluid_j, Real limiter_parameter)
    : compressible_fluid_i_(compressible_fluid_i), compressible_fluid_j_(compressible_fluid_j){};
//=================================================================================================//
CompressibleFluidStarState HLLCRiemannSolver::
    getInterfaceState(const CompressibleFluidState &state_i, const CompressibleFluidState &state_j, const Vecd &e_ij)
{
    Real ul = -e_ij.dot(state_i.vel_);
    Real ur = -e_ij.dot(state_j.vel_);
    Real s_l = ul - compressible_fluid_i_.getSoundSpeed(state_i.p_, state_i.rho_);
    Real s_r = ur + compressible_fluid_j_.getSoundSpeed(state_j.p_, state_j.rho_);
    Real s_star = (state_j.rho_ * ur * (s_r - ur) + state_i.rho_ * ul * (ul - s_l) + state_i.p_ - state_j.p_) / (state_j.rho_ * (s_r - ur) + state_i.rho_ * (ul - s_l));
    Real p_star = 0.0;
    Vecd v_star = Vecd::Zero();
    Real rho_star = 0.0;
    Real energy_star = 0.0;
    if (0.0 < s_l)
    {
        p_star = state_i.p_;
        v_star = state_i.vel_;
        rho_star = state_i.rho_;
        energy_star = state_i.E_;
    }
    if (s_l <= 0.0 && 0.0 <= s_star)
    {
        p_star = state_i.p_ + state_i.rho_ * (s_l - ul) * (s_star - ul);
        v_star = state_i.vel_ - e_ij * (s_star - ul);
        rho_star = state_i.rho_ * (s_l - ul) / (s_l - s_star);
        energy_star = state_i.rho_ * (s_l - ul) / (s_l - s_star) * (state_i.E_ / state_i.rho_ + (s_star - ul) * (s_star + state_i.p_ / state_i.rho_ / (s_l - ul)));
    }
    if (s_star <= 0.0 && 0.0 <= s_r)
    {
        p_star = state_i.p_ + state_i.rho_ * (s_l - ul) * (s_star - ul);
        v_star = state_j.vel_ - e_ij * (s_star - ur);
        rho_star = state_j.rho_ * (s_r - ur) / (s_r - s_star);
        energy_star = state_j.rho_ * (s_r - ur) / (s_r - s_star) * (state_j.E_ / state_j.rho_ + (s_star - ur) * (s_star + state_j.p_ / state_j.rho_ / (s_r - ur)));
    }
    if (s_r < 0.0)
    {
        p_star = state_j.p_;
        v_star = state_j.vel_;
        rho_star = state_j.rho_;
        energy_star = state_j.E_;
    }
    return CompressibleFluidStarState(rho_star, v_star, p_star, energy_star);
}
//=================================================================================================//
HLLCWithLimiterRiemannSolver::HLLCWithLimiterRiemannSolver(CompressibleFluid &compressible_fluid_i, CompressibleFluid &compressible_fluid_j, Real limiter_parameter)
    : compressible_fluid_i_(compressible_fluid_i), compressible_fluid_j_(compressible_fluid_j), limiter_parameter_(limiter_parameter){};
//=================================================================================================//
CompressibleFluidStarState HLLCWithLimiterRiemannSolver::
    getInterfaceState(const CompressibleFluidState &state_i, const CompressibleFluidState &state_j, const Vecd &e_ij)
{
    Real ul = -e_ij.dot(state_i.vel_);
    Real ur = -e_ij.dot(state_j.vel_);
    Real s_l = ul - compressible_fluid_i_.getSoundSpeed(state_i.p_, state_i.rho_);
    Real s_r = ur + compressible_fluid_j_.getSoundSpeed(state_j.p_, state_j.rho_);
    Real rhol_cl = compressible_fluid_i_.getSoundSpeed(state_i.p_, state_i.rho_) * state_i.rho_;
    Real rhor_cr = compressible_fluid_j_.getSoundSpeed(state_j.p_, state_j.rho_) * state_j.rho_;
    Real clr = (rhol_cl + rhor_cr) / (state_i.rho_ + state_j.rho_);
    Real s_star = (state_j.p_ - state_i.p_) * pow(SMIN(limiter_parameter_ * SMAX((ul - ur) / clr, Real(0)), Real(1)), 2) / (state_i.rho_ * (s_l - ul) - state_j.rho_ * (s_r - ur)) +
                  (state_i.rho_ * (s_l - ul) * ul - state_j.rho_ * (s_r - ur) * ur) / (state_i.rho_ * (s_l - ul) - state_j.rho_ * (s_r - ur));
    Real p_star = 0.0;
    Vecd v_star = Vecd::Zero();
    Real rho_star = 0.0;
    Real energy_star = 0.0;
    if (0.0 < s_l)
    {
        p_star = state_i.p_;
        v_star = state_i.vel_;
        rho_star = state_i.rho_;
        energy_star = state_i.E_;
    }
    if (s_l <= 0.0 && 0.0 <= s_star)
    {
        p_star = 0.5 * (state_i.p_ + state_j.p_) +
                 0.5 * (state_i.rho_ * (s_l - ul) * (s_star - ul) + state_j.rho_ * (s_r - ur) * (s_star - ur)) * SMIN(limiter_parameter_ * SMAX((ul - ur) / clr, Real(0)), Real(1));
        v_star = state_i.vel_ - e_ij * (s_star - ul);
        rho_star = state_i.rho_ * (s_l - ul) / (s_l - s_star);
        energy_star = ((s_l - ul) * state_i.E_ - state_i.p_ * ul + p_star * s_star) / (s_l - s_star);
    }
    if (s_star <= 0.0 && 0.0 <= s_r)
    {
        p_star = 0.5 * (state_i.p_ + state_j.p_) +
                 0.5 * (state_i.rho_ * (s_l - ul) * (s_star - ul) + state_j.rho_ * (s_r - ur) * (s_star - ur)) * SMIN(limiter_parameter_ * SMAX((ul - ur) / clr, Real(0)), Real(1));
        v_star = state_j.vel_ - e_ij * (s_star - ur);
        rho_star = state_j.rho_ * (s_r - ur) / (s_r - s_star);
        energy_star = ((s_r - ur) * state_j.E_ - state_j.p_ * ur + p_star * s_star) / (s_r - s_star);
    }
    if (s_r < 0.0)
    {
        p_star = state_j.p_;
        v_star = state_j.vel_;
        rho_star = state_j.rho_;
        energy_star = state_j.E_;
    }
    return CompressibleFluidStarState(rho_star, v_star, p_star, energy_star);
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
BaseIntegrationInCompressible::BaseIntegrationInCompressible(BaseInnerRelation &inner_relation) : BaseIntegration(inner_relation),
                                                                                                  compressible_fluid_(CompressibleFluid(1.0, 1.4)),
                                                                                                  Vol_(particles_->Vol_), E_(*particles_->getVariableByName<Real>("TotalEnergy")),
                                                                                                  dE_dt_(*particles_->getVariableByName<Real>("TotalEnergyChangeRate")),
                                                                                                  dE_dt_prior_(*particles_->getVariableByName<Real>("OtherEnergyChangeRate")),
                                                                                                  mom_(*particles_->getVariableByName<Vecd>("Momentum")),
                                                                                                  dmom_dt_(*particles_->getVariableByName<Vecd>("MomentumChangeRate")),
                                                                                                  dmom_dt_prior_(*particles_->getVariableByName<Vecd>("OtherMomentumChangeRate")){};
//=================================================================================================//
} // namespace SPH
  //=================================================================================================//