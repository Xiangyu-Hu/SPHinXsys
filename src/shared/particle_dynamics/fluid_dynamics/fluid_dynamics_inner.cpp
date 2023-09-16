#include "fluid_dynamics_inner.hpp"

namespace SPH
{
//=====================================================================================================//
namespace fluid_dynamics
{
//=================================================================================================//
FluidInitialCondition::
    FluidInitialCondition(SPHBody &sph_body)
    : LocalDynamics(sph_body), FluidDataSimple(sph_body),
      pos_(particles_->pos_), vel_(particles_->vel_) {}
//=================================================================================================//
void BaseDensitySummationInner::update(size_t index_i, Real dt)
{
    rho_[index_i] = rho_sum_[index_i];
}
//=================================================================================================//
DensitySummationInner::DensitySummationInner(BaseInnerRelation &inner_relation)
    : BaseDensitySummationInner(inner_relation),
      W0_(sph_body_.sph_adaptation_->getKernel()->W0(ZeroVecd)) {}
//=================================================================================================//
void DensitySummationInner::interaction(size_t index_i, Real dt)
{
    Real sigma = W0_;
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        sigma += inner_neighborhood.W_ij_[n];

    rho_sum_[index_i] = sigma * rho0_ * inv_sigma0_;
}
//=================================================================================================//
DensitySummationInnerAdaptive::DensitySummationInnerAdaptive(BaseInnerRelation &inner_relation)
    : BaseDensitySummationInner(inner_relation),
      sph_adaptation_(*sph_body_.sph_adaptation_),
      kernel_(*sph_adaptation_.getKernel()),
      h_ratio_(*particles_->getVariableByName<Real>("SmoothingLengthRatio")) {}
//=================================================================================================//
void DensitySummationInnerAdaptive::interaction(size_t index_i, Real dt)
{
    Real sigma_i = mass_[index_i] * kernel_.W0(h_ratio_[index_i], ZeroVecd);
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        sigma_i += inner_neighborhood.W_ij_[n] * mass_[inner_neighborhood.j_[n]];

    rho_sum_[index_i] = sigma_i * rho0_ * inv_sigma0_ / mass_[index_i] /
                        sph_adaptation_.NumberDensityScaleFactor(h_ratio_[index_i]);
}
//=================================================================================================//
void ViscousAccelerationInner::interaction(size_t index_i, Real dt)
{
    Vecd acceleration = Vecd::Zero();
    Vecd vel_derivative = Vecd::Zero();
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];

        // viscous force
        vel_derivative = (vel_[index_i] - vel_[index_j]) / (inner_neighborhood.r_ij_[n] + 0.01 * smoothing_length_);
        acceleration += 2.0 * mu_ * vel_derivative * inner_neighborhood.dW_ijV_j_[n];
    }

    acc_prior_[index_i] += acceleration / rho_[index_i];
}
//=================================================================================================//
void AngularConservativeViscousAccelerationInner::interaction(size_t index_i, Real dt)
{
    Vecd acceleration = Vecd::Zero();
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Vecd &e_ij = inner_neighborhood.e_ij_[n];
        Real r_ij = inner_neighborhood.r_ij_[n];

        /** The following viscous force is given in Monaghan 2005 (Rep. Prog. Phys.), it seems that
         * this formulation is more accurate than the previous one for Taylor-Green-Vortex flow. */
        Real v_r_ij = (vel_[index_i] - vel_[index_j]).dot(r_ij * e_ij);
        Real eta_ij = 8.0 * mu_ * v_r_ij / (r_ij * r_ij + 0.01 * smoothing_length_);
        acceleration += eta_ij * inner_neighborhood.dW_ijV_j_[n] * e_ij;
    }

    acc_prior_[index_i] += acceleration / rho_[index_i];
}
//=================================================================================================//
AcousticTimeStepSize::AcousticTimeStepSize(SPHBody &sph_body, Real acousticCFL)
    : LocalDynamicsReduce<Real, ReduceMax>(sph_body, Real(0)),
      FluidDataSimple(sph_body), fluid_(DynamicCast<Fluid>(this, particles_->getBaseMaterial())),
      rho_(particles_->rho_), p_(*particles_->getVariableByName<Real>("Pressure")), vel_(particles_->vel_),
      smoothing_length_min_(sph_body.sph_adaptation_->MinimumSmoothingLength()),
      acousticCFL_(acousticCFL) {}
//=================================================================================================//
Real AcousticTimeStepSize::reduce(size_t index_i, Real dt)
{
    return fluid_.getSoundSpeed(p_[index_i], rho_[index_i]) + vel_[index_i].norm();
}
//=================================================================================================//
Real AcousticTimeStepSize::outputResult(Real reduced_value)
{
    // since the particle does not change its configuration in pressure relaxation step
    // I chose a time-step size according to Eulerian method
    return acousticCFL_ * smoothing_length_min_ / (reduced_value + TinyReal);
}
//=================================================================================================//
AdvectionTimeStepSizeForImplicitViscosity::
    AdvectionTimeStepSizeForImplicitViscosity(SPHBody &sph_body, Real U_ref, Real advectionCFL)
    : LocalDynamicsReduce<Real, ReduceMax>(sph_body, U_ref * U_ref),
      FluidDataSimple(sph_body), vel_(particles_->vel_),
      smoothing_length_min_(sph_body.sph_adaptation_->MinimumSmoothingLength()),
      speed_ref_(U_ref), advectionCFL_(advectionCFL) {}
//=================================================================================================//
Real AdvectionTimeStepSizeForImplicitViscosity::reduce(size_t index_i, Real dt)
{
    return vel_[index_i].squaredNorm();
}
//=================================================================================================//
Real AdvectionTimeStepSizeForImplicitViscosity::outputResult(Real reduced_value)
{
    Real speed_max = sqrt(reduced_value);
    return advectionCFL_ * smoothing_length_min_ / (SMAX(speed_max, speed_ref_) + TinyReal);
}
//=================================================================================================//
AdvectionTimeStepSize::AdvectionTimeStepSize(SPHBody &sph_body, Real U_ref, Real advectionCFL)
    : AdvectionTimeStepSizeForImplicitViscosity(sph_body, U_ref, advectionCFL),
      fluid_(DynamicCast<Fluid>(this, particles_->getBaseMaterial()))
{
    Real viscous_speed = fluid_.ReferenceViscosity() / fluid_.ReferenceDensity() / smoothing_length_min_;
    speed_ref_ = SMAX(viscous_speed, speed_ref_);
}
//=================================================================================================//
Real AdvectionTimeStepSize::reduce(size_t index_i, Real dt)
{
    return AdvectionTimeStepSizeForImplicitViscosity::reduce(index_i, dt);
}
//=================================================================================================//
VorticityInner::VorticityInner(BaseInnerRelation &inner_relation)
    : LocalDynamics(inner_relation.getSPHBody()), FluidDataInner(inner_relation),
      vel_(particles_->vel_)
{
    particles_->registerVariable(vorticity_, "VorticityInner");
    particles_->addVariableToWrite<AngularVecd>("VorticityInner");
}
//=================================================================================================//
void VorticityInner::interaction(size_t index_i, Real dt)
{
    AngularVecd vorticity = ZeroData<AngularVecd>::value;
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];

        Vecd vel_diff = vel_[index_i] - vel_[index_j];
        vorticity += getCrossProduct(vel_diff, inner_neighborhood.e_ij_[n]) * inner_neighborhood.dW_ijV_j_[n];
    }

    vorticity_[index_i] = vorticity;
}
//=================================================================================================//
Oldroyd_BIntegration1stHalf ::Oldroyd_BIntegration1stHalf(BaseInnerRelation &inner_relation)
    : Integration1stHalfInnerDissipativeRiemann(inner_relation)
{
    particles_->registerVariable(tau_, "ElasticStress");
    particles_->registerVariable(dtau_dt_, "ElasticStressChangeRate");
    particles_->registerSortableVariable<Matd>("ElasticStress");
    particles_->addVariableToRestart<Matd>("ElasticStress");
}
//=================================================================================================//
void Oldroyd_BIntegration1stHalf::initialization(size_t index_i, Real dt)
{
    Integration1stHalfInnerDissipativeRiemann::initialization(index_i, dt);

    tau_[index_i] += dtau_dt_[index_i] * dt * 0.5;
}
//=================================================================================================//
void Oldroyd_BIntegration1stHalf::interaction(size_t index_i, Real dt)
{
    Integration1stHalfInnerDissipativeRiemann::interaction(index_i, dt);

    Vecd acceleration = Vecd::Zero();
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Vecd nablaW_ijV_j = inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];

        // elastic force
        acceleration += (tau_[index_i] + tau_[index_j]) * nablaW_ijV_j;
    }

    acc_[index_i] += acceleration / rho_[index_i];
}
//=================================================================================================//
Oldroyd_BIntegration2ndHalf::
    Oldroyd_BIntegration2ndHalf(BaseInnerRelation &inner_relation)
    : Integration2ndHalfInnerDissipativeRiemann(inner_relation),
      oldroyd_b_fluid_(DynamicCast<Oldroyd_B_Fluid>(this, particles_->getBaseMaterial())),
      tau_(*particles_->getVariableByName<Matd>("ElasticStress")),
      dtau_dt_(*particles_->getVariableByName<Matd>("ElasticStressChangeRate"))
{
    mu_p_ = oldroyd_b_fluid_.ReferencePolymericViscosity();
    lambda_ = oldroyd_b_fluid_.getReferenceRelaxationTime();
}
//=================================================================================================//
void Oldroyd_BIntegration2ndHalf::update(size_t index_i, Real dt)
{
    Integration2ndHalfInnerDissipativeRiemann::update(index_i, dt);

    tau_[index_i] += dtau_dt_[index_i] * dt * 0.5;
}
//=================================================================================================//
void Oldroyd_BIntegration2ndHalf::
    interaction(size_t index_i, Real dt)
{
    Integration2ndHalfInnerDissipativeRiemann::interaction(index_i, dt);

    Matd tau_i = tau_[index_i];
    Matd stress_rate = Matd::Zero();
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Vecd nablaW_ijV_j = inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];

        Matd velocity_gradient = -(vel_[index_i] - vel_[index_j]) * nablaW_ijV_j.transpose();
        stress_rate += velocity_gradient.transpose() * tau_i + tau_i * velocity_gradient - tau_i / lambda_ +
                       (velocity_gradient.transpose() + velocity_gradient) * mu_p_ / lambda_;
    }

    dtau_dt_[index_i] = stress_rate;
}
//=================================================================================================//
} // namespace fluid_dynamics
  //=====================================================================================================//
} // namespace SPH
  //=========================================================================================================//