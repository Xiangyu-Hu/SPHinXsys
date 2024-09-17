#include "continuum_integration.hpp"

namespace SPH
{
//=================================================================================================//
namespace continuum_dynamics
{
ContinuumInitialCondition::ContinuumInitialCondition(SPHBody &sph_body)
    : LocalDynamics(sph_body),
      pos_(particles_->getVariableDataByName<Vecd>("Position")),
      vel_(particles_->registerStateVariable<Vecd>("Velocity")),
      stress_tensor_3D_(particles_->registerStateVariable<Mat3d>("StressTensor3D")) {}
//=================================================================================================//
AcousticTimeStep::AcousticTimeStep(SPHBody &sph_body, Real acousticCFL)
    : LocalDynamicsReduce<ReduceMax>(sph_body),
      fluid_(DynamicCast<Fluid>(this, particles_->getBaseMaterial())),
      rho_(particles_->getVariableDataByName<Real>("Density")),
      p_(particles_->getVariableDataByName<Real>("Pressure")),
      vel_(particles_->getVariableDataByName<Vecd>("Velocity")),
      smoothing_length_min_(sph_body.sph_adaptation_->MinimumSmoothingLength()),
      acousticCFL_(acousticCFL) {}
//=================================================================================================//
Real AcousticTimeStep::reduce(size_t index_i, Real dt)
{
    return fluid_.getSoundSpeed(p_[index_i], rho_[index_i]) + vel_[index_i].norm();
}
//=================================================================================================//
Real AcousticTimeStep::outputResult(Real reduced_value)
{
    return acousticCFL_ * smoothing_length_min_ / (reduced_value + TinyReal);
}
//=================================================================================================//
ShearAccelerationRelaxation::ShearAccelerationRelaxation(BaseInnerRelation &inner_relation)
    : fluid_dynamics::BaseIntegration<DataDelegateInner>(inner_relation),
      continuum_(DynamicCast<GeneralContinuum>(this, particles_->getBaseMaterial())),
      G_(continuum_.getShearModulus(continuum_.getYoungsModulus(), continuum_.getPoissonRatio())),
      smoothing_length_(sph_body_.sph_adaptation_->ReferenceSmoothingLength()),
      acc_shear_(particles_->registerStateVariable<Vecd>("AccelerationByShear"))
{
    particles_->addVariableToSort<Vecd>("AccelerationByShear");
}
//=================================================================================================//
void ShearAccelerationRelaxation::interaction(size_t index_i, Real dt)
{
    Real rho_i = rho_[index_i];
    Vecd acceleration = Vecd::Zero();
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Real r_ij = inner_neighborhood.r_ij_[n];
        Real dW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j];
        Vecd &e_ij = inner_neighborhood.e_ij_[n];
        Real eta_ij = 2 * (0.7 * (Real)Dimensions + 2.1) * (vel_[index_i] - vel_[index_j]).dot(e_ij) / (r_ij + TinyReal);
        acceleration += eta_ij * dW_ijV_j * e_ij;
    }
    acc_shear_[index_i] += G_ * acceleration * dt / rho_i;
}
//=================================================================================================//
ShearStressRelaxation::ShearStressRelaxation(BaseInnerRelation &inner_relation)
    : fluid_dynamics::BaseIntegration<DataDelegateInner>(inner_relation),
      continuum_(DynamicCast<GeneralContinuum>(this, particles_->getBaseMaterial())),
      shear_stress_(particles_->registerStateVariable<Matd>("ShearStress")),
      shear_stress_rate_(particles_->registerStateVariable<Matd>("ShearStressRate")),
      velocity_gradient_(particles_->registerStateVariable<Matd>("VelocityGradient")),
      strain_tensor_(particles_->registerStateVariable<Matd>("StrainTensor")),
      strain_tensor_rate_(particles_->registerStateVariable<Matd>("StrainTensorRate")),
      B_(particles_->getVariableDataByName<Matd>("LinearGradientCorrectionMatrix")),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure"))
{
    particles_->addVariableToSort<Matd>("ShearStress");
    particles_->addVariableToSort<Matd>("ShearStressRate");
    particles_->addVariableToSort<Matd>("VelocityGradient");
    particles_->addVariableToSort<Matd>("StrainTensor");
    particles_->addVariableToSort<Matd>("StrainTensorRate");
}
//====================================================================================//
void ShearStressRelaxation::initialization(size_t index_i, Real dt)
{
    strain_tensor_[index_i] += strain_tensor_rate_[index_i] * 0.5 * dt;
    shear_stress_[index_i] += shear_stress_rate_[index_i] * 0.5 * dt;
}
//====================================================================================//
void ShearStressRelaxation::interaction(size_t index_i, Real dt)
{
    Matd velocity_gradient = Matd::Zero();
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Real dW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_i];
        Vecd &e_ij = inner_neighborhood.e_ij_[n];
        Vecd v_ij = vel_[index_i] - vel_[index_j];
        velocity_gradient -= v_ij * (B_[index_i] * e_ij * dW_ijV_j).transpose();
    }
    velocity_gradient_[index_i] = velocity_gradient;
    strain_tensor_rate_[index_i] = 0.5 * (velocity_gradient + velocity_gradient.transpose());
    strain_tensor_[index_i] += strain_tensor_rate_[index_i] * 0.5 * dt;
}
//====================================================================================//
void ShearStressRelaxation::update(size_t index_i, Real dt)
{
    shear_stress_rate_[index_i] = continuum_.ConstitutiveRelationShearStress(velocity_gradient_[index_i], shear_stress_[index_i]);
    shear_stress_[index_i] += shear_stress_rate_[index_i] * dt * 0.5;
}
//====================================================================================//
StressDiffusion::StressDiffusion(BaseInnerRelation &inner_relation)
    : BasePlasticIntegration<DataDelegateInner>(inner_relation),
      phi_(DynamicCast<PlasticContinuum>(this, plastic_continuum_).getFrictionAngle()),
      smoothing_length_(sph_body_.sph_adaptation_->ReferenceSmoothingLength()),
      sound_speed_(plastic_continuum_.ReferenceSoundSpeed()) {}
//====================================================================================//
void StressDiffusion::interaction(size_t index_i, Real dt)
{
    Vecd acc_prior_i = force_prior_[index_i] / mass_[index_i];
    Real gravity = abs(acc_prior_i(1, 0));
    Real density = plastic_continuum_.getDensity();
    Mat3d diffusion_stress_rate_ = Mat3d::Zero();
    Mat3d diffusion_stress_ = Mat3d::Zero();
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Real r_ij = inner_neighborhood.r_ij_[n];
        Real dW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j];
        Real y_ij = pos_[index_i](1, 0) - pos_[index_j](1, 0);
        diffusion_stress_ = stress_tensor_3D_[index_i] - stress_tensor_3D_[index_j];
        diffusion_stress_(0, 0) -= (1 - sin(phi_)) * density * gravity * y_ij;
        diffusion_stress_(1, 1) -= density * gravity * y_ij;
        diffusion_stress_(2, 2) -= (1 - sin(phi_)) * density * gravity * y_ij;
        diffusion_stress_rate_ += 2 * zeta_ * smoothing_length_ * sound_speed_ *
                                  diffusion_stress_ * r_ij * dW_ijV_j / (r_ij * r_ij + 0.01 * smoothing_length_);
    }
    stress_rate_3D_[index_i] = diffusion_stress_rate_;
}
//====================================================================================//
ShearStressRelaxationHourglassControl1stHalf ::
    ShearStressRelaxationHourglassControl1stHalf(BaseInnerRelation &inner_relation, Real xi)
    : fluid_dynamics::BaseIntegration<DataDelegateInner>(inner_relation),
      continuum_(DynamicCast<GeneralContinuum>(this, particles_->getBaseMaterial())),
      shear_stress_(particles_->registerStateVariable<Matd>("ShearStress")),
      shear_stress_rate_(particles_->registerStateVariable<Matd>("ShearStressRate")),
      velocity_gradient_(particles_->registerStateVariable<Matd>("VelocityGradient")),
      strain_tensor_(particles_->registerStateVariable<Matd>("StrainTensor")),
      strain_tensor_rate_(particles_->registerStateVariable<Matd>("StrainTensorRate")),
      B_(particles_->getVariableDataByName<Matd>("LinearGradientCorrectionMatrix")),
      scale_penalty_force_(particles_->registerStateVariable<Real>("ScalePenaltyForce")), xi_(xi)
{
    particles_->addVariableToSort<Matd>("ShearStress");
    particles_->addVariableToSort<Matd>("ShearStressRate");
    particles_->addVariableToSort<Matd>("VelocityGradient");
    particles_->addVariableToSort<Matd>("StrainTensor");
    particles_->addVariableToSort<Matd>("StrainTensorRate");
    particles_->addVariableToSort<Real>("ScalePenaltyForce");
}
//====================================================================================//
void ShearStressRelaxationHourglassControl1stHalf::interaction(size_t index_i, Real dt)
{
    Matd velocity_gradient = Matd::Zero();
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    Vecd vel_i = vel_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Real dW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j];
        Vecd &e_ij = inner_neighborhood.e_ij_[n];
        Matd velocity_gradient_ij;
        velocity_gradient_ij = -(vel_i - vel_[index_j]) * (B_[index_i] * e_ij * dW_ijV_j).transpose();
        velocity_gradient += velocity_gradient_ij;
    }
    velocity_gradient_[index_i] = velocity_gradient;
}
//====================================================================================//
void ShearStressRelaxationHourglassControl1stHalf::update(size_t index_i, Real dt)
{
    scale_penalty_force_[index_i] = xi_;
    shear_stress_rate_[index_i] = continuum_.ConstitutiveRelationShearStress(velocity_gradient_[index_i], shear_stress_[index_i]);
    shear_stress_[index_i] += shear_stress_rate_[index_i] * dt;
    strain_tensor_rate_[index_i] = 0.5 * (velocity_gradient_[index_i] + velocity_gradient_[index_i].transpose());
    strain_tensor_[index_i] += strain_tensor_rate_[index_i] * dt;
}
//====================================================================================//
ShearStressRelaxationHourglassControl2ndHalf ::
    ShearStressRelaxationHourglassControl2ndHalf(BaseInnerRelation &inner_relation)
    : fluid_dynamics::BaseIntegration<DataDelegateInner>(inner_relation),
      continuum_(DynamicCast<GeneralContinuum>(this, particles_->getBaseMaterial())),
      shear_stress_(particles_->getVariableDataByName<Matd>("ShearStress")),
      velocity_gradient_(particles_->getVariableDataByName<Matd>("VelocityGradient")),
      acc_shear_(particles_->registerStateVariable<Vecd>("AccelerationByShear")),
      acc_hourglass_(particles_->registerStateVariable<Vecd>("AccelerationHourglass")),
      scale_penalty_force_(particles_->getVariableDataByName<Real>("ScalePenaltyForce")),
      G_(continuum_.getShearModulus(continuum_.getYoungsModulus(), continuum_.getPoissonRatio()))
{
    particles_->addVariableToSort<Vecd>("AccelerationByShear");
    particles_->addVariableToSort<Vecd>("AccelerationHourglass");
}
//====================================================================================//
void ShearStressRelaxationHourglassControl2ndHalf::interaction(size_t index_i, Real dt)
{
    Real rho_i = rho_[index_i];
    Matd shear_stress_i = shear_stress_[index_i];
    Vecd acceleration = Vecd::Zero();
    Vecd acceleration_hourglass = Vecd::Zero();
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Real dW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j];
        Vecd &e_ij = inner_neighborhood.e_ij_[n];
        acceleration += ((shear_stress_i + shear_stress_[index_j]) / rho_i) * dW_ijV_j * e_ij;
        Vecd v_ij = vel_[index_i] - vel_[index_j];
        Real r_ij = inner_neighborhood.r_ij_[n];
        Vecd v_ij_correction = v_ij - 0.5 * (velocity_gradient_[index_i] + velocity_gradient_[index_j]) * r_ij * e_ij;
        acceleration_hourglass += 0.5 * (scale_penalty_force_[index_i] + scale_penalty_force_[index_j]) * G_ * v_ij_correction * dW_ijV_j * dt / (rho_i * r_ij);
    }
    acc_hourglass_[index_i] += acceleration_hourglass;
    acc_shear_[index_i] = acceleration + acc_hourglass_[index_i];
}
//====================================================================================//
ShearStressRelaxationHourglassControl1stHalfJ2Plasticity ::
    ShearStressRelaxationHourglassControl1stHalfJ2Plasticity(BaseInnerRelation &inner_relation, Real xi)
    : ShearStressRelaxationHourglassControl1stHalf(inner_relation, xi),
      J2_plasticity_(DynamicCast<J2Plasticity>(this, particles_->getBaseMaterial())),
      hardening_factor_(particles_->registerStateVariable<Real>("HardeningFactor"))
{
    particles_->addVariableToSort<Real>("HardeningFactor");
}
//====================================================================================//
void ShearStressRelaxationHourglassControl1stHalfJ2Plasticity::update(size_t index_i, Real dt)
{
    shear_stress_rate_[index_i] = J2_plasticity_.ConstitutiveRelationShearStress(velocity_gradient_[index_i], shear_stress_[index_i], hardening_factor_[index_i]);
    Matd shear_stress_try = shear_stress_[index_i] + shear_stress_rate_[index_i] * dt;
    Real hardening_factor_increment = J2_plasticity_.HardeningFactorRate(shear_stress_try, hardening_factor_[index_i]);
    hardening_factor_[index_i] += sqrt(2.0 / 3.0) * hardening_factor_increment;
    scale_penalty_force_[index_i] = xi_ * J2_plasticity_.ScalePenaltyForce(shear_stress_try, hardening_factor_[index_i]);
    shear_stress_[index_i] = J2_plasticity_.ReturnMappingShearStress(shear_stress_try, hardening_factor_[index_i]);
    Matd strain_rate = 0.5 * (velocity_gradient_[index_i] + velocity_gradient_[index_i].transpose());
    strain_tensor_[index_i] += strain_rate * dt;
}
} // namespace continuum_dynamics
} // namespace SPH