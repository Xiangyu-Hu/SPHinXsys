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
      smoothing_length_min_(sph_body.getSPHAdaptation().MinimumSmoothingLength()),
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
//====================================================================================//
StressDiffusion::StressDiffusion(BaseInnerRelation &inner_relation)
    : BasePlasticIntegration<DataDelegateInner>(inner_relation),
      phi_(DynamicCast<PlasticContinuum>(this, plastic_continuum_).getFrictionAngle()),
      smoothing_length_(sph_body_.getSPHAdaptation().ReferenceSmoothingLength()),
      sound_speed_(plastic_continuum_.ReferenceSoundSpeed()) {}
//====================================================================================//
void StressDiffusion::interaction(size_t index_i, Real dt)
{
    Vecd acc_prior_i = force_prior_[index_i] / mass_[index_i];
    Real gravity = abs(acc_prior_i(1, 0));
    Real density = plastic_continuum_.getDensity();
    Mat3d diffusion_stress_rate = Mat3d::Zero();
    Mat3d diffusion_stress = Mat3d::Zero();
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Real r_ij = inner_neighborhood.r_ij_[n];
        Real dW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j];
        Real y_ij = pos_[index_i](1, 0) - pos_[index_j](1, 0);
        diffusion_stress = stress_tensor_3D_[index_i] - stress_tensor_3D_[index_j];
        diffusion_stress(0, 0) -= (1 - sin(phi_)) * density * gravity * y_ij;
        diffusion_stress(1, 1) -= density * gravity * y_ij;
        diffusion_stress(2, 2) -= (1 - sin(phi_)) * density * gravity * y_ij;
        diffusion_stress_rate += 2 * zeta_ * smoothing_length_ * sound_speed_ *
                                 diffusion_stress * r_ij * dW_ijV_j / (r_ij * r_ij + 0.01 * smoothing_length_);
    }
    stress_rate_3D_[index_i] = diffusion_stress_rate;
}
//====================================================================================//
ShearStressRelaxationHourglassControl1stHalf ::
    ShearStressRelaxationHourglassControl1stHalf(BaseInnerRelation &inner_relation, Real xi)
    : fluid_dynamics::BaseIntegration<DataDelegateInner>(inner_relation),
      continuum_(DynamicCast<GeneralContinuum>(this, particles_->getBaseMaterial())),
      shear_stress_(particles_->registerStateVariable<Matd>("ShearStress")),
      velocity_gradient_(particles_->registerStateVariable<Matd>("VelocityGradient")),
      strain_tensor_(particles_->registerStateVariable<Matd>("StrainTensor")),
      B_(particles_->getVariableDataByName<Matd>("LinearGradientCorrectionMatrix")),
      scale_penalty_force_(particles_->registerStateVariable<Real>("ScalePenaltyForce")), xi_(xi)
{
    particles_->addEvolvingVariable<Matd>("ShearStress");
    particles_->addEvolvingVariable<Matd>("VelocityGradient");
    particles_->addEvolvingVariable<Matd>("StrainTensor");
    particles_->addEvolvingVariable<Real>("ScalePenaltyForce");
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
        Matd velocity_gradient_ij = -(vel_i - vel_[index_j]) * (B_[index_i] * e_ij * dW_ijV_j).transpose();
        velocity_gradient += velocity_gradient_ij;
    }
    velocity_gradient_[index_i] = velocity_gradient;
}
//====================================================================================//
void ShearStressRelaxationHourglassControl1stHalf::update(size_t index_i, Real dt)
{
    Matd shear_stress_rate = continuum_.ConstitutiveRelationShearStress(velocity_gradient_[index_i], shear_stress_[index_i]);
    shear_stress_[index_i] += shear_stress_rate * dt;
    scale_penalty_force_[index_i] = xi_;
    Matd strain_tensor_rate = 0.5 * (velocity_gradient_[index_i] + velocity_gradient_[index_i].transpose());
    strain_tensor_[index_i] += strain_tensor_rate * dt;
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
    particles_->addEvolvingVariable<Vecd>("AccelerationByShear");
    particles_->addEvolvingVariable<Vecd>("AccelerationHourglass");
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
    particles_->addEvolvingVariable<Real>("HardeningFactor");
}
//====================================================================================//
void ShearStressRelaxationHourglassControl1stHalfJ2Plasticity::update(size_t index_i, Real dt)
{
    Matd shear_stress_rate = J2_plasticity_.ConstitutiveRelationShearStressWithHardening(
        velocity_gradient_[index_i], shear_stress_[index_i], hardening_factor_[index_i]);
    Matd shear_stress_try = shear_stress_[index_i] + shear_stress_rate * dt;
    Real hardening_factor_increment = J2_plasticity_.HardeningFactorRate(shear_stress_try, hardening_factor_[index_i]);
    hardening_factor_[index_i] += sqrt(2.0 / 3.0) * hardening_factor_increment;
    scale_penalty_force_[index_i] = xi_ * J2_plasticity_.ScalePenaltyForce(shear_stress_try, hardening_factor_[index_i]);
    shear_stress_[index_i] = J2_plasticity_.ReturnMappingShearStress(shear_stress_try, hardening_factor_[index_i]);
    Matd strain_rate = 0.5 * (velocity_gradient_[index_i] + velocity_gradient_[index_i].transpose());
    strain_tensor_[index_i] += strain_rate * dt;
}
} // namespace continuum_dynamics
} // namespace SPH