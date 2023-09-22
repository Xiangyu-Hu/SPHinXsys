#include "shear_dynamics.hpp"
namespace SPH
{
namespace continuum_dynamics
{
//=================================================================================================//
ContinuumInitialCondition::ContinuumInitialCondition(SPHBody &sph_body)
    : LocalDynamics(sph_body), PlasticContinuumDataSimple(sph_body),
      pos_(particles_->pos_), vel_(particles_->vel_), stress_tensor_3D_(particles_->stress_tensor_3D_) {}
//=================================================================================================//
ContinuumAcousticTimeStepSize::ContinuumAcousticTimeStepSize(SPHBody &sph_body, Real acousticCFL)
    : fluid_dynamics::AcousticTimeStepSize(sph_body, acousticCFL) {}
//=================================================================================================//
Real ContinuumAcousticTimeStepSize::reduce(size_t index_i, Real dt)
{
    return fluid_.getSoundSpeed(p_[index_i], rho_[index_i]) + vel_[index_i].norm();
}
//=================================================================================================//
Real ContinuumAcousticTimeStepSize::outputResult(Real reduced_value)
{
    return acousticCFL_ * smoothing_length_min_ / (fluid_.ReferenceSoundSpeed() + TinyReal);
}
//=================================================================================================//
ShearStressIntegration::ShearStressIntegration(BaseInnerRelation &inner_relation)
    : BaseShearStressIntegration<ContinuumDataInner>(inner_relation),
      continuum_(particles_->continuum_){};
//=================================================================================================//
void ShearStressIntegration::update(size_t index_i, Real dt)
{
    Matd shear_stress_rate = continuum_.ConstitutiveRelationShearStress(velocity_gradient_[index_i], shear_stress_[index_i]);
    shear_stress_[index_i] += shear_stress_rate * dt;
}
//=================================================================================================//
PlasticShearStressIntegration::PlasticShearStressIntegration(BaseInnerRelation &inner_relation)
    : BaseShearStressIntegration<PlasticContinuumDataInner>(inner_relation),
      plastic_continuum_(particles_->plastic_continuum_),
      stress_tensor_3D_(particles_->stress_tensor_3D_), strain_tensor_3D_(particles_->strain_tensor_3D_),
      stress_rate_3D_(particles_->stress_rate_3D_), strain_rate_3D_(particles_->strain_rate_3D_),
      elastic_strain_tensor_3D_(particles_->elastic_strain_tensor_3D_),
      elastic_strain_rate_3D_(particles_->elastic_strain_rate_3D_),
      E_(plastic_continuum_.getYoungsModulus()), nu_(plastic_continuum_.getPoissonRatio()),
      shear_stress_(*particles_->getVariableByName<Matd>("ShearStress")) {}
//=================================================================================================//
void PlasticShearStressIntegration::update(size_t index_i, Real dt)
{
    Mat3d velocity_gradient = upgradeToMat3d(velocity_gradient_[index_i]);
    stress_rate_3D_[index_i] += plastic_continuum_.ConstitutiveRelation(velocity_gradient, stress_tensor_3D_[index_i]);
    stress_tensor_3D_[index_i] += stress_rate_3D_[index_i] * dt;
    stress_tensor_3D_[index_i] = plastic_continuum_.ReturnMapping(stress_tensor_3D_[index_i]);
    strain_rate_3D_[index_i] = 0.5 * (velocity_gradient + velocity_gradient.transpose());
    strain_tensor_3D_[index_i] += strain_rate_3D_[index_i] * dt * 0.5;
    shear_stress_[index_i] = degradeToMatd(strain_tensor_3D_[index_i]);
    // calculate elastic strain for output visualization
    Mat3d deviatoric_stress = stress_tensor_3D_[index_i] - (1 / 3) * stress_tensor_3D_[index_i].trace() * Mat3d::Identity();
    Real hydrostatic_pressure = (1 / 3) * stress_tensor_3D_[index_i].trace();
    elastic_strain_tensor_3D_[index_i] = deviatoric_stress / (2 * plastic_continuum_.getShearModulus(E_, nu_)) +
                                         hydrostatic_pressure * Mat3d::Identity() / (9 * plastic_continuum_.getBulkModulus(E_, nu_));
}
//=================================================================================================//
void ShearStressAcceleration::interaction(size_t index_i, Real dt)
{
    Real rho_i = rho_[index_i];
    Matd shear_stress_i = shear_stress_[index_i];
    Vecd acceleration = Vecd::Zero();
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Vecd nablaW_ijV_j = inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
        acceleration += rho_[index_j] * (shear_stress_i / (rho_i * rho_i) + shear_stress_[index_j] / (rho_[index_j] * rho_[index_j])) * nablaW_ijV_j;
    }
    acc_prior_[index_i] += acceleration;
}
//=================================================================================================//
void ShearStressAccelerationWithWall::interaction(size_t index_i, Real dt)
{
    Matd shear_stress_i = shear_stress_[index_i];
    Real rho_i = rho_[index_i];
    Real rho_in_wall = rho_i;
    Vecd acceleration = Vecd::Zero();
    for (size_t k = 0; k < contact_configuration_.size(); ++k) // There may be several wall bodies.
    {
        Neighborhood &wall_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
        {
            size_t index_j = wall_neighborhood.j_[n];
            Vecd nablaW_ijV_j = wall_neighborhood.dW_ijV_j_[n] * wall_neighborhood.e_ij_[n];

            Matd shear_stress_in_wall = shear_stress_i;
            acceleration += rho_[index_j] * (shear_stress_i / (rho_i * rho_i) + shear_stress_in_wall / (rho_in_wall * rho_in_wall)) * nablaW_ijV_j;
        }
    }
    acc_prior_[index_i] += acceleration;
}
//=================================================================================================//
ShearAccelerationIntegration::ShearAccelerationIntegration(BaseInnerRelation &inner_relation)
    : LocalDynamics(inner_relation.getSPHBody()), ContinuumDataInner(inner_relation),
      continuum_(particles_->continuum_),
      G_(continuum_.getShearModulus(continuum_.getYoungsModulus(), continuum_.getPoissonRatio())),
      smoothing_length_(sph_body_.sph_adaptation_->ReferenceSmoothingLength()),
      vel_(particles_->vel_), acc_prior_(particles_->acc_prior_), rho_(particles_->rho_)
{
    particles_->registerVariable(acc_shear_, "AccumulatedShearAcceleration");
}
//=================================================================================================//
void ShearAccelerationIntegration::interaction(size_t index_i, Real dt)
{
    Vecd acceleration = Vecd::Zero();
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Real r_ij = inner_neighborhood.r_ij_[n];
        Real dW_ijV_j = inner_neighborhood.dW_ijV_j_[n];
        Vecd &e_ij = inner_neighborhood.e_ij_[n];
        Real eta_ij = 2 * (0.7 * (Real)Dimensions + 2.1) * (vel_[index_i] - vel_[index_j]).dot(e_ij) / (r_ij + TinyReal);
        acceleration += eta_ij * dW_ijV_j * e_ij;
    }
    acc_shear_[index_i] += G_ * acceleration * dt / rho_[index_i];
    acc_prior_[index_i] += acc_shear_[index_i];
}
//=================================================================================================//
StressDiffusion::StressDiffusion(BaseInnerRelation &inner_relation, SharedPtr<Gravity> gravity__ptr, int axis)
    : LocalDynamics(inner_relation.getSPHBody()), PlasticContinuumDataInner(inner_relation),
      plastic_continuum_(particles_->plastic_continuum_),
      axis_(axis), rho0_(plastic_continuum_.ReferenceDensity()),
      gravity_(gravity__ptr->InducedAcceleration()[axis]),
      smoothing_length_(sph_body_.sph_adaptation_->ReferenceSmoothingLength()),
      phi_(plastic_continuum_.getFrictionAngle()),
      diffusion_coeff_(zeta_ * smoothing_length_ * plastic_continuum_.ReferenceSoundSpeed()),
      pos_(particles_->pos_), stress_tensor_3D_(particles_->stress_tensor_3D_),
      stress_rate_3D_(particles_->stress_rate_3D_) {}
//=================================================================================================//
void StressDiffusion::interaction(size_t index_i, Real dt)
{
    Mat3d diffusion_stress_rate_ = Mat3d::Zero();
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Real r_ij = inner_neighborhood.r_ij_[n];
        Real dW_ijV_j = inner_neighborhood.dW_ijV_j_[n];
        Real y_ij = (pos_[index_i] - pos_[index_j])[axis_];
        Mat3d difference = stress_tensor_3D_[index_i] - stress_tensor_3D_[index_j];
        difference(0, 0) -= (1 - sin(phi_)) * rho0_ * gravity_ * y_ij;
        difference(1, 1) -= rho0_ * gravity_ * y_ij;
        difference(2, 2) -= (1 - sin(phi_)) * rho0_ * gravity_ * y_ij;
        diffusion_stress_rate_ += 2.0 * diffusion_coeff_ * difference * dW_ijV_j /
                                  (r_ij + 0.01 * smoothing_length_);
    }
    stress_rate_3D_[index_i] = diffusion_stress_rate_;
}
//=================================================================================================//
} // namespace continuum_dynamics
} // namespace SPH
