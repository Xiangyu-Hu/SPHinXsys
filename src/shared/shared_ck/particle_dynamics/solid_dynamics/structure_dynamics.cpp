#include "structure_dynamics.h"

namespace SPH
{
namespace solid_dynamics
{
//=================================================================================================//
AcousticTimeStepCK::AcousticTimeStepCK(SPHBody &sph_body, Real acousticCFL)
    : LocalDynamicsReduce<ReduceMax>(sph_body), acousticCFL_(acousticCFL),
      h_min_(sph_body.getSPHAdaptation().MinimumSmoothingLength()),
      c0_(DynamicCast<ElasticSolid>(this, sph_body.getBaseMaterial()).ReferenceSoundSpeed()),
      dv_mass_(particles_->getVariableByName<Real>("Mass")),
      dv_vel_(particles_->getVariableByName<Vecd>("Velocity")),
      dv_force_(particles_->getVariableByName<Vecd>("Force")),
      dv_force_prior_(particles_->getVariableByName<Vecd>("ForcePrior")) {}
//=================================================================================================//
AcousticTimeStepCK::FinishDynamics::FinishDynamics(AcousticTimeStepCK &encloser)
    : h_min_(encloser.h_min_), acousticCFL_(encloser.acousticCFL_) {}
//=================================================================================================//
Real AcousticTimeStepCK::FinishDynamics::Result(Real reduced_value)
{
    // since the particle does not change its configuration in the acoustic time steps
    // I chose a time-step size according to Eulerian method
    return acousticCFL_ * h_min_ / (reduced_value + TinyReal);
}
//=================================================================================================//
StructureDynamicsVariables::StructureDynamicsVariables(BaseParticles *particles)
    : dv_rho_(particles->getVariableByName<Real>("Density")),
      dv_mass_(particles->getVariableByName<Real>("Mass")),
      dv_pos_(particles->getVariableByName<Vecd>("Position")),
      dv_vel_(particles->getVariableByName<Vecd>("Velocity")),
      dv_force_(particles->registerStateVariable<Vecd>("Force")),
      dv_B_(particles->getVariableByName<Matd>("LinearCorrectionMatrix")),
      dv_F_(particles->registerStateVariable<Matd>(
          "DeformationGradient", IdentityMatrix<Matd>::value)),
      dv_dF_dt_(particles->registerStateVariable<Matd>("DeformationRate")),
      dv_inverse_F_(particles->registerStateVariable<Matd>(
          "InverseDeformationGradient", IdentityMatrix<Matd>::value)),
      dv_stress_on_particle_(particles->registerStateVariable<Matd>("StressOnParticle")),
      dv_scaling_matrix_(particles->registerStateVariable<Matd>("ScalingMatrix"))
{
    particles->addEvolvingVariable<Vecd>(dv_pos_);
    particles->addEvolvingVariable<Vecd>(dv_vel_);
    particles->addEvolvingVariable<Matd>(dv_F_);
    particles->addEvolvingVariable<Matd>(dv_dF_dt_);
}
//=================================================================================================//
UpdateElasticNormalDirectionCK::UpdateElasticNormalDirectionCK(SPHBody &sph_body)
    : LocalDynamics(sph_body),
      dv_n_(particles_->getVariableByName<Vecd>("NormalDirection")),
      dv_n0_(particles_->registerStateVariableFrom<Vecd>("InitialNormalDirection", "NormalDirection")),
      dv_phi_(particles_->getVariableByName<Real>("SignedDistance")),
      dv_phi0_(particles_->getVariableByName<Real>("InitialSignedDistance")),
      dv_F_(particles_->getVariableByName<Matd>("DeformationGradient")) {}
//=================================================================================================//
UpdateAnisotropicMeasure::UpdateAnisotropicMeasure(SPHBody &sph_body)
    : LocalDynamics(sph_body),
      dv_scaling_(particles_->getVariableByName<Vecd>("AnisotropicScaling")),
      dv_scaling0_(particles_->registerStateVariableFrom<Vecd>("InitialAnisotropicScaling", "AnisotropicScaling")),
      dv_orientation_(particles_->getVariableByName<Vecd>("AnisotropicOrientation")),
      dv_orientation0_(particles_->registerStateVariableFrom<Vecd>("InitialAnisotropicOrientation", "AnisotropicOrientation")),
      dv_F_(particles_->getVariableByName<Matd>("DeformationGradient")) {}
//=================================================================================================//
} // namespace solid_dynamics
} // namespace SPH
