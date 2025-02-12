#include "elastic_dynamics.h"
#include "base_general_dynamics.h"

namespace SPH
{
//=========================================================================================================//
namespace solid_dynamics
{
//=================================================================================================//
AcousticTimeStep::AcousticTimeStep(SPHBody &sph_body, Real CFL)
    : LocalDynamicsReduce<ReduceMin>(sph_body),
      CFL_(CFL),
      elastic_solid_(DynamicCast<ElasticSolid>(this, sph_body.getBaseMaterial())),
      vel_(particles_->getVariableDataByName<Vecd>("Velocity")),
      force_(particles_->getVariableDataByName<Vecd>("Force")),
      force_prior_(particles_->getVariableDataByName<Vecd>("ForcePrior")),
      mass_(particles_->getVariableDataByName<Real>("Mass")),
      smoothing_length_min_(sph_body.getSPHAdaptation().MinimumSmoothingLength()),
      c0_(elastic_solid_.ReferenceSoundSpeed()) {}
//=================================================================================================//
Real AcousticTimeStep::reduce(size_t index_i, Real dt)
{
    // since the particle does not change its configuration in pressure relaxation step
    // I chose a time-step size according to Eulerian method
    Real acceleration_norm = ((force_[index_i] + force_prior_[index_i]) / mass_[index_i]).norm();
    return CFL_ * SMIN((Real)sqrt(smoothing_length_min_ / (acceleration_norm + TinyReal)),
                       smoothing_length_min_ / (c0_ + vel_[index_i].norm()));
}
//=================================================================================================//
ElasticDynamicsInitialCondition::ElasticDynamicsInitialCondition(SPHBody &sph_body)
    : LocalDynamics(sph_body),
      pos_(particles_->getVariableDataByName<Vecd>("Position")),
      vel_(particles_->registerStateVariable<Vecd>("Velocity")) {}
//=================================================================================================//
UpdateElasticNormalDirection::UpdateElasticNormalDirection(SPHBody &sph_body)
    : LocalDynamics(sph_body),
      n_(particles_->getVariableDataByName<Vecd>("NormalDirection")),
      n0_(particles_->registerStateVariableFrom<Vecd>("InitialNormalDirection", "NormalDirection")),
      phi_(particles_->getVariableDataByName<Real>("SignedDistance")),
      phi0_(particles_->getVariableDataByName<Real>("InitialSignedDistance")),
      F_(particles_->getVariableDataByName<Matd>("DeformationGradient")) {}
//=================================================================================================//
void UpdateElasticNormalDirection::update(size_t index_i, Real dt)
{
    // Still used polar decomposition to update the normal direction
    n_[index_i] = getRotatedNormalDirection(F_[index_i], n0_[index_i]);
    // Nanson's relation is used to update the distance to surface
    Vecd current_normal = F_[index_i].inverse().transpose() * n0_[index_i];
    phi_[index_i] = phi0_[index_i] / (current_normal.norm() + SqrtEps); // todo: check this
}
//=================================================================================================//
DeformationGradientBySummation::
    DeformationGradientBySummation(BaseInnerRelation &inner_relation)
    : LocalDynamics(inner_relation.getSPHBody()), DataDelegateInner(inner_relation),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      pos_(particles_->getVariableDataByName<Vecd>("Position")),
      B_(particles_->getVariableDataByName<Matd>("LinearGradientCorrectionMatrix")),
      F_(particles_->registerStateVariable<Matd>("DeformationGradient", IdentityMatrix<Matd>::value)) {}
//=================================================================================================//
BaseElasticIntegration::
    BaseElasticIntegration(BaseInnerRelation &inner_relation)
    : LocalDynamics(inner_relation.getSPHBody()), DataDelegateInner(inner_relation),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      pos_(particles_->getVariableDataByName<Vecd>("Position")),
      vel_(particles_->registerStateVariable<Vecd>("Velocity")),
      force_(particles_->registerStateVariable<Vecd>("Force")),
      B_(particles_->getVariableDataByName<Matd>("LinearGradientCorrectionMatrix")),
      F_(particles_->registerStateVariable<Matd>("DeformationGradient", IdentityMatrix<Matd>::value)),
      dF_dt_(particles_->registerStateVariable<Matd>("DeformationRate")) {}
//=================================================================================================//
BaseIntegration1stHalf::
    BaseIntegration1stHalf(BaseInnerRelation &inner_relation)
    : BaseElasticIntegration(inner_relation),
      elastic_solid_(DynamicCast<ElasticSolid>(this, sph_body_.getBaseMaterial())),
      rho0_(elastic_solid_.ReferenceDensity()), inv_rho0_(1.0 / rho0_),
      rho_(particles_->getVariableDataByName<Real>("Density")),
      mass_(particles_->getVariableDataByName<Real>("Mass")),
      force_prior_(particles_->registerStateVariable<Vecd>("ForcePrior")),
      smoothing_length_(sph_body_.getSPHAdaptation().ReferenceSmoothingLength()) {}
//=================================================================================================//
void BaseIntegration1stHalf::update(size_t index_i, Real dt)
{
    vel_[index_i] += (force_prior_[index_i] + force_[index_i]) / mass_[index_i] * dt;
}
//=================================================================================================//
Integration1stHalf::Integration1stHalf(BaseInnerRelation &inner_relation)
    : BaseIntegration1stHalf(inner_relation),
      stress_PK1_B_(particles_->registerStateVariable<Matd>("StressPK1OnParticle")),
      numerical_dissipation_factor_(0.25) {}
//=================================================================================================//
Integration1stHalfPK2::Integration1stHalfPK2(BaseInnerRelation &inner_relation)
    : Integration1stHalf(inner_relation) {};
//=================================================================================================//
void Integration1stHalfPK2::initialization(size_t index_i, Real dt)
{
    pos_[index_i] += vel_[index_i] * dt * 0.5;
    F_[index_i] += dF_dt_[index_i] * dt * 0.5;
    rho_[index_i] = rho0_ / F_[index_i].determinant();
    // obtain the first Piola-Kirchhoff stress from the second Piola-Kirchhoff stress
    // it seems using reproducing correction here increases convergence rate near the free surface, note that the correction matrix is in a form of transpose
    stress_PK1_B_[index_i] = elastic_solid_.StressPK1(F_[index_i], index_i) * B_[index_i].transpose();
}
//=================================================================================================//
Integration1stHalfKirchhoff::
    Integration1stHalfKirchhoff(BaseInnerRelation &inner_relation)
    : Integration1stHalf(inner_relation) {};
//=================================================================================================//
void Integration1stHalfKirchhoff::initialization(size_t index_i, Real dt)
{
    pos_[index_i] += vel_[index_i] * dt * 0.5;
    F_[index_i] += dF_dt_[index_i] * dt * 0.5;
    rho_[index_i] = rho0_ / F_[index_i].determinant();
    Real J = F_[index_i].determinant();
    Real one_over_J = 1.0 / J;
    rho_[index_i] = rho0_ * one_over_J;
    Real J_to_minus_2_over_dimension = pow(one_over_J, 2.0 * OneOverDimensions);
    Matd normalized_b = (F_[index_i] * F_[index_i].transpose()) * J_to_minus_2_over_dimension;
    Matd deviatoric_b = normalized_b - Matd::Identity() * normalized_b.trace() * OneOverDimensions;
    Matd inverse_F_T = F_[index_i].inverse().transpose();
    // obtain the first Piola-Kirchhoff stress from the Kirchhoff stress
    // it seems using reproducing correction here increases convergence rate
    // near the free surface however, this correction is not used for the numerical dissipation
    stress_PK1_B_[index_i] = (Matd::Identity() * elastic_solid_.VolumetricKirchhoff(J) +
                              elastic_solid_.DeviatoricKirchhoff(deviatoric_b)) *
                             inverse_F_T * B_[index_i];
}
//=================================================================================================//
Integration1stHalfCauchy::
    Integration1stHalfCauchy(BaseInnerRelation &inner_relation)
    : Integration1stHalf(inner_relation) {}
//=================================================================================================//
void Integration1stHalfCauchy::initialization(size_t index_i, Real dt)
{
    pos_[index_i] += vel_[index_i] * dt * 0.5;
    F_[index_i] += dF_dt_[index_i] * dt * 0.5;
    Real J = F_[index_i].determinant();
    Matd F_T = F_[index_i].transpose();
    rho_[index_i] = rho0_ / J;
    Matd inverse_F_T = F_T.inverse();
    Matd almansi_strain = 0.5 * (Matd::Identity() - (F_[index_i] * F_T).inverse());
    // obtain the first Piola-Kirchhoff stress from the  Cauchy stress
    stress_PK1_B_[index_i] = J * elastic_solid_.StressCauchy(almansi_strain, index_i) *
                             inverse_F_T * B_[index_i];
}
//=================================================================================================//
DecomposedIntegration1stHalf::
    DecomposedIntegration1stHalf(BaseInnerRelation &inner_relation)
    : BaseIntegration1stHalf(inner_relation),
      J_to_minus_2_over_dimension_(particles_->registerStateVariable<Real>("DeterminantTerm")),
      stress_on_particle_(particles_->registerStateVariable<Matd>("StressOnParticle")),
      inverse_F_T_(particles_->registerStateVariable<Matd>("InverseTransposedDeformation")) {}
//=================================================================================================//
void DecomposedIntegration1stHalf::initialization(size_t index_i, Real dt)
{
    pos_[index_i] += vel_[index_i] * dt * 0.5;
    F_[index_i] += dF_dt_[index_i] * dt * 0.5;
    Real J = F_[index_i].determinant();
    Real one_over_J = 1.0 / J;
    rho_[index_i] = rho0_ * one_over_J;
    J_to_minus_2_over_dimension_[index_i] = pow(one_over_J * one_over_J, OneOverDimensions);

    inverse_F_T_[index_i] = F_[index_i].inverse().transpose();
    stress_on_particle_[index_i] =
        inverse_F_T_[index_i] * (elastic_solid_.VolumetricKirchhoff(J) -
                                 correction_factor_ * elastic_solid_.ShearModulus() * J_to_minus_2_over_dimension_[index_i] *
                                     (F_[index_i] * F_[index_i].transpose()).trace() * OneOverDimensions) +
        elastic_solid_.NumericalDampingLeftCauchy(F_[index_i], dF_dt_[index_i], smoothing_length_, index_i) * inverse_F_T_[index_i];
}
//=================================================================================================//
void Integration2ndHalf::initialization(size_t index_i, Real dt)
{
    pos_[index_i] += vel_[index_i] * dt * 0.5;
}
//=================================================================================================//
void Integration2ndHalf::update(size_t index_i, Real dt)
{
    F_[index_i] += dF_dt_[index_i] * dt * 0.5;
}
//=================================================================================================//
void Integration1stHalfPK2RightCauchy::initialization(size_t index_i, Real dt)
{
    Integration1stHalfPK2::initialization(index_i, dt);
    // add damping stress
    const Matd numerical_damping_stress = elastic_solid_.NumericalDampingRightCauchy(F_[index_i], dF_dt_[index_i], smoothing_length_ / h_ratio_[index_i], index_i);
    stress_PK1_B_[index_i] += F_[index_i] * 0.5 * numerical_dissipation_factor_ * numerical_damping_stress;
}
//=================================================================================================//
} // namespace solid_dynamics
} // namespace SPH
