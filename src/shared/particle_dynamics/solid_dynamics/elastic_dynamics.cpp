#include "elastic_dynamics.h"
#include "general_dynamics.h"

#include <numeric>

namespace SPH
{
//=========================================================================================================//
namespace solid_dynamics
{
//=================================================================================================//
AcousticTimeStepSize::AcousticTimeStepSize(SPHBody &sph_body, Real CFL)
    : LocalDynamicsReduce<Real, ReduceMin>(sph_body, Real(MaxRealNumber)),
      ElasticSolidDataSimple(sph_body), CFL_(CFL),
      vel_(particles_->vel_), acc_(particles_->acc_), acc_prior_(particles_->acc_prior_),
      smoothing_length_(sph_body.sph_adaptation_->ReferenceSmoothingLength()),
      c0_(particles_->elastic_solid_.ReferenceSoundSpeed()) {}
//=================================================================================================//
Real AcousticTimeStepSize::reduce(size_t index_i, Real dt)
{
    // since the particle does not change its configuration in pressure relaxation step
    // I chose a time-step size according to Eulerian method
    return CFL_ * SMIN((Real)sqrt(smoothing_length_ / ((acc_[index_i] + acc_prior_[index_i]).norm() + TinyReal)),
                       smoothing_length_ / (c0_ + vel_[index_i].norm()));
}
//=================================================================================================//
ElasticDynamicsInitialCondition::ElasticDynamicsInitialCondition(SPHBody &sph_body)
    : LocalDynamics(sph_body),
      ElasticSolidDataSimple(sph_body),
      pos_(particles_->pos_), vel_(particles_->vel_) {}
//=================================================================================================//
UpdateElasticNormalDirection::UpdateElasticNormalDirection(SPHBody &sph_body)
    : LocalDynamics(sph_body),
      ElasticSolidDataSimple(sph_body),
      n_(particles_->n_), n0_(particles_->n0_), F_(particles_->F_) {}
//=================================================================================================//
DeformationGradientBySummation::
    DeformationGradientBySummation(BaseInnerRelation &inner_relation)
    : LocalDynamics(inner_relation.getSPHBody()), ElasticSolidDataInner(inner_relation),
      pos_(particles_->pos_), B_(particles_->B_), F_(particles_->F_) {}
//=================================================================================================//
BaseElasticIntegration::
    BaseElasticIntegration(BaseInnerRelation &inner_relation)
    : LocalDynamics(inner_relation.getSPHBody()), ElasticSolidDataInner(inner_relation),
      rho_(particles_->rho_), mass_(particles_->mass_),
      pos_(particles_->pos_), vel_(particles_->vel_), acc_(particles_->acc_),
      B_(particles_->B_), F_(particles_->F_), dF_dt_(particles_->dF_dt_) {}
//=================================================================================================//
BaseIntegration1stHalf::
    BaseIntegration1stHalf(BaseInnerRelation &inner_relation)
    : BaseElasticIntegration(inner_relation),
      elastic_solid_(particles_->elastic_solid_),
      acc_prior_(particles_->acc_prior_)
{
    rho0_ = particles_->elastic_solid_.ReferenceDensity();
    inv_rho0_ = 1.0 / rho0_;
    smoothing_length_ = sph_body_.sph_adaptation_->ReferenceSmoothingLength();
}
//=================================================================================================//
void BaseIntegration1stHalf::update(size_t index_i, Real dt)
{
    vel_[index_i] += (acc_prior_[index_i] + acc_[index_i]) * dt;
}
//=================================================================================================//
Integration1stHalf::
    Integration1stHalf(BaseInnerRelation &inner_relation)
    : BaseIntegration1stHalf(inner_relation)
{
    particles_->registerVariable(stress_PK1_B_, "CorrectedStressPK1");
    numerical_dissipation_factor_ = 0.25;
}
//=================================================================================================//
Integration1stHalfPK2::Integration1stHalfPK2(BaseInnerRelation &inner_relation)
    : Integration1stHalf(inner_relation){};
//=================================================================================================//
void Integration1stHalfPK2::initialization(size_t index_i, Real dt)
{
    pos_[index_i] += vel_[index_i] * dt * 0.5;
    F_[index_i] += dF_dt_[index_i] * dt * 0.5;
    rho_[index_i] = rho0_ / F_[index_i].determinant();
    // obtain the first Piola-Kirchhoff stress from the second Piola-Kirchhoff stress
    // it seems using reproducing correction here increases convergence rate near the free surface
    stress_PK1_B_[index_i] = elastic_solid_.StressPK1(F_[index_i], index_i) * B_[index_i];
}
//=================================================================================================//
Integration1stHalfKirchhoff::
    Integration1stHalfKirchhoff(BaseInnerRelation &inner_relation)
    : Integration1stHalf(inner_relation){};
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
    stress_PK1_B_[index_i] = J * elastic_solid_.StressCauchy(almansi_strain, F_[index_i], index_i) *
                             inverse_F_T * B_[index_i];
}
//=================================================================================================//
DecomposedIntegration1stHalf::
    DecomposedIntegration1stHalf(BaseInnerRelation &inner_relation)
    : BaseIntegration1stHalf(inner_relation)
{
    particles_->registerVariable(J_to_minus_2_over_dimension_, "DeterminantTerm");
    particles_->registerVariable(stress_on_particle_, "StressOnParticle");
    particles_->registerVariable(inverse_F_T_, "InverseTransposedDeformation");
};
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
} // namespace solid_dynamics
} // namespace SPH
