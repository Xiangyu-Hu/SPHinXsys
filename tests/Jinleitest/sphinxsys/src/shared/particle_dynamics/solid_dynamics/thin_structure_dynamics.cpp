#include "thin_structure_dynamics.h"
#include "base_particles.hpp"

namespace SPH
{
//=====================================================================================================//
namespace thin_structure_dynamics
{
//=================================================================================================//
UpdateShellNormalDirection::UpdateShellNormalDirection(SPHBody &sph_body)
    : LocalDynamics(sph_body),
      n_(particles_->getVariableDataByName<Vecd>("NormalDirection")),
      F_(particles_->getVariableDataByName<Matd>("DeformationGradient")),
      transformation_matrix0_(particles_->getVariableDataByName<Matd>("TransformationMatrix")) {}
//=========================================================================================//
void UpdateShellNormalDirection::update(size_t index_i, Real dt)
{
    /** Calculate the current normal direction of mid-surface. */
    n_[index_i] = transformation_matrix0_[index_i].transpose() * getNormalFromDeformationGradientTensor(F_[index_i]);
}
//=================================================================================================//
ShellAcousticTimeStepSize::ShellAcousticTimeStepSize(SPHBody &sph_body, Real CFL)
    : LocalDynamicsReduce<ReduceMin>(sph_body),
      CFL_(CFL),
      elastic_solid_(DynamicCast<ElasticSolid>(this, sph_body.getBaseMaterial())),
      vel_(particles_->getVariableDataByName<Vecd>("Velocity")),
      force_(particles_->getVariableDataByName<Vecd>("Force")),
      angular_vel_(particles_->getVariableDataByName<Vecd>("AngularVelocity")),
      dangular_vel_dt_(particles_->getVariableDataByName<Vecd>("AngularAcceleration")),
      force_prior_(particles_->getVariableDataByName<Vecd>("ForcePrior")),
      thickness_(particles_->getVariableDataByName<Real>("Thickness")),
      mass_(particles_->getVariableDataByName<Real>("Mass")),
      rho0_(elastic_solid_.ReferenceDensity()),
      E0_(elastic_solid_.YoungsModulus()),
      nu_(elastic_solid_.PoissonRatio()),
      c0_(elastic_solid_.ReferenceSoundSpeed()),
      smoothing_length_(sph_body.getSPHAdaptation().ReferenceSmoothingLength()) {}
//=================================================================================================//
Real ShellAcousticTimeStepSize::reduce(size_t index_i, Real dt)
{
    // Since the particle does not change its configuration in pressure relaxation step,
    // I chose a time-step size according to Eulerian method.
    Real time_setp_0 = SMIN((Real)sqrt(smoothing_length_ / ((force_[index_i] + force_prior_[index_i]).norm() / mass_[index_i] + TinyReal)),
                            smoothing_length_ / (c0_ + vel_[index_i].norm()));
    Real time_setp_1 = SMIN((Real)sqrt(1.0 / (dangular_vel_dt_[index_i].norm() + TinyReal)),
                            Real(1.0) / (angular_vel_[index_i].norm() + TinyReal));
    Real time_setp_2 = smoothing_length_ * (Real)sqrt(rho0_ * (1.0 - nu_ * nu_) / E0_ /
                                                      (2.0 + (Pi * Pi / 12.0) * (1.0 - nu_) *
                                                                 (1.0 + 1.5 * pow(smoothing_length_ / thickness_[index_i], 2))));
    return CFL_ * SMIN(time_setp_0, time_setp_1, time_setp_2);
}
//=================================================================================================//
ShellCorrectConfiguration::
    ShellCorrectConfiguration(BaseInnerRelation &inner_relation)
    : LocalDynamics(inner_relation.getSPHBody()), DataDelegateInner(inner_relation),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      B_(particles_->registerStateVariable<Matd>("LinearGradientCorrectionMatrix", IdentityMatrix<Matd>::value)),
      n0_(particles_->registerStateVariableFrom<Vecd>("InitialNormalDirection", "NormalDirection")),
      transformation_matrix0_(particles_->getVariableDataByName<Matd>("TransformationMatrix")) {}
//=================================================================================================//
ShellDeformationGradientTensor::
    ShellDeformationGradientTensor(BaseInnerRelation &inner_relation)
    : LocalDynamics(inner_relation.getSPHBody()), DataDelegateInner(inner_relation),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      pos_(particles_->getVariableDataByName<Vecd>("Position")),
      pseudo_n_(particles_->registerStateVariableFrom<Vecd>("PseudoNormal", "NormalDirection")),
      n0_(particles_->registerStateVariableFrom<Vecd>("InitialNormalDirection", "NormalDirection")),
      B_(particles_->getVariableDataByName<Matd>("LinearGradientCorrectionMatrix")),
      F_(particles_->registerStateVariable<Matd>("DeformationGradient", IdentityMatrix<Matd>::value)),
      F_bending_(particles_->registerStateVariable<Matd>("BendingDeformationGradient")),
      transformation_matrix0_(particles_->getVariableDataByName<Matd>("TransformationMatrix")) {}
//=================================================================================================//
BaseShellRelaxation::BaseShellRelaxation(BaseInnerRelation &inner_relation)
    : LocalDynamics(inner_relation.getSPHBody()), DataDelegateInner(inner_relation),
      thickness_(particles_->getVariableDataByName<Real>("Thickness")),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      pos_(particles_->getVariableDataByName<Vecd>("Position")),
      vel_(particles_->registerStateVariable<Vecd>("Velocity")),
      force_(particles_->registerStateVariable<Vecd>("Force")),
      force_prior_(particles_->registerStateVariable<Vecd>("ForcePrior")),
      n0_(particles_->registerStateVariableFrom<Vecd>("InitialNormalDirection", "NormalDirection")),
      pseudo_n_(particles_->registerStateVariableFrom<Vecd>("PseudoNormal", "NormalDirection")),
      dpseudo_n_dt_(particles_->registerStateVariable<Vecd>("PseudoNormalChangeRate")),
      dpseudo_n_d2t_(particles_->registerStateVariable<Vecd>("PseudoNormal2ndOrderTimeDerivative")),
      rotation_(particles_->registerStateVariable<Vecd>("Rotation")),
      angular_vel_(particles_->registerStateVariable<Vecd>("AngularVelocity")),
      dangular_vel_dt_(particles_->registerStateVariable<Vecd>("AngularAcceleration")),
      transformation_matrix0_(particles_->getVariableDataByName<Matd>("TransformationMatrix")),
      B_(particles_->getVariableDataByName<Matd>("LinearGradientCorrectionMatrix")),
      F_(particles_->registerStateVariable<Matd>("DeformationGradient", IdentityMatrix<Matd>::value)),
      dF_dt_(particles_->registerStateVariable<Matd>("DeformationRate")),
      F_bending_(particles_->registerStateVariable<Matd>("BendingDeformationGradient")),
      dF_bending_dt_(particles_->registerStateVariable<Matd>("BendingDeformationRate")) {}
//=================================================================================================//
ShellStressRelaxationFirstHalf::
    ShellStressRelaxationFirstHalf(BaseInnerRelation &inner_relation,
                                   int number_of_gaussian_points, bool hourglass_control,
                                   Real hourglass_control_factor)
    : BaseShellRelaxation(inner_relation),
      elastic_solid_(DynamicCast<ElasticSolid>(this, sph_body_.getBaseMaterial())),
      rho0_(elastic_solid_.ReferenceDensity()),
      inv_rho0_(1.0 / rho0_),
      smoothing_length_(sph_body_.getSPHAdaptation().ReferenceSmoothingLength()),
      numerical_damping_scaling_matrix_(Matd::Identity() * smoothing_length_),
      rho_(particles_->getVariableDataByName<Real>("Density")),
      mass_(particles_->getVariableDataByName<Real>("Mass")),
      global_stress_(particles_->registerStateVariable<Matd>("GlobalStress")),
      global_moment_(particles_->registerStateVariable<Matd>("GlobalMoment")),
      mid_surface_cauchy_stress_(particles_->registerStateVariable<Matd>("MidSurfaceCauchyStress")),
      global_shear_stress_(particles_->registerStateVariable<Vecd>("GlobalShearStress")),
      global_F_(particles_->registerStateVariable<Matd>("GlobalDeformationGradient")),
      global_F_bending_(particles_->registerStateVariable<Matd>("GlobalBendingDeformationGradient")),
      E0_(elastic_solid_.YoungsModulus()),
      G0_(elastic_solid_.ShearModulus()),
      nu_(elastic_solid_.PoissonRatio()),
      hourglass_control_factor_(hourglass_control_factor),
      hourglass_control_(hourglass_control),
      number_of_gaussian_points_(number_of_gaussian_points)
{
    /** Note that, only three-point and five-point Gaussian quadrature rules are defined. */
    switch (number_of_gaussian_points)
    {
    case 1:
        gaussian_point_ = one_gaussian_point_;
        gaussian_weight_ = one_gaussian_weight_;
        break;
    case 5:
        gaussian_point_ = five_gaussian_points_;
        gaussian_weight_ = five_gaussian_weights_;
        break;
    default:
        gaussian_point_ = three_gaussian_points_;
        gaussian_weight_ = three_gaussian_weights_;
    }
}
//=================================================================================================//
void ShellStressRelaxationFirstHalf::initialization(size_t index_i, Real dt)
{
    // Note that F_[index_i], F_bending_[index_i], dF_dt_[index_i], dF_bending_dt_[index_i]
    // and rotation_[index_i], angular_vel_[index_i], dangular_vel_dt_[index_i], B_[index_i]
    // are defined in local coordinates, while others in global coordinates.
    pos_[index_i] += vel_[index_i] * dt * 0.5;
    rotation_[index_i] += angular_vel_[index_i] * dt * 0.5;
    pseudo_n_[index_i] += dpseudo_n_dt_[index_i] * dt * 0.5;

    F_[index_i] += dF_dt_[index_i] * dt * 0.5;
    F_bending_[index_i] += dF_bending_dt_[index_i] * dt * 0.5;

    global_F_[index_i] = transformation_matrix0_[index_i].transpose() * F_[index_i] * transformation_matrix0_[index_i];
    global_F_bending_[index_i] = transformation_matrix0_[index_i].transpose() * F_bending_[index_i] * transformation_matrix0_[index_i];

    Real J = F_[index_i].determinant();
    Matd inverse_transpose_global_F = global_F_[index_i].inverse().transpose();

    rho_[index_i] = rho0_ / J;

    /** Get transformation matrix from global coordinates to current local coordinates. */
    Matd current_transformation_matrix = getTransformationMatrix(pseudo_n_[index_i]);
    /** Get transformation matrix from initial local coordinates to current local coordinates. */
    Matd transformation_matrix_0_to_current = current_transformation_matrix * transformation_matrix0_[index_i].transpose();

    Matd resultant_stress = Matd::Zero();
    Matd resultant_moment = Matd::Zero();
    Vecd resultant_shear_stress = Vecd::Zero();

    for (int i = 0; i != number_of_gaussian_points_; ++i)
    {
        Matd F_gaussian_point = F_[index_i] + gaussian_point_[i] * F_bending_[index_i] * thickness_[index_i] * 0.5;
        Matd dF_gaussian_point_dt = dF_dt_[index_i] + gaussian_point_[i] * dF_bending_dt_[index_i] * thickness_[index_i] * 0.5;
        Matd inverse_F_gaussian_point = F_gaussian_point.inverse();
        Matd current_local_almansi_strain = transformation_matrix_0_to_current * 0.5 *
                                            (Matd::Identity() - inverse_F_gaussian_point.transpose() * inverse_F_gaussian_point) *
                                            transformation_matrix_0_to_current.transpose();

        /** correct Almansi strain tensor according to plane stress problem. */
        current_local_almansi_strain = getCorrectedAlmansiStrain(current_local_almansi_strain, nu_);

        /** correct out-plane numerical damping. */
        numerical_damping_scaling_matrix_(Dimensions - 1, Dimensions - 1) = thickness_[i] < smoothing_length_ ? thickness_[i] : smoothing_length_;
        Matd cauchy_stress = elastic_solid_.StressCauchy(current_local_almansi_strain, index_i) +
                             transformation_matrix_0_to_current * F_gaussian_point *
                                 elastic_solid_.NumericalDampingRightCauchy(F_gaussian_point, dF_gaussian_point_dt, numerical_damping_scaling_matrix_, index_i) *
                                 F_gaussian_point.transpose() * transformation_matrix_0_to_current.transpose() / F_gaussian_point.determinant();

        /** Impose modeling assumptions. */
        cauchy_stress.col(Dimensions - 1) *= shear_correction_factor_;
        cauchy_stress.row(Dimensions - 1) *= shear_correction_factor_;
        cauchy_stress(Dimensions - 1, Dimensions - 1) = 0.0;

        if (i == 0)
        {
            mid_surface_cauchy_stress_[index_i] = cauchy_stress;
        }

        /** Integrate Cauchy stress along thickness. */
        resultant_stress +=
            0.5 * thickness_[index_i] * gaussian_weight_[i] * cauchy_stress;
        resultant_moment +=
            0.5 * thickness_[index_i] * gaussian_weight_[i] * (cauchy_stress * gaussian_point_[i] * thickness_[index_i] * 0.5);
        resultant_shear_stress -=
            0.5 * thickness_[index_i] * gaussian_weight_[i] * cauchy_stress.col(Dimensions - 1);

        resultant_stress.col(Dimensions - 1) = Vecd::Zero();
        resultant_moment.col(Dimensions - 1) = Vecd::Zero();
    }

    /** stress and moment in global coordinates for pair interaction */
    global_stress_[index_i] = J * current_transformation_matrix.transpose() *
                              resultant_stress * current_transformation_matrix * inverse_transpose_global_F;
    global_moment_[index_i] = J * current_transformation_matrix.transpose() *
                              resultant_moment * current_transformation_matrix * inverse_transpose_global_F;
    global_shear_stress_[index_i] = J * current_transformation_matrix.transpose() * resultant_shear_stress;
}
//=================================================================================================//
void ShellStressRelaxationFirstHalf::update(size_t index_i, Real dt)
{
    vel_[index_i] += (force_prior_[index_i] + force_[index_i]) / mass_[index_i] * dt;
    angular_vel_[index_i] += dangular_vel_dt_[index_i] * dt;
}
//=================================================================================================//
void ShellStressRelaxationSecondHalf::initialization(size_t index_i, Real dt)
{
    pos_[index_i] += vel_[index_i] * dt * 0.5;
    rotation_[index_i] += angular_vel_[index_i] * dt * 0.5;
    dpseudo_n_dt_[index_i] = transformation_matrix0_[index_i].transpose() *
                             getVectorChangeRateAfterThinStructureRotation(local_pseudo_n_0, rotation_[index_i], angular_vel_[index_i]);
    pseudo_n_[index_i] += dpseudo_n_dt_[index_i] * dt * 0.5;
}
//=================================================================================================//
void ShellStressRelaxationSecondHalf::update(size_t index_i, Real dt)
{
    F_[index_i] += dF_dt_[index_i] * dt * 0.5;
    F_bending_[index_i] += dF_bending_dt_[index_i] * dt * 0.5;
}
//=================================================================================================//
ConstrainShellBodyRegion::
    ConstrainShellBodyRegion(BodyPartByParticle &body_part)
    : BaseLocalDynamics<BodyPartByParticle>(body_part),
      vel_(particles_->getVariableDataByName<Vecd>("Velocity")),
      angular_vel_(particles_->getVariableDataByName<Vecd>("AngularVelocity")) {}
//=================================================================================================//
void ConstrainShellBodyRegion::update(size_t index_i, Real dt)
{
    vel_[index_i] = Vecd::Zero();
    angular_vel_[index_i] = Vecd::Zero();
}
//=================================================================================================//
ConstrainShellBodyRegionAlongAxis::ConstrainShellBodyRegionAlongAxis(BodyPartByParticle &body_part, int axis)
    : BaseLocalDynamics<BodyPartByParticle>(body_part),
      axis_(axis), pos_(particles_->getVariableDataByName<Vecd>("Position")),
      pos0_(particles_->registerStateVariableFrom<Vecd>("InitialPosition", "Position")),
      vel_(particles_->getVariableDataByName<Vecd>("Velocity")),
      force_(particles_->getVariableDataByName<Vecd>("Force")),
      rotation_(particles_->getVariableDataByName<Vecd>("Rotation")),
      angular_vel_(particles_->getVariableDataByName<Vecd>("AngularVelocity")),
      dangular_vel_dt_(particles_->getVariableDataByName<Vecd>("AngularAcceleration")),
      mass_(particles_->getVariableDataByName<Real>("Mass")) {}
//=================================================================================================//
void ConstrainShellBodyRegionAlongAxis::update(size_t index_i, Real dt)
{
    vel_[index_i][axis_] = 0.0;
    vel_[index_i][2] = 0.0;
    force_[index_i][axis_] = 0.0;
    force_[index_i][2] = 0.0;

    angular_vel_[index_i][1 - axis_] = 0.0;
    dangular_vel_dt_[index_i][1 - axis_] = 0.0;
}
//=================================================================================================//
InitialShellCurvature::InitialShellCurvature(BaseInnerRelation &inner_relation)
    : LocalDynamics(inner_relation.getSPHBody()), DataDelegateInner(inner_relation),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      n0_(particles_->registerStateVariableFrom<Vecd>("InitialNormalDirection", "NormalDirection")),
      B_(particles_->getVariableDataByName<Matd>("LinearGradientCorrectionMatrix")),
      transformation_matrix0_(particles_->getVariableDataByName<Matd>("TransformationMatrix")),
      n_(particles_->getVariableDataByName<Vecd>("NormalDirection")),
      F_(particles_->getVariableDataByName<Matd>("DeformationGradient")),
      F_bending_(particles_->getVariableDataByName<Matd>("BendingDeformationGradient")),
      k1_(particles_->registerStateVariable<Real>("1stPrincipleCurvature")),
      k2_(particles_->registerStateVariable<Real>("2ndPrincipleCurvature")),
      dn_0_(particles_->registerStateVariable<Matd>("InitialNormalGradient")) {};
//=================================================================================================//
void InitialShellCurvature::update(size_t index_i, Real)
{

    Matd dn_0_i = Matd::Zero();
    // transform initial local B_ to global B_
    const Matd B_global_i = transformation_matrix0_[index_i].transpose() * B_[index_i] * transformation_matrix0_[index_i];
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        const size_t index_j = inner_neighborhood.j_[n];
        const Vecd gradW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j] * inner_neighborhood.e_ij_[n];
        dn_0_i -= (n0_[index_i] - n0_[index_j]) * gradW_ijV_j.transpose();
    }
    dn_0_[index_i] = dn_0_i * B_global_i;
    auto [k1, k2] = get_principle_curvatures(dn_0_[index_i]);
    k1_[index_i] = k1;
    k2_[index_i] = k2;
}
//=================================================================================================//
ShellCurvatureUpdate::ShellCurvatureUpdate(SPHBody &sph_body)
    : LocalDynamics(sph_body),
      transformation_matrix0_(particles_->getVariableDataByName<Matd>("TransformationMatrix")),
      F_(particles_->getVariableDataByName<Matd>("DeformationGradient")),
      F_bending_(particles_->getVariableDataByName<Matd>("BendingDeformationGradient")),
      k1_(particles_->registerStateVariable<Real>("1stPrincipleCurvature")),
      k2_(particles_->registerStateVariable<Real>("2ndPrincipleCurvature")),
      dn_0_(particles_->registerStateVariable<Matd>("InitialNormalGradient")) {};
//=================================================================================================//
void ShellCurvatureUpdate::update(size_t index_i, Real)
{
    Matd dn_0_i = dn_0_[index_i] + transformation_matrix0_[index_i].transpose() *
                                       F_bending_[index_i] * transformation_matrix0_[index_i];
    Matd inverse_F = F_[index_i].inverse();
    Matd dn_i = dn_0_i * transformation_matrix0_[index_i].transpose() * inverse_F * transformation_matrix0_[index_i];
    auto [k1, k2] = get_principle_curvatures(dn_i);
    k1_[index_i] = k1;
    k2_[index_i] = k2;
}
//=================================================================================================//
AverageShellCurvature::AverageShellCurvature(BaseInnerRelation &inner_relation)
    : LocalDynamics(inner_relation.getSPHBody()), DataDelegateInner(inner_relation),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      n_(particles_->getVariableDataByName<Vecd>("NormalDirection")),
      k1_ave_(particles_->registerStateVariable<Real>("Average1stPrincipleCurvature")),
      k2_ave_(particles_->registerStateVariable<Real>("Average2ndPrincipleCurvature")) {};
//=================================================================================================//
void AverageShellCurvature::update(size_t index_i, Real)
{
    Matd dn_i = Matd::Zero();
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        const size_t index_j = inner_neighborhood.j_[n];
        const Vecd gradW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j] * inner_neighborhood.e_ij_[n];
        dn_i -= (n_[index_i] - n_[index_j]) * gradW_ijV_j.transpose();
    }
    auto [k1, k2] = get_principle_curvatures(dn_i);
    k1_ave_[index_i] = k1;
    k2_ave_[index_i] = k2;
}
//=================================================================================================//
} // namespace thin_structure_dynamics
} // namespace SPH