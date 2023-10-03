#include "relax_dynamics.h"
#include "base_particles.hpp"
#include "level_set_shape.h"
#include "particle_generator_lattice.h"

//========================================================================================================//
namespace SPH
{
//=====================================================================================================//
namespace relax_dynamics
{
//=================================================================================================//
GetTimeStepSizeSquare::GetTimeStepSizeSquare(SPHBody &sph_body)
    : LocalDynamicsReduce<Real, ReduceMax>(sph_body, Real(0)),
      RelaxDataDelegateSimple(sph_body), acc_(particles_->acc_),
      h_ref_(sph_body.sph_adaptation_->ReferenceSmoothingLength()) {}
//=================================================================================================//
Real GetTimeStepSizeSquare::reduce(size_t index_i, Real dt)
{
    return acc_[index_i].norm();
}
//=================================================================================================//
Real GetTimeStepSizeSquare::outputResult(Real reduced_value)
{
    return 0.0625 * h_ref_ / (reduced_value + TinyReal);
}
//=================================================================================================//
RelaxationAccelerationInner::RelaxationAccelerationInner(BaseInnerRelation &inner_relation)
    : LocalDynamics(inner_relation.getSPHBody()), RelaxDataDelegateInner(inner_relation),
      acc_(particles_->acc_), pos_(particles_->pos_) {}
//=================================================================================================//
RelaxationAccelerationInnerWithLevelSetCorrection::
    RelaxationAccelerationInnerWithLevelSetCorrection(BaseInnerRelation &inner_relation)
    : RelaxationAccelerationInner(inner_relation), sph_adaptation_(sph_body_.sph_adaptation_)
{
    level_set_shape_ = DynamicCast<LevelSetShape>(this, sph_body_.body_shape_);
}
//=================================================================================================//
UpdateParticlePosition::UpdateParticlePosition(SPHBody &sph_body)
    : LocalDynamics(sph_body), RelaxDataDelegateSimple(sph_body),
      sph_adaptation_(sph_body.sph_adaptation_),
      pos_(particles_->pos_), acc_(particles_->acc_) {}
//=================================================================================================//
void UpdateParticlePosition::update(size_t index_i, Real dt_square)
{
    pos_[index_i] += acc_[index_i] * dt_square * 0.5 / sph_adaptation_->SmoothingLengthRatio(index_i);
}
//=================================================================================================//
UpdateSmoothingLengthRatioByShape::
    UpdateSmoothingLengthRatioByShape(SPHBody &sph_body, Shape &target_shape)
    : LocalDynamics(sph_body), RelaxDataDelegateSimple(sph_body),
      h_ratio_(*particles_->getVariableByName<Real>("SmoothingLengthRatio")),
      Vol_(particles_->Vol_), pos_(particles_->pos_), target_shape_(target_shape),
      particle_adaptation_(DynamicCast<ParticleRefinementByShape>(this, sph_body.sph_adaptation_)),
      reference_spacing_(particle_adaptation_->ReferenceSpacing()) {}
//=================================================================================================//
UpdateSmoothingLengthRatioByShape::UpdateSmoothingLengthRatioByShape(SPHBody &sph_body)
    : UpdateSmoothingLengthRatioByShape(sph_body, *sph_body.body_shape_) {}
//=================================================================================================//
void UpdateSmoothingLengthRatioByShape::update(size_t index_i, Real dt_square)
{
    Real local_spacing = particle_adaptation_->getLocalSpacing(target_shape_, pos_[index_i]);
    h_ratio_[index_i] = reference_spacing_ / local_spacing;
    Vol_[index_i] = pow(local_spacing, Dimensions);
}
//=================================================================================================//
RelaxationAccelerationComplex::
    RelaxationAccelerationComplex(ComplexRelation &complex_relation)
    : LocalDynamics(complex_relation.getSPHBody()),
      RelaxDataDelegateComplex(complex_relation),
      acc_(particles_->acc_), pos_(particles_->pos_) {}
//=================================================================================================//
ShapeSurfaceBounding::ShapeSurfaceBounding(NearShapeSurface &near_shape_surface)
    : BaseLocalDynamics<BodyPartByCell>(near_shape_surface),
      RelaxDataDelegateSimple(sph_body_), pos_(particles_->pos_),
      constrained_distance_(0.5 * sph_body_.sph_adaptation_->MinimumSpacing())
{
    level_set_shape_ = &near_shape_surface.level_set_shape_;
}
//=================================================================================================//
void ShapeSurfaceBounding::update(size_t index_i, Real dt)
{
    Real phi = level_set_shape_->findSignedDistance(pos_[index_i]);

    if (phi > -constrained_distance_)
    {
        Vecd unit_normal = level_set_shape_->findNormalDirection(pos_[index_i]);
        pos_[index_i] -= (phi + constrained_distance_) * unit_normal;
    }
}
//=================================================================================================//
NearSurfaceVolumeCorrection::NearSurfaceVolumeCorrection(NearShapeSurface& near_shape_surface)
    : BaseLocalDynamics<BodyPartByCell>(near_shape_surface),
    RelaxDataDelegateSimple(sph_body_), pos_(particles_->pos_), Vol_(particles_->Vol_),
    level_set_shape_(&near_shape_surface.level_set_shape_),
    sph_adaptation_(sph_body_.sph_adaptation_) {}
//=================================================================================================//
void NearSurfaceVolumeCorrection::update(size_t index_i, Real dt)
{
    Real particle_spacing = sph_adaptation_->getLocalSpacing(*level_set_shape_, pos_[index_i]);
    Real overlap = level_set_shape_->computeKernelIntegral(pos_[index_i], sph_adaptation_->SmoothingLengthRatio(index_i));
    Vol_[index_i] = pow(particle_spacing, Dimensions) * SMAX(0.5, (1.0 - overlap));
}
//=================================================================================================//
RelaxationStepInner::
    RelaxationStepInner(BaseInnerRelation &inner_relation, bool level_set_correction)
    : BaseDynamics<void>(inner_relation.getSPHBody()), real_body_(inner_relation.real_body_),
      inner_relation_(inner_relation), near_shape_surface_(*real_body_),
      get_time_step_square_(*real_body_), update_particle_position_(*real_body_),
      surface_bounding_(near_shape_surface_)
{
    if (!level_set_correction)
    {
        relaxation_acceleration_inner_ =
            makeUnique<InteractionDynamics<RelaxationAccelerationInner>>(inner_relation);
		surface_correction_ = makeShared<SimpleDynamics<ShapeSurfaceBounding>>(near_shape_surface_);
    }
    else
    {
        relaxation_acceleration_inner_ =
            makeUnique<InteractionDynamics<RelaxationAccelerationInnerWithLevelSetCorrection>>(inner_relation);
		surface_correction_ = makeShared<SimpleDynamics<NearSurfaceVolumeCorrection>>(near_shape_surface_);
	}
}
//=================================================================================================//
void RelaxationStepInner::exec(Real dt)
{
    real_body_->updateCellLinkedList();
    inner_relation_.updateConfiguration();
    relaxation_acceleration_inner_->exec();
    Real dt_square = get_time_step_square_.exec();
    update_particle_position_.exec(dt_square);
    surface_bounding_.exec();
}
//=================================================================================================//
RelaxationAccelerationComplexWithLevelSetCorrection::
    RelaxationAccelerationComplexWithLevelSetCorrection(ComplexRelation &complex_relation, const std::string &shape_name)
    : RelaxationAccelerationComplex(complex_relation),
      sph_adaptation_(sph_body_.sph_adaptation_)
{
    ComplexShape &complex_shape = DynamicCast<ComplexShape>(this, *sph_body_.body_shape_);
    level_set_shape_ = DynamicCast<LevelSetShape>(this, complex_shape.getShapeByName(shape_name));
}
//=================================================================================================//
RelaxationStepComplex::RelaxationStepComplex(ComplexRelation &complex_relation,
                                             const std::string &shape_name, bool level_set_correction)
    : BaseDynamics<void>(complex_relation.getSPHBody()),
      real_body_(complex_relation.getInnerRelation().real_body_),
      complex_relation_(complex_relation),
      near_shape_surface_(*real_body_, shape_name),
      get_time_step_square_(*real_body_), update_particle_position_(*real_body_),
      surface_bounding_(near_shape_surface_)
{
    if (!level_set_correction)
    {
        relaxation_acceleration_complex_ =
            makeUnique<InteractionDynamics<RelaxationAccelerationComplex>>(complex_relation);
    }
    else
    {
        relaxation_acceleration_complex_ =
            makeUnique<InteractionDynamics<RelaxationAccelerationComplexWithLevelSetCorrection>>(complex_relation, shape_name);
    }
}
//=================================================================================================//
void RelaxationStepComplex::exec(Real dt)
{
    real_body_->updateCellLinkedList();
    complex_relation_.updateConfiguration();
    relaxation_acceleration_complex_->exec();
    Real dt_square = get_time_step_square_.exec();
    update_particle_position_.exec(dt_square);
    surface_bounding_.exec();
}
//=================================================================================================//
RelaxationAccelerationByCMInner::RelaxationAccelerationByCMInner(BaseInnerRelation& inner_relation)
    : LocalDynamics(inner_relation.getSPHBody()), RelaxDataDelegateInner(inner_relation),
    acc_(particles_->acc_), pos_(particles_->pos_),
    B_(*particles_->template getVariableByName<Matd>("CorrectionMatrix")) {}
//=================================================================================================//
RelaxationAccelerationByCMInnerWithLevelSetCorrection::
RelaxationAccelerationByCMInnerWithLevelSetCorrection(BaseInnerRelation& inner_relation)
    : RelaxationAccelerationByCMInner(inner_relation), sph_adaptation_(sph_body_.sph_adaptation_)
{
    level_set_shape_ = DynamicCast<LevelSetShape>(this, sph_body_.body_shape_);
}
//=================================================================================================//
RelaxationStepByCMInner::RelaxationStepByCMInner(BaseInnerRelation& inner_relation, bool level_set_correction)
    : BaseDynamics<void>(inner_relation.getSPHBody()), real_body_(inner_relation.real_body_),
    inner_relation_(inner_relation),
    near_shape_surface_(*real_body_), get_time_step_square_(*real_body_),
    update_particle_position_(*real_body_), surface_bounding_(near_shape_surface_)
{
    if (!level_set_correction)
    {
        relaxation_acceleration_inner_ =
            makeUnique<InteractionDynamics<RelaxationAccelerationByCMInner>>(inner_relation);
        surface_correction_ = makeShared<SimpleDynamics<ShapeSurfaceBounding>>(near_shape_surface_);
    }
    else
    {
        relaxation_acceleration_inner_ =
            makeUnique<InteractionDynamics<RelaxationAccelerationByCMInnerWithLevelSetCorrection>>(inner_relation);
        surface_correction_ = makeShared<SimpleDynamics<NearSurfaceVolumeCorrection>>(near_shape_surface_);
    }
}
//=================================================================================================//
void RelaxationStepByCMInner::exec(Real dt)
{
    relaxation_acceleration_inner_->exec();
    Real dt_square = get_time_step_square_.exec();
    update_particle_position_.exec(dt_square);
    //surface_correction_->exec();
    //real_body_->updateCellLinkedList();
    //inner_relation_.updateConfiguration();
}
//=================================================================================================//
RelaxationAccelerationByCMComplex::
RelaxationAccelerationByCMComplex(ComplexRelation& complex_relation)
    : LocalDynamics(complex_relation.getSPHBody()), 
    RelaxDataDelegateComplex(complex_relation),
    acc_(particles_->acc_), pos_(particles_->pos_),
    B_(*particles_->template getVariableByName<Matd>("CorrectionMatrix"))
{
    contact_B_.resize(this->contact_particles_.size());
    for (size_t k = 0; k != this->contact_particles_.size(); ++k)
    {
        contact_B_[k] = this->contact_particles_[k]->template registerSharedVariable<Matd>("CorrectionMatrix");
    }
}
//=================================================================================================//
RelaxationAccelerationByCMComplexWithLevelSetCorrection::
RelaxationAccelerationByCMComplexWithLevelSetCorrection(ComplexRelation& complex_relation, const std::string& shape_name)
    : RelaxationAccelerationByCMComplex(complex_relation),
    sph_adaptation_(sph_body_.sph_adaptation_)
{
    ComplexShape& complex_shape = DynamicCast<ComplexShape>(this, *sph_body_.body_shape_);
    level_set_shape_ = DynamicCast<LevelSetShape>(this, complex_shape.getShapeByName(shape_name));
}
//=================================================================================================//
RelaxationStepByCMComplex::RelaxationStepByCMComplex(ComplexRelation& complex_relation,
    const std::string& shape_name, bool level_set_correction)
    : BaseDynamics<void>(complex_relation.getSPHBody()),
    real_body_(complex_relation.getInnerRelation().real_body_),
    complex_relation_(complex_relation),
    near_shape_surface_(*real_body_, shape_name),
    get_time_step_square_(*real_body_), update_particle_position_(*real_body_),
    surface_bounding_(near_shape_surface_)
{
    if (!level_set_correction)
    {
        relaxation_acceleration_complex_ =
            makeUnique<InteractionDynamics<RelaxationAccelerationByCMComplex>>(complex_relation);
    }
    else
    {
        relaxation_acceleration_complex_ =
            makeUnique<InteractionDynamics<RelaxationAccelerationByCMComplexWithLevelSetCorrection>>(complex_relation, shape_name);
    }
}
//=================================================================================================//
void RelaxationStepByCMComplex::exec(Real dt)
{
    real_body_->updateCellLinkedList();
    complex_relation_.updateConfiguration();
    relaxation_acceleration_complex_->exec();
    Real dt_square = get_time_step_square_.exec();
    update_particle_position_.exec(dt_square);
    surface_bounding_.exec();
}
//****************************************IMPLICIT SCHEME******************************************//
//=================================================================================================//
RelaxationImplicitInner::RelaxationImplicitInner(BaseInnerRelation& inner_relation)
	: LocalDynamics(inner_relation.getSPHBody()), RelaxDataDelegateInner(inner_relation),
	kernel_(inner_relation.getSPHBody().sph_adaptation_->getKernel()), target_residual_pressure_(1.0),
	Vol_(particles_->Vol_), pos_(particles_->pos_), acc_(particles_->acc_),
	sph_adaptation_(sph_body_.sph_adaptation_)
{
    particles_->registerVariable(residual_pressure_, "ResidualPressure");
    particles_->addVariableToWrite<Real>("ResidualPressure");
    particles_->registerVariable(implicit_residual_pressure_, "ImplicitResidualPressure");
    particles_->addVariableToWrite<Real>("ImplicitResidualPressure");
    particles_->addVariableToWrite<Real>("VolumetricMeasure");
    level_set_shape_ = DynamicCast<LevelSetShape>(this, sph_body_.body_shape_);
};
//=================================================================================================//
ErrorAndParameters<Vecd, Matd, Matd> RelaxationImplicitInner::computeErrorAndParameters(size_t index_i, Real dt)
{
	ErrorAndParameters<Vecd, Matd, Matd> error_and_parameters;
	Neighborhood& inner_neighborhood = inner_configuration_[index_i];

	for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
	{
		size_t index_j = inner_neighborhood.j_[n];
		Matd parameter_b = 2.0 * inner_neighborhood.e_ij_[n] * inner_neighborhood.e_ij_[n].transpose() *
			               kernel_->d2W(inner_neighborhood.r_ij_[n], inner_neighborhood.e_ij_[n]) * Vol_[index_j] * dt * dt;

		error_and_parameters.error_ += 2.0 * inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n] * dt * dt;
		error_and_parameters.a_ -= parameter_b;
		error_and_parameters.c_ += parameter_b * parameter_b;
	}

	Matd evolution = Matd::Identity();
	error_and_parameters.a_ -= evolution;
	return error_and_parameters;
}
//=================================================================================================//
void RelaxationImplicitInner::updateStates(size_t index_i, Real dt,
	const ErrorAndParameters<Vecd, Matd, Matd>& error_and_parameters)
{
	Matd parameter_l = error_and_parameters.a_ * error_and_parameters.a_ + error_and_parameters.c_;
	Vecd parameter_k = parameter_l.inverse() * error_and_parameters.error_;

	pos_[index_i] += error_and_parameters.a_ * parameter_k;

	Neighborhood& inner_neighborhood = inner_configuration_[index_i];
	for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
	{
		size_t index_j = inner_neighborhood.j_[n];
		Matd parameter_b = 2.0 * inner_neighborhood.e_ij_[n] * inner_neighborhood.e_ij_[n].transpose() *
			               kernel_->d2W(inner_neighborhood.r_ij_[n], inner_neighborhood.e_ij_[n]) * Vol_[index_j] * dt * dt;
		pos_[index_j] -= parameter_b * parameter_k;
	}
}
//=================================================================================================//
void RelaxationImplicitInner::interaction(size_t index_i, Real dt)
{
	ErrorAndParameters<Vecd, Matd, Matd> error_and_parameters = computeErrorAndParameters(index_i, dt);
	updateStates(index_i, dt, error_and_parameters);
	acc_[index_i] = -error_and_parameters.error_ / dt / dt;
	implicit_residual_pressure_[index_i] = (error_and_parameters.error_ / dt / dt).norm();
    residual_pressure_[index_i] = error_and_parameters.error_.norm();
}
//=================================================================================================//
RelaxationImplicitInnerWithLevelSetCorrection::RelaxationImplicitInnerWithLevelSetCorrection(BaseInnerRelation& inner_relation)
	: RelaxationImplicitInner(inner_relation) {}
//=================================================================================================//
ErrorAndParameters<Vecd, Matd, Matd> RelaxationImplicitInnerWithLevelSetCorrection::computeErrorAndParameters(size_t index_i, Real dt)
{
	ErrorAndParameters<Vecd, Matd, Matd> error_and_parameters = RelaxationImplicitInner::computeErrorAndParameters(index_i, dt);
	error_and_parameters.error_ += 2.0 * level_set_shape_->computeKernelGradientIntegral(pos_[index_i],
		sph_adaptation_->SmoothingLengthRatio(index_i)) * dt * dt;
	error_and_parameters.a_ -= 2.0 * level_set_shape_->computeKernelSecondGradientIntegral(pos_[index_i],
		sph_adaptation_->SmoothingLengthRatio(index_i)) * dt * dt;
	return error_and_parameters;
}
//=================================================================================================//
RelaxationStepImplicitInner::RelaxationStepImplicitInner(BaseInnerRelation& inner_relation, bool level_set_correction)
	: BaseDynamics<void>(inner_relation.getSPHBody()), real_body_(inner_relation.real_body_),
	inner_relation_(inner_relation), time_step_size_(0.01), target_residual_pressure_(1.0),
	near_shape_surface_(*real_body_), get_time_step_(*real_body_),
	relaxation_evolution_inner_(inner_relation), surface_bounding_(near_shape_surface_),
	surface_correction_(near_shape_surface_), update_averaged_error_(inner_relation.getSPHBody(), "ResidualPressure") {}
//=================================================================================================//
void RelaxationStepImplicitInner::exec(Real dt)
{
	//real_body_->updateCellLinkedList();
	//inner_relation_.updateConfiguration();
    relaxation_evolution_inner_.exec(time_step_size_);
    time_step_size_ = 50 * sqrt(get_time_step_.exec());
	//surface_correction_.exec(); 
	//surface_bounding_.exec();
	//target_residual_pressure_ = update_averaged_error_.exec();
	//relaxation_evolution_inner_.updateTargetError(target_residual_pressure_);
}
//=================================================================================================//
RelaxationByCMImplicitInner::RelaxationByCMImplicitInner(BaseInnerRelation& inner_relation)
	: LocalDynamics(inner_relation.getSPHBody()), RelaxDataDelegateInner(inner_relation),
    kernel_(inner_relation.getSPHBody().sph_adaptation_->getKernel()),
	Vol_(particles_->Vol_), pos_(particles_->pos_), acc_(particles_->acc_),
	B_(*particles_->template getVariableByName<Matd>("CorrectionMatrix")),
	sph_adaptation_(sph_body_.sph_adaptation_)
{
	particles_->registerVariable(implicit_residual_cm_, "ImplicitResidualCM");
	particles_->addVariableToWrite<Real>("ImplicitResidualCM");
    particles_->addVariableToWrite<Real>("VolumetricMeasure");
	level_set_shape_ = DynamicCast<LevelSetShape>(this, sph_body_.body_shape_);
};
//=================================================================================================//
ErrorAndParameters<Vecd, Matd, Matd> RelaxationByCMImplicitInner::computeErrorAndParameters(size_t index_i, Real dt)
{
	ErrorAndParameters<Vecd, Matd, Matd> error_and_parameters;
	Neighborhood& inner_neighborhood = inner_configuration_[index_i];

	for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
	{
		size_t index_j = inner_neighborhood.j_[n];
		Matd parameter_b = (B_[index_i] + B_[index_j]) * inner_neighborhood.e_ij_[n] * inner_neighborhood.e_ij_[n].transpose() *
			                kernel_->d2W(inner_neighborhood.r_ij_[n], inner_neighborhood.e_ij_[n]) * Vol_[index_j] * dt * dt;

		error_and_parameters.error_ += (B_[index_i] + B_[index_j]) * inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n] * dt * dt;
		error_and_parameters.a_ -= parameter_b;
		error_and_parameters.c_ += parameter_b * parameter_b;
	}

	Matd evolution = Matd::Identity();
	error_and_parameters.a_ -= evolution;
	return error_and_parameters;
}
//=================================================================================================//
void RelaxationByCMImplicitInner::updateStates(size_t index_i, Real dt,
	const ErrorAndParameters<Vecd, Matd, Matd>& error_and_parameters)
{
	Matd parameter_l = error_and_parameters.a_ * error_and_parameters.a_ + error_and_parameters.c_;
	Vecd parameter_k = parameter_l.inverse() * error_and_parameters.error_;

	pos_[index_i] += error_and_parameters.a_ * parameter_k;

	Neighborhood& inner_neighborhood = inner_configuration_[index_i];
	for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
	{
		size_t index_j = inner_neighborhood.j_[n];
		Matd parameter_b = (B_[index_i] + B_[index_j]) * inner_neighborhood.e_ij_[n] * inner_neighborhood.e_ij_[n].transpose() *
			                kernel_->d2W(inner_neighborhood.r_ij_[n], inner_neighborhood.e_ij_[n]) * Vol_[index_j] * dt * dt;
		pos_[index_j] -= parameter_b * parameter_k;
	}
}
//=================================================================================================//
void RelaxationByCMImplicitInner::interaction(size_t index_i, Real dt)
{
	ErrorAndParameters<Vecd, Matd, Matd> error_and_parameters = computeErrorAndParameters(index_i, dt);
	updateStates(index_i, dt, error_and_parameters);
	acc_[index_i] = -error_and_parameters.error_ / dt / dt;
	implicit_residual_cm_[index_i] = (error_and_parameters.error_ / dt / dt).norm();
}
//=================================================================================================//
RelaxationByCMImplicitInnerWithLevelSetCorrection::
RelaxationByCMImplicitInnerWithLevelSetCorrection(BaseInnerRelation& inner_relation)
	: RelaxationByCMImplicitInner(inner_relation) {};
//=================================================================================================//
ErrorAndParameters<Vecd, Matd, Matd> RelaxationByCMImplicitInnerWithLevelSetCorrection::
computeErrorAndParameters(size_t index_i, Real dt)
{
	ErrorAndParameters<Vecd, Matd, Matd> error_and_parameters = RelaxationByCMImplicitInner::computeErrorAndParameters(index_i, dt);

	error_and_parameters.error_ += (B_[index_i] + B_[index_i]) * level_set_shape_->computeKernelGradientIntegral(
		pos_[index_i], sph_adaptation_->SmoothingLengthRatio(index_i)) * dt * dt;
	error_and_parameters.a_ -= (B_[index_i] + B_[index_i]) * level_set_shape_->computeKernelSecondGradientIntegral(
		pos_[index_i], sph_adaptation_->SmoothingLengthRatio(index_i)) * dt * dt;
	return error_and_parameters;
}
//=================================================================================================//
RelaxationStepByCMImplicitInner::RelaxationStepByCMImplicitInner(BaseInnerRelation &inner_relation, bool level_set_correction)
    : BaseDynamics<void>(inner_relation.getSPHBody()),
      time_step_size_(0.01), real_body_(inner_relation.real_body_),
      inner_relation_(inner_relation), near_shape_surface_(*real_body_),
      get_time_step_(*real_body_), relaxation_evolution_inner_(inner_relation),
      surface_bounding_(near_shape_surface_), surface_correction_(near_shape_surface_){};
//=================================================================================================//
void RelaxationStepByCMImplicitInner::exec(Real dt)
{
	relaxation_evolution_inner_.exec(time_step_size_);
	time_step_size_ = 50 * sqrt(get_time_step_.exec());
	surface_correction_.exec();
	//surface_bounding_.exec();
    real_body_->updateCellLinkedList();
    inner_relation_.updateConfiguration();
}
//=================================================================================================//
CorrectedConfigurationInnerWithLevelSet::
    CorrectedConfigurationInnerWithLevelSet(BaseInnerRelation &inner_relation, bool level_set_correction)
    : LocalDynamics(inner_relation.getSPHBody()),
      RelaxDataDelegateInner(inner_relation),
      pos_(particles_->pos_),
      B_(*particles_->getVariableByName<Matd>("CorrectionMatrix")),
      level_set_correction_(level_set_correction),
      sph_adaptation_(sph_body_.sph_adaptation_)
{
    level_set_shape_ = DynamicCast<LevelSetShape>(this, sph_body_.body_shape_);
}
 //=================================================================================================//     
void CorrectedConfigurationInnerWithLevelSet::interaction(size_t index_i, Real dt)
{
    Matd local_configuration = Eps * Matd::Identity();

    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        Vecd gradW_ij = inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
        Vecd r_ji = inner_neighborhood.r_ij_[n] * inner_neighborhood.e_ij_[n];
        local_configuration -= r_ji * gradW_ij.transpose();
    }

    if (level_set_correction_)
    {
        local_configuration -= level_set_shape_->computeDisplacementKernelGradientIntegral(pos_[index_i],
                               sph_adaptation_->SmoothingLengthRatio(index_i));
    }

    B_[index_i] = local_configuration;

}
//=================================================================================================//
void CorrectedConfigurationInnerWithLevelSet::update(size_t index_i, Real dt)
{
    B_[index_i] = B_[index_i].inverse();
}
//=================================================================================================//
UpdateParticleKineticEnergy::
UpdateParticleKineticEnergy(BaseInnerRelation& inner_relation) :
	LocalDynamics(inner_relation.getSPHBody()), GeneralDataDelegateInner(inner_relation),
	mass_(particles_->mass_), acc_(particles_->acc_)
{
	particles_->registerVariable(particle_kinetic_energy, "ParticleKineticEnergy");
	particles_->addVariableToWrite<Real>("ParticleKineticEnergy");
};
//=================================================================================================//
void UpdateParticleKineticEnergy::interaction(size_t index_i, Real dt)
{
	particle_kinetic_energy[index_i] = acc_[index_i].norm(); /* L2 norm. */
};
//=================================================================================================//
CheckCorrectedZeroOrderConsistency::
CheckCorrectedZeroOrderConsistency(BaseInnerRelation& inner_relation, bool level_set_correction) :
	LocalDynamics(inner_relation.getSPHBody()), GeneralDataDelegateInner(inner_relation),
	level_set_correction_(level_set_correction), pos_(particles_->pos_),
	B_(*particles_->template getVariableByName<Matd>("CorrectionMatrix")),
	sph_adaptation_(sph_body_.sph_adaptation_)
{
	particles_->registerVariable(corrected_zero_order_error_norm_, "CorrectedZeroOrderErrorNorm");
	particles_->addVariableToWrite<Real>("CorrectedZeroOrderErrorNorm");
	particles_->registerVariable(corrected_zero_order_error_, "CorrectedZeroOrderError");
	particles_->addVariableToWrite<Vecd>("CorrectedZeroOrderError");
    level_set_shape_ = DynamicCast<LevelSetShape>(this, sph_body_.body_shape_);
}
//=================================================================================================//
void CheckCorrectedZeroOrderConsistency::interaction(size_t index_i, Real dt)
{
	Vecd acceleration = Vecd::Zero();
	Neighborhood& inner_neighborhood = inner_configuration_[index_i];
	for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
	{
		size_t index_j = inner_neighborhood.j_[n];
		acceleration -= (B_[index_i] + B_[index_j]) * inner_neighborhood.e_ij_[n] * inner_neighborhood.dW_ijV_j_[n];
	}

	/*if (level_set_correction_)
	{
		acceleration -= 0.5 * (B_[index_i]  + B_[index_i]) * level_set_shape_->computeKernelGradientIntegral(
						pos_[index_i], sph_adaptation_->SmoothingLengthRatio(index_i));
	}*/

	corrected_zero_order_error_[index_i] = acceleration;
	corrected_zero_order_error_norm_[index_i] = acceleration.norm();
}
//=================================================================================================//
CheckCorrectedFirstOrderConsistency::
CheckCorrectedFirstOrderConsistency(BaseInnerRelation& inner_relation, bool level_set_correction) :
	LocalDynamics(inner_relation.getSPHBody()), GeneralDataDelegateInner(inner_relation),
	level_set_correction_(level_set_correction), pos_(particles_->pos_),
	B_(*particles_->template getVariableByName<Matd>("CorrectionMatrix")),
	sph_adaptation_(sph_body_.sph_adaptation_)
{
	particles_->registerVariable(corrected_first_order_error_norm_, "CorrectedFirstOrderErrorNorm");
	particles_->addVariableToWrite<Real>("CorrectedFirstOrderErrorNorm");
	particles_->registerVariable(corrected_first_order_error_, "CorrectedFirstOrderError");
	particles_->addVariableToWrite<Matd>("CorrectedFirstOrderError");
    level_set_shape_ = DynamicCast<LevelSetShape>(this, sph_body_.body_shape_);
}
//=================================================================================================//
void CheckCorrectedFirstOrderConsistency::interaction(size_t index_i, Real dt)
{
	Matd acceleration = Matd::Zero();
	Neighborhood& inner_neighborhood = inner_configuration_[index_i];
	for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
	{
		size_t index_j = inner_neighborhood.j_[n];
		acceleration -= 0.5 * (B_[index_i] + B_[index_j]) * inner_neighborhood.r_ij_[n] * inner_neighborhood.dW_ijV_j_[n] *
			inner_neighborhood.e_ij_[n] * (inner_neighborhood.e_ij_[n]).transpose();
	}

	/*if (level_set_correction_)
	{
		acceleration -= 0.5 * (B_[index_i] + B_[index_i]) * level_set_shape_->computeDisplacementKernelGradientIntegral(
						pos_[index_i], sph_adaptation_->SmoothingLengthRatio(index_i));
	}*/

	corrected_first_order_error_[index_i] = acceleration;
	corrected_first_order_error_norm_[index_i] = (acceleration - Matd::Identity()).norm();
}
//=================================================================================================//
CheckConsistencyRealization::
CheckConsistencyRealization(BaseInnerRelation& inner_relation, bool level_set_correction) :
    LocalDynamics(inner_relation.getSPHBody()), GeneralDataDelegateInner(inner_relation),
    level_set_correction_(level_set_correction), pos_(particles_->pos_),
    pressure_(*particles_->template getVariableByName<Real>("Pressure")),
    B_(*particles_->template getVariableByName<Matd>("CorrectionMatrix")),
    sph_adaptation_(sph_body_.sph_adaptation_)
{
    particles_->registerVariable(relaxation_error_, "RelaxationErrorNorm");
    particles_->addVariableToWrite<Real>("RelaxationErrorNorm");
    particles_->registerVariable(pressure_gradient_error_norm_, "PressureGradientErrorNorm");
    particles_->addVariableToWrite<Real>("PressureGradientErrorNorm");
    particles_->registerVariable(pressure_gradient_, "PressureGradient");
    particles_->addVariableToWrite<Vecd>("PressureGradient");
    particles_->registerVariable(zero_order_error_norm_, "ZeroOrderErrorNorm");
    particles_->addVariableToWrite<Real>("ZeroOrderErrorNorm");
    particles_->registerVariable(zero_order_error_, "ZeroOrderError");
    particles_->addVariableToWrite<Vecd>("ZeroOrderError");
    particles_->registerVariable(reproduce_gradient_error_norm_, "ReproduceGradientErrorNorm");
    particles_->addVariableToWrite<Real>("ReproduceGradientErrorNorm");
    particles_->registerVariable(reproduce_gradient_, "ReproduceGradient");
    particles_->addVariableToWrite<Vecd>("ReproduceGradient");
    particles_->addVariableToWrite<Real>("Pressure");

    level_set_shape_ = DynamicCast<LevelSetShape>(this, sph_body_.body_shape_);
}
//=================================================================================================//
void CheckConsistencyRealization::interaction(size_t index_i, Real dt)
{
    Vecd relaxation_error = Vecd::Zero();
    Vecd pressure_gradient = Vecd::Zero();
    Vecd zero_order_error = Vecd::Zero();
    Vecd reproduce_gradient = Vecd::Zero();
    Neighborhood& inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        relaxation_error += pressure_[index_i] * inner_neighborhood.e_ij_[n] * inner_neighborhood.dW_ijV_j_[n];
        pressure_gradient += (pressure_[index_i] * B_[index_i] + pressure_[index_j] * B_[index_j]) * 
                              inner_neighborhood.e_ij_[n] * inner_neighborhood.dW_ijV_j_[n];
        zero_order_error += pressure_[index_i] * (B_[index_i] + B_[index_j]) *
                            inner_neighborhood.e_ij_[n] * inner_neighborhood.dW_ijV_j_[n];
        reproduce_gradient -= (pressure_[index_i] - pressure_[index_j]) * B_[index_j] *
                               inner_neighborhood.e_ij_[n] * inner_neighborhood.dW_ijV_j_[n];
    }
    
    relaxation_error_[index_i] = relaxation_error.norm();
    pressure_gradient_[index_i] = pressure_gradient;
    //pressure_gradient_error_norm_[index_i] = sqrt(pow((pressure_gradient[0] - 1), 2) + pow((pressure_gradient[1] - 0), 2));
    pressure_gradient_error_norm_[index_i] = sqrt(pow((pressure_gradient[0] - (cos(pos_[index_i][0] * 2 * Pi))), 2) + pow((pressure_gradient[1] - 0), 2));
    zero_order_error_[index_i] = zero_order_error;
    zero_order_error_norm_[index_i] = zero_order_error.norm();
    reproduce_gradient_[index_i] = reproduce_gradient;
    //reproduce_gradient_error_norm_[index_i] = sqrt(pow((reproduce_gradient[0] - 1), 2) + pow((reproduce_gradient[1] - 0), 2));
    reproduce_gradient_error_norm_[index_i] = sqrt(pow((reproduce_gradient[0] - (cos(pos_[index_i][0] * 2 * Pi))), 2) + pow((reproduce_gradient[1] - 0), 2));
}
//=================================================================================================//
CheckReverseConsistencyRealization::
CheckReverseConsistencyRealization(BaseInnerRelation& inner_relation, bool level_set_correction) :
    LocalDynamics(inner_relation.getSPHBody()), GeneralDataDelegateInner(inner_relation),
    level_set_correction_(level_set_correction), pos_(particles_->pos_),
    pressure_(*particles_->template getVariableByName<Real>("Pressure")),
    B_(*particles_->template getVariableByName<Matd>("CorrectionMatrix")),
    sph_adaptation_(sph_body_.sph_adaptation_)
{
    particles_->registerVariable(pressure_gradient_error_norm_, "PressureGradientErrorNormReverse");
    particles_->addVariableToWrite<Real>("PressureGradientErrorNormReverse");
    particles_->registerVariable(pressure_gradient_, "PressureGradientReverse");
    particles_->addVariableToWrite<Vecd>("PressureGradientReverse");
    particles_->registerVariable(zero_order_error_norm_, "ZeroOrderErrorNormReverse");
    particles_->addVariableToWrite<Real>("ZeroOrderErrorNormReverse");
    particles_->registerVariable(zero_order_error_, "ZeroOrderErrorReverse");
    particles_->addVariableToWrite<Vecd>("ZeroOrderErrorReverse");
    particles_->registerVariable(reproduce_gradient_error_norm_, "ReproduceGradientErrorNormReverse");
    particles_->addVariableToWrite<Real>("ReproduceGradientErrorNormReverse");
    particles_->registerVariable(reproduce_gradient_, "ReproduceGradientReverse");
    particles_->addVariableToWrite<Vecd>("ReproduceGradientReverse");
    level_set_shape_ = DynamicCast<LevelSetShape>(this, sph_body_.body_shape_);
}
//=================================================================================================//
void CheckReverseConsistencyRealization::interaction(size_t index_i, Real dt)
{
    Neighborhood& inner_neighborhood = inner_configuration_[index_i];
    Vecd pressure_gradient = Vecd::Zero();
    Vecd zero_order_error = Vecd::Zero();
    Vecd reproduce_gradient = Vecd::Zero();
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Vecd gradW_ij = inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
        pressure_gradient += (pressure_[index_i] * B_[index_j] + pressure_[index_j] * B_[index_i]) * 
                             inner_neighborhood.e_ij_[n] * inner_neighborhood.dW_ijV_j_[n];
        zero_order_error += pressure_[index_i] * (B_[index_i] + B_[index_j]) *
                            inner_neighborhood.e_ij_[n] * inner_neighborhood.dW_ijV_j_[n];
        reproduce_gradient -= (pressure_[index_i] - pressure_[index_j]) * B_[index_i] * gradW_ij;
    }

    pressure_gradient_[index_i] = pressure_gradient;
    //pressure_gradient_error_norm_[index_i] = sqrt(pow((pressure_gradient[0] - 1), 2) + pow((pressure_gradient[1] - 0), 2));
    pressure_gradient_error_norm_[index_i] = sqrt(pow((pressure_gradient[0] - (cos(pos_[index_i][0] * 2 * Pi))), 2) + pow((pressure_gradient[1] - 0), 2));
    zero_order_error_[index_i] = zero_order_error;
    zero_order_error_norm_[index_i] = zero_order_error.norm();
    reproduce_gradient_[index_i] = reproduce_gradient;
    //reproduce_gradient_error_norm_[index_i] = sqrt(pow((reproduce_gradient[0] - 1), 2) + pow((reproduce_gradient[1] - 0), 2));
    reproduce_gradient_error_norm_[index_i] = sqrt(pow((reproduce_gradient[0] - (cos(pos_[index_i][0] * 2 * Pi))), 2) + pow((reproduce_gradient[1] - 0), 2));
}
//=================================================================================================//
ShellMidSurfaceBounding::
    ShellMidSurfaceBounding(NearShapeSurface &body_part, BaseInnerRelation &inner_relation,
                            Real thickness, Real level_set_refinement_ratio)
    : BaseLocalDynamics<BodyPartByCell>(body_part), RelaxDataDelegateInner(inner_relation),
      pos_(particles_->pos_), constrained_distance_(0.5 * sph_body_.sph_adaptation_->MinimumSpacing()),
      particle_spacing_ref_(sph_body_.sph_adaptation_->MinimumSpacing()),
      thickness_(thickness), level_set_refinement_ratio_(level_set_refinement_ratio),
      level_set_shape_(DynamicCast<LevelSetShape>(this, sph_body_.body_shape_)) {}
//=================================================================================================//
void ShellMidSurfaceBounding::update(size_t index_i, Real dt)
{
    Vecd none_normalized_normal = level_set_shape_->findLevelSetGradient(pos_[index_i]);
    Vecd normal = level_set_shape_->findNormalDirection(pos_[index_i]);
    Real factor = none_normalized_normal.squaredNorm() / level_set_refinement_ratio_;
    pos_[index_i] -= factor * constrained_distance_ * normal;
}
//=================================================================================================//
ShellNormalDirectionPrediction::
    ShellNormalDirectionPrediction(BaseInnerRelation &inner_relation,
                                   Real thickness, Real consistency_criterion)
    : BaseDynamics<void>(inner_relation.getSPHBody()),
      convergence_criterion_(cos(0.01 * Pi)),
      consistency_criterion_(consistency_criterion),
      normal_prediction_(inner_relation.getSPHBody(), thickness),
      normal_prediction_convergence_check_(inner_relation.getSPHBody(), convergence_criterion_),
      consistency_correction_(inner_relation, consistency_criterion_),
      consistency_updated_check_(inner_relation.getSPHBody()),
      smoothing_normal_(inner_relation) {}
//=================================================================================================//
void ShellNormalDirectionPrediction::exec(Real dt)
{
    predictNormalDirection();
    correctNormalDirection();
    predictNormalDirection();
    smoothing_normal_.exec();
}
//=================================================================================================//
void ShellNormalDirectionPrediction::predictNormalDirection()
{
    bool prediction_convergence = false;
    size_t ite_predict = 0;
    while (!prediction_convergence)
    {
        normal_prediction_.exec();
        prediction_convergence = normal_prediction_convergence_check_.exec();
        if (ite_predict > 100)
        {
            std::cout << "\n Error: class ShellNormalDirectionPrediction normal prediction not converged after 100 iterations." << std::endl;
            std::cout << __FILE__ << ':' << __LINE__ << std::endl;
            exit(1);
        }

        ite_predict++;
    }
    std::cout << "\n Information: normal direction prediction converged after '" << ite_predict << "' steps." << std::endl;
}
//=================================================================================================//
void ShellNormalDirectionPrediction::correctNormalDirection()
{
    bool consistency_updated = false;
    size_t ite_updated = 0;
    while (!consistency_updated)
    {
        consistency_correction_.exec();
        consistency_updated = consistency_updated_check_.exec();
        if (ite_updated > 100)
        {
            std::cout << "\n Error: class ShellNormalDirectionPrediction normal consistency not updated  after 100 iterations." << std::endl;
            std::cout << __FILE__ << ':' << __LINE__ << std::endl;
            exit(1);
        }
        ite_updated++;
    }
    std::cout << "\n Information: normal consistency updated after '" << ite_updated << "' steps." << std::endl;
}
//=================================================================================================//
ShellNormalDirectionPrediction::NormalPrediction::NormalPrediction(SPHBody &sph_body, Real thickness)
    : RelaxDataDelegateSimple(sph_body), LocalDynamics(sph_body), thickness_(thickness),
      level_set_shape_(DynamicCast<LevelSetShape>(this, sph_body.body_shape_)),
      pos_(particles_->pos_), n_(*particles_->getVariableByName<Vecd>("NormalDirection"))
{
    particles_->registerVariable(n_temp_, "PreviousNormalDirection", [&](size_t i) -> Vecd
                                 { return n_[i]; });
}
//=================================================================================================//
void ShellNormalDirectionPrediction::NormalPrediction::update(size_t index_i, Real dt)
{
    n_temp_[index_i] = n_[index_i];
    n_[index_i] = level_set_shape_->findNormalDirection(pos_[index_i] + 0.3 * thickness_ * n_temp_[index_i]);
}
//=================================================================================================//
ShellNormalDirectionPrediction::PredictionConvergenceCheck::
    PredictionConvergenceCheck(SPHBody &sph_body, Real convergence_criterion)
    : LocalDynamicsReduce<bool, ReduceAND>(sph_body, true), RelaxDataDelegateSimple(sph_body),
      convergence_criterion_(convergence_criterion), n_(*particles_->getVariableByName<Vecd>("NormalDirection")),
      n_temp_(*particles_->getVariableByName<Vecd>("PreviousNormalDirection")) {}
//=================================================================================================//
bool ShellNormalDirectionPrediction::PredictionConvergenceCheck::reduce(size_t index_i, Real dt)
{
    return n_[index_i].dot(n_temp_[index_i]) > convergence_criterion_;
}
//=================================================================================================//
ShellNormalDirectionPrediction::ConsistencyCorrection::
    ConsistencyCorrection(BaseInnerRelation &inner_relation, Real consistency_criterion)
    : LocalDynamics(inner_relation.getSPHBody()), RelaxDataDelegateInner(inner_relation),
      consistency_criterion_(consistency_criterion),
      n_(*particles_->getVariableByName<Vecd>("NormalDirection"))
{
    particles_->registerVariable(updated_indicator_, "UpdatedIndicator", [&](size_t i) -> int
                                 { return 0; });
    updated_indicator_[particles_->total_real_particles_ / 3] = 1;
}
//=================================================================================================//
ShellNormalDirectionPrediction::ConsistencyUpdatedCheck::ConsistencyUpdatedCheck(SPHBody &sph_body)
    : LocalDynamicsReduce<bool, ReduceAND>(sph_body, true),
      RelaxDataDelegateSimple(sph_body),
      updated_indicator_(*particles_->getVariableByName<int>("UpdatedIndicator")) {}
//=================================================================================================//
bool ShellNormalDirectionPrediction::ConsistencyUpdatedCheck::reduce(size_t index_i, Real dt)
{
    return updated_indicator_[index_i] != 0;
}
//=================================================================================================//
ShellNormalDirectionPrediction::SmoothingNormal::
    SmoothingNormal(BaseInnerRelation &inner_relation)
    : ParticleSmoothing<Vecd>(inner_relation, "NormalDirection"){};
//=================================================================================================//
void ShellNormalDirectionPrediction::SmoothingNormal::update(size_t index_i, Real dt)
{
    ParticleSmoothing<Vecd>::update(index_i, dt);
    smoothed_[index_i] /= temp_[index_i].norm() + TinyReal;
}
//=================================================================================================//
ShellRelaxationStepInner::
    ShellRelaxationStepInner(BaseInnerRelation &inner_relation, Real thickness,
                             Real level_set_refinement_ratio, bool level_set_correction)
    : RelaxationStepInner(inner_relation, level_set_correction),
      update_shell_particle_position_(*real_body_),
      mid_surface_bounding_(near_shape_surface_, inner_relation,
                            thickness, level_set_refinement_ratio) {}
//=================================================================================================//
void ShellRelaxationStepInner::exec(Real ite_p)
{
    real_body_->updateCellLinkedList();
    inner_relation_.updateConfiguration();
    relaxation_acceleration_inner_->exec();
    Real dt_square = get_time_step_square_.exec();
    update_shell_particle_position_.exec(dt_square);
    mid_surface_bounding_.exec();
}
//=================================================================================================//
} // namespace relax_dynamics
  //=================================================================================================//
} // namespace SPH
//=================================================================================================//
