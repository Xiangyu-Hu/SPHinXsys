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
    //real_body_->updateCellLinkedList();
    //inner_relation_.updateConfiguration();
    relaxation_acceleration_inner_->exec();
    Real dt_square = get_time_step_square_.exec();
    update_particle_position_.exec(dt_square);
    //surface_bounding_.exec();
}
//=================================================================================================//
RelaxationAccelerationByStressInner::RelaxationAccelerationByStressInner(BaseInnerRelation& inner_relation)
    : LocalDynamics(inner_relation.getSPHBody()), RelaxDataDelegateInner(inner_relation),
    acc_(particles_->acc_), pos_(particles_->pos_),
    B_(*particles_->template getVariableByName<Matd>("CorrectionMatrix")) {}
//=================================================================================================//
RelaxationAccelerationByStressInnerWithLevelSetCorrection::
RelaxationAccelerationByStressInnerWithLevelSetCorrection(BaseInnerRelation& inner_relation)
    : RelaxationAccelerationByStressInner(inner_relation), sph_adaptation_(sph_body_.sph_adaptation_)
{
    level_set_shape_ = DynamicCast<LevelSetShape>(this, sph_body_.body_shape_);
}
//=================================================================================================//
RelaxationStepByStressInner::RelaxationStepByStressInner(BaseInnerRelation& inner_relation, bool level_set_correction)
    : BaseDynamics<void>(inner_relation.getSPHBody()), real_body_(inner_relation.real_body_),
    inner_relation_(inner_relation),
    near_shape_surface_(*real_body_), get_time_step_square_(*real_body_),
    update_particle_position_(*real_body_), surface_bounding_(near_shape_surface_)
{
    if (!level_set_correction)
    {
        relaxation_acceleration_inner_ =
            makeUnique<InteractionDynamics<RelaxationAccelerationByStressInner>>(inner_relation);
        surface_correction_ = makeShared<SimpleDynamics<ShapeSurfaceBounding>>(near_shape_surface_);
    }
    else
    {
        relaxation_acceleration_inner_ =
            makeUnique<InteractionDynamics<RelaxationAccelerationByStressInnerWithLevelSetCorrection>>(inner_relation);
        surface_correction_ = makeShared<SimpleDynamics<NearSurfaceVolumeCorrection>>(near_shape_surface_);
    }
}
//=================================================================================================//
void RelaxationStepByStressInner::exec(Real dt)
{
    real_body_->updateCellLinkedList();
    inner_relation_.updateConfiguration();
    relaxation_acceleration_inner_->exec();
    Real dt_square = get_time_step_square_.exec();
    update_particle_position_.exec(dt_square);
    surface_correction_->exec();
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

/****This is the implicit scheme for relaxation****/

//=================================================================================================//
RelaxationImplicitInner::RelaxationImplicitInner(BaseInnerRelation& inner_relation)
	: LocalDynamics(inner_relation.getSPHBody()), RelaxDataDelegateInner(inner_relation),
	kernel_(inner_relation.getSPHBody().sph_adaptation_->getKernel()), target_error_p_(1.0),
	Vol_(particles_->Vol_), pos_(particles_->pos_), acc_(particles_->acc_),
	sph_adaptation_(sph_body_.sph_adaptation_)
{
	particles_->registerVariable(implicit_residue_p_, "implicit_residue_p_");
	particles_->addVariableToWrite<Real>("implicit_residue_p_");
	particles_->addVariableToWrite<Real>("VolumetricMeasure");
	level_set_shape_ = DynamicCast<LevelSetShape>(this, sph_body_.body_shape_);

	particles_->registerVariable(error_p_, "error_p_");
	particles_->addVariableToWrite<Real>("error_p_");
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
	implicit_residue_p_[index_i] = (error_and_parameters.error_ / dt / dt).norm();

	error_p_[index_i] = error_and_parameters.error_.norm();
}
//=================================================================================================//
RelaxationImplicitInnerWithLevelSetCorrection::
RelaxationImplicitInnerWithLevelSetCorrection(BaseInnerRelation& inner_relation)
	: RelaxationImplicitInner(inner_relation) {}
//=================================================================================================//
ErrorAndParameters<Vecd, Matd, Matd> RelaxationImplicitInnerWithLevelSetCorrection::
computeErrorAndParameters(size_t index_i, Real dt)
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
	inner_relation_(inner_relation), time_step_size_(0.01), target_error_p_(1.0),
	near_shape_surface_(*real_body_), get_time_step_(*real_body_),
	relaxation_evolution_inner_(inner_relation), surface_bounding_(near_shape_surface_),
	surface_correction_(near_shape_surface_),
	update_averaged_error_(inner_relation.getSPHBody(), "error_p_") {}
//=================================================================================================//
void RelaxationStepImplicitInner::exec(Real dt)
{
	//real_body_->updateCellLinkedList();
	//inner_relation_.updateConfiguration();
	relaxation_evolution_inner_.exec(dt);
	time_step_size_ = 20 * sqrt(get_time_step_.exec());
	//surface_correction_.exec(); 
	//surface_bounding_.exec();
	target_error_p_ = update_averaged_error_.exec();
	relaxation_evolution_inner_.updateTargetError(target_error_p_);
}
//=================================================================================================//
RelaxationByStressImplicitInner::RelaxationByStressImplicitInner(BaseInnerRelation& inner_relation)
	: LocalDynamics(inner_relation.getSPHBody()), RelaxDataDelegateInner(inner_relation),
	kernel_(inner_relation.getSPHBody().sph_adaptation_->getKernel()), target_error_s_(1.0),
	Vol_(particles_->Vol_), pos_(particles_->pos_), acc_(particles_->acc_),
	B_(*particles_->template getVariableByName<Matd>("CorrectionMatrix")),
	sph_adaptation_(sph_body_.sph_adaptation_)
{
	particles_->registerVariable(implicit_residue_s_, "implicit_residue_s_");
	particles_->addVariableToWrite<Real>("implicit_residue_s_");
	level_set_shape_ = DynamicCast<LevelSetShape>(this, sph_body_.body_shape_);

	particles_->registerVariable(error_s_, "error_s_");
	particles_->addVariableToWrite<Real>("error_s_");
};
//=================================================================================================//
ErrorAndParameters<Vecd, Matd, Matd> RelaxationByStressImplicitInner::computeErrorAndParameters(size_t index_i, Real dt)
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
void RelaxationByStressImplicitInner::updateStates(size_t index_i, Real dt,
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
void RelaxationByStressImplicitInner::interaction(size_t index_i, Real dt)
{
	ErrorAndParameters<Vecd, Matd, Matd> error_and_parameters = computeErrorAndParameters(index_i, dt);
	updateStates(index_i, dt, error_and_parameters);
	acc_[index_i] = -error_and_parameters.error_ / dt / dt;
	implicit_residue_s_[index_i] = (error_and_parameters.error_ / dt / dt).norm();

	error_s_[index_i] = error_and_parameters.error_.norm();
	error_s_x_[index_i] = abs(error_and_parameters.error_[0]);
	error_s_y_[index_i] = abs(error_and_parameters.error_[1]);
}
//=================================================================================================//
RelaxationByStressImplicitInnerWithLevelSetCorrection::
RelaxationByStressImplicitInnerWithLevelSetCorrection(BaseInnerRelation& inner_relation)
	: RelaxationByStressImplicitInner(inner_relation) {};
//=================================================================================================//
ErrorAndParameters<Vecd, Matd, Matd> RelaxationByStressImplicitInnerWithLevelSetCorrection::
computeErrorAndParameters(size_t index_i, Real dt)
{
	ErrorAndParameters<Vecd, Matd, Matd> error_and_parameters = RelaxationByStressImplicitInner::computeErrorAndParameters(index_i, dt);

	error_and_parameters.error_ += (B_[index_i] + B_[index_i]) * level_set_shape_->computeKernelGradientIntegral(
		pos_[index_i], sph_adaptation_->SmoothingLengthRatio(index_i)) * dt * dt;
	error_and_parameters.a_ -= (B_[index_i] + B_[index_i]) * level_set_shape_->computeKernelSecondGradientIntegral(
		pos_[index_i], sph_adaptation_->SmoothingLengthRatio(index_i)) * dt * dt;
	return error_and_parameters;
}
//=================================================================================================//
RelaxationStepByStressImplicitInner::RelaxationStepByStressImplicitInner(BaseInnerRelation& inner_relation, bool level_set_correction)
	: BaseDynamics<void>(inner_relation.getSPHBody()), real_body_(inner_relation.real_body_),
	inner_relation_(inner_relation), time_step_size_(0.01), target_error_s_(1.0),
	near_shape_surface_(*real_body_), get_time_step_(*real_body_),
	relaxation_evolution_inner_(inner_relation), surface_bounding_(near_shape_surface_),
	surface_correction_(near_shape_surface_),
	update_averaged_error_(inner_relation.getSPHBody(), "error_s_") {};
//=================================================================================================//
void RelaxationStepByStressImplicitInner::exec(Real dt)
{
	//real_body_->updateCellLinkedList();
	//inner_relation_.updateConfiguration();
	relaxation_evolution_inner_.exec(dt);
	//time_step_size_ = 50 * sqrt(get_time_step_.exec());
	//surface_correction_.exec();
	//surface_bounding_.exec();
	target_error_s_ = update_averaged_error_.exec();
	relaxation_evolution_inner_.updateTargetError(target_error_s_);
}
//=================================================================================================//
CalculateParticleStress::CalculateParticleStress(BaseInnerRelation& inner_relation, bool level_set_correction)
	: LocalDynamics(inner_relation.getSPHBody()), RelaxDataDelegateInner(inner_relation),
	pos_(particles_->pos_), level_set_correction_(level_set_correction),
	sph_adaptation_(sph_body_.sph_adaptation_)
{
	level_set_shape_ = DynamicCast<LevelSetShape>(this, sph_body_.body_shape_);
	particles_->registerVariable(stress_, "ParticleStress", [&](size_t i) -> Matd { return Matd::Identity(); });
	particles_->addVariableToWrite<Matd>("ParticleStress");
	particles_->registerVariable(B_, "CorrectionMatrix", [&](size_t i) -> Matd {return Matd::Identity(); });
	particles_->addVariableToWrite<Matd>("CorrectionMatrix");
}
//=================================================================================================//
void CalculateParticleStress::interaction(size_t index_i, Real dt)
{
	Matd particle_stress_ = Matd::Zero();
	Neighborhood& inner_neighborhood = inner_configuration_[index_i];
	for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
	{
		size_t index_j = inner_neighborhood.j_[n];
		particle_stress_ -= inner_neighborhood.r_ij_[n] * inner_neighborhood.e_ij_[n] *
			inner_neighborhood.e_ij_[n].transpose() * inner_neighborhood.dW_ijV_j_[n];
	}

	if (level_set_correction_)
	{
		//The calculation of the particle stress is not related to the level set correction.
		particle_stress_ -= level_set_shape_->computeDisplacementKernelGradientIntegral(pos_[index_i],
			sph_adaptation_->SmoothingLengthRatio(index_i));
	}

	stress_[index_i] = particle_stress_;
	B_[index_i] = (Eps * Matd::Identity() + particle_stress_).inverse();
}
//=================================================================================================//
ZeroOrderConsistencyEvolution::
ZeroOrderConsistencyEvolution(BaseInnerRelation& inner_relation, bool level_set_correction)
	:LocalDynamics(inner_relation.getSPHBody()), RelaxDataDelegateInner(inner_relation),
	kernel_(inner_relation.getSPHBody().sph_adaptation_->getKernel()),
	Vol_(particles_->Vol_), pos_(particles_->pos_), acc_(particles_->acc_),
	B_(*particles_->template getVariableByName<Matd>("CorrectionMatrix")),
	level_set_correction_(level_set_correction), target_error_(1.0),
	sph_adaptation_(sph_body_.sph_adaptation_)
{
	level_set_shape_ = DynamicCast<LevelSetShape>(this, sph_body_.body_shape_);
	particles_->registerVariable(zero_order_evolution_residue_, "zero_order_evolution_residue_");
	particles_->addVariableToWrite<Real>("zero_order_evolution_residue_");
}
//=================================================================================================//
ErrorAndParameters<Vecd, Matd, Matd> ZeroOrderConsistencyEvolution::computeErrorAndParameters(size_t index_i, Real dt)
{
	ErrorAndParameters<Vecd, Matd, Matd> error_and_parameters;
	Neighborhood& inner_neighborhood = inner_configuration_[index_i];

	for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
	{
		size_t index_j = inner_neighborhood.j_[n];
		Matd parameter_b = (B_[index_i] + B_[index_j]) * inner_neighborhood.e_ij_[n] * inner_neighborhood.e_ij_[n].transpose() *
			kernel_->d2W(inner_neighborhood.r_ij_[n], inner_neighborhood.e_ij_[n]) * Vol_[index_j] * dt * dt;
		error_and_parameters.error_ += (B_[index_i] + B_[index_j]) * inner_neighborhood.dW_ijV_j_[n] *
			inner_neighborhood.e_ij_[n] * dt * dt;
		error_and_parameters.a_ -= parameter_b;
		error_and_parameters.c_ += parameter_b * parameter_b;
	}

	/*if (level_set_correction_)
	{
		error_and_parameters.error_ += (B_[index_i] + B_[index_i]) * level_set_shape_->computeKernelGradientIntegral(
			pos_[index_i], sph_adaptation_->SmoothingLengthRatio(index_i)) * dt * dt;
		error_and_parameters.a_ -= (B_[index_i] + B_[index_i]) * level_set_shape_->computeKernelSecondGradientIntegral(
			pos_[index_i], sph_adaptation_->SmoothingLengthRatio(index_i)) * dt * dt;
	}*/

	Matd evolution = Matd::Identity();
	error_and_parameters.a_ -= evolution;
	return error_and_parameters;
}
//=================================================================================================//
void ZeroOrderConsistencyEvolution::updateStates(size_t index_i, Real dt,
	const ErrorAndParameters<Vecd, Matd, Matd>& error_and_parameters)
{
	Matd parameter_l = error_and_parameters.a_ * error_and_parameters.a_ + error_and_parameters.c_;
	Vecd parameter_k = parameter_l.inverse() * error_and_parameters.error_;

	//parameter_k = parameter_k * sqrt(sqrt(abs(error_and_parameters.error_.norm()) / target_error_));
	pos_[index_i] += error_and_parameters.a_ * parameter_k;

	Neighborhood& inner_neighborhood = inner_configuration_[index_i];
	for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
	{
		size_t index_j = inner_neighborhood.j_[n];
		Matd parameter_b = (B_[index_i] + B_[index_j]) *
			inner_neighborhood.e_ij_[n] * inner_neighborhood.e_ij_[n].transpose() *
			kernel_->d2W(inner_neighborhood.r_ij_[n], inner_neighborhood.e_ij_[n]) *
			Vol_[index_j] * dt * dt;
		pos_[index_j] -= parameter_b * parameter_k;
	}
}
//=================================================================================================//
void ZeroOrderConsistencyEvolution::interaction(size_t index_i, Real dt)
{
	ErrorAndParameters<Vecd, Matd, Matd> error_and_parameters = computeErrorAndParameters(index_i, dt);
	updateStates(index_i, dt, error_and_parameters);
	acc_[index_i] = -error_and_parameters.error_ / dt / dt;
	zero_order_evolution_residue_[index_i] = (error_and_parameters.error_ / dt / dt).norm();
}
//=================================================================================================//
ZeroOrderEvolutionStep::ZeroOrderEvolutionStep(BaseInnerRelation& inner_relation, bool level_set_correction)
	: BaseDynamics<void>(inner_relation.getSPHBody()), real_body_(inner_relation.real_body_),
	inner_relation_(inner_relation), time_step_size_(0.01), target_error_(1.0),
	near_shape_surface_(*real_body_), get_time_step_(*real_body_),
	zero_order_consistency_evolution_(inner_relation, false),
	surface_bounding_(near_shape_surface_), surface_correction_(near_shape_surface_),
	update_averaged_error_(inner_relation.getSPHBody(), "zero_order_evolution_residue_") {}
//=================================================================================================//
void ZeroOrderEvolutionStep::exec(Real dt)
{
	//real_body_->updateCellLinkedList();
	//inner_relation_.updateConfiguration();
	zero_order_consistency_evolution_.exec(dt);
	//time_step_size_ = 25 * sqrt(get_time_step_.exec());
	//target_error_ = update_averaged_error_.exec();
	//zero_order_consistency_evolution_.updateTargetError(target_error_);
}
//=================================================================================================//
FirstOrderConsistencyEvolution::
FirstOrderConsistencyEvolution(BaseInnerRelation& inner_relation, bool level_set_correction)
	: LocalDynamics(inner_relation.getSPHBody()), RelaxDataDelegateInner(inner_relation),
	kernel_(inner_relation.getSPHBody().sph_adaptation_->getKernel()),
	Vol_(particles_->Vol_), pos_(particles_->pos_),
	B_(*particles_->template getVariableByName<Matd>("CorrectionMatrix")),
	level_set_correction_(level_set_correction),
	sph_adaptation_(sph_body_.sph_adaptation_),
	B_backup_(*particles_->template getVariableByName<Matd>("CorrectionMatrix"))
{
	level_set_shape_ = DynamicCast<LevelSetShape>(this, sph_body_.body_shape_);
	particles_->registerVariable(first_order_evolution_residue_, "first_order_evolution_residue_");
	particles_->addVariableToWrite<Real>("first_order_evolution_residue_");
	particles_->registerVariable(evolution_indicator_, "evolution_indicator_");
	particles_->addVariableToWrite<Real>("evolution_indicator_");
}
//=================================================================================================//
ErrorAndParameters<Matd, Matd, Matd> FirstOrderConsistencyEvolution::
computeErrorAndParameters(size_t index_i, Real dt)
{
	ErrorAndParameters<Matd, Matd, Matd> error_and_parameters;
	Neighborhood& inner_neighborhood = inner_configuration_[index_i];

	for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
	{
		size_t index_j = inner_neighborhood.j_[n];
		Matd variable_derivative = B_[index_i] + B_[index_j];
		Matd parameter_b = -0.5 * inner_neighborhood.r_ij_[n] * inner_neighborhood.dW_ijV_j_[n] *
			inner_neighborhood.e_ij_[n] * (inner_neighborhood.e_ij_[n]).transpose() * dt;

		error_and_parameters.error_ -= variable_derivative * parameter_b; /* correction matrix is left multiplied. */
		error_and_parameters.a_ += parameter_b;
		error_and_parameters.c_ += parameter_b * parameter_b;
	}

	/*if (level_set_correction_)
	{
		error_and_parameters.a_ -= level_set_shape_->computeDisplacementKernelGradientIntegral(
								   pos_[index_i], sph_adaptation_->SmoothingLengthRatio(index_i)) * dt;
		error_and_parameters.error_ += B_[index_i] * level_set_shape_->computeDisplacementKernelGradientIntegral(
									   pos_[index_i], sph_adaptation_->SmoothingLengthRatio(index_i)) * dt; // correction matrix is left multiplied. //
	}*/

	error_and_parameters.a_ -= Matd::Identity();
	error_and_parameters.error_ += Matd::Identity() * dt;
	return error_and_parameters;
}
//=================================================================================================//
void FirstOrderConsistencyEvolution::updateStates(size_t index_i, Real dt,
	const ErrorAndParameters<Matd, Matd, Matd>& error_and_parameters)
{
	Matd parameter_l = error_and_parameters.a_ * error_and_parameters.a_ + error_and_parameters.c_;
	Matd parameter_k = parameter_l.inverse() * error_and_parameters.error_; /* the learning rate is right multiplied. */
	B_backup_[index_i] = B_[index_i];
	B_[index_i] -= error_and_parameters.a_ * parameter_k;

	Neighborhood& inner_neighborhood = inner_configuration_[index_i];
	for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
	{
		size_t index_j = inner_neighborhood.j_[n];
		Matd parameter_b = -0.5 * inner_neighborhood.r_ij_[n] * inner_neighborhood.dW_ijV_j_[n] *
			inner_neighborhood.e_ij_[n] * (inner_neighborhood.e_ij_[n]).transpose() * dt;
		B_backup_[index_j] = B_[index_i];
		B_[index_j] -= parameter_b * parameter_b;
	}
}
//=================================================================================================//
void FirstOrderConsistencyEvolution::interaction(size_t index_i, Real dt)
{
	ErrorAndParameters<Matd, Matd, Matd> error_and_parameters = computeErrorAndParameters(index_i, dt);
	updateStates(index_i, dt, error_and_parameters);
	first_order_evolution_residue_[index_i] = error_and_parameters.error_.norm();

	ErrorAndParameters<Matd, Matd, Matd> error_and_parameters_after = computeErrorAndParameters(index_i, dt);
	if (error_and_parameters_after.error_.norm() > error_and_parameters.error_.norm())
	{
		B_[index_i] = B_backup_[index_i];
		Neighborhood& inner_neighborhood = inner_configuration_[index_i];
		for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
		{
			size_t index_j = inner_neighborhood.j_[n];
			B_[index_j] = B_backup_[index_j];
		}

		ErrorAndParameters<Matd, Matd, Matd> error_and_parameters_inverse = computeErrorAndParameters(index_i, -dt);
		updateStates(index_i, -dt, error_and_parameters_inverse);

		ErrorAndParameters<Matd, Matd, Matd> error_and_parameters_inverse_after = computeErrorAndParameters(index_i, -dt);
		if (error_and_parameters_inverse_after.error_.norm() > error_and_parameters_inverse.error_.norm())
		{
			B_[index_i] = B_backup_[index_i];
			Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				B_[index_j] = B_backup_[index_j];
			}
		}
		else
		{
			evolution_indicator_[index_i] = 2;
		}
	}
	else
	{
		evolution_indicator_[index_i] = 1;
	}
}
//=================================================================================================//
CorrectionMatrixRegularization::CorrectionMatrixRegularization(BaseInnerRelation& inner_relation, Real eta)
	: LocalDynamics(inner_relation.getSPHBody()), RelaxDataDelegateInner(inner_relation),
	eta_(eta), mass_(particles_->mass_), Vol_(particles_->Vol_),
	B_(*particles_->template getVariableByName<Matd>("CorrectionMatrix"))
{
	particles_->registerVariable(correction_matrix_variation_, "correction_matrix_variation_");
	particles_->addVariableToWrite<Real>("correction_matrix_variation_");
}
//=================================================================================================//
ErrorAndParameters<Matd, Matd, Matd> CorrectionMatrixRegularization::
computeErrorAndParameters(size_t index_i, Real dt)
{
	ErrorAndParameters<Matd, Matd, Matd> error_and_parameters;
	Neighborhood& inner_neighborhood = this->inner_configuration_[index_i];
	for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
	{
		size_t& index_j = inner_neighborhood.j_[n];

		Matd variable_derivative = B_[index_i] - B_[index_j];

		Matd parameter_b = 2.0 * eta_ * inner_neighborhood.e_ij_[n] * inner_neighborhood.e_ij_[n].transpose() *
			inner_neighborhood.dW_ijV_j_[n] * Vol_[index_i] * dt / inner_neighborhood.r_ij_[n];

		error_and_parameters.error_ -= variable_derivative * parameter_b;
		error_and_parameters.a_ += parameter_b;
		error_and_parameters.c_ += parameter_b * parameter_b;
	}
	error_and_parameters.a_ -= mass_[index_i] * Matd::Identity();
	return error_and_parameters;
}
//=================================================================================================//
void CorrectionMatrixRegularization::
updateStates(size_t index_i, Real dt, const ErrorAndParameters<Matd, Matd, Matd>& error_and_parameters)
{
	Matd parameter_l = error_and_parameters.a_ * error_and_parameters.a_ + error_and_parameters.c_;
	Matd parameter_k = parameter_l.inverse() * error_and_parameters.error_;

	parameter_k = 0.1 * parameter_k;

	B_[index_i] += error_and_parameters.a_ * parameter_k;

	Neighborhood& inner_neighborhood = inner_configuration_[index_i];
	for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
	{
		size_t index_j = inner_neighborhood.j_[n];
		Matd parameter_b = 2.0 * eta_ * inner_neighborhood.e_ij_[n] * inner_neighborhood.e_ij_[n].transpose() *
			inner_neighborhood.dW_ijV_j_[n] * Vol_[index_i] * dt / inner_neighborhood.r_ij_[n];

		B_[index_j] -= parameter_b * parameter_k;
	}
}
//=================================================================================================//
void CorrectionMatrixRegularization::interaction(size_t index_i, Real dt)
{
	ErrorAndParameters<Matd, Matd, Matd> error_and_parameters = computeErrorAndParameters(index_i, dt);
	updateStates(index_i, dt, error_and_parameters);
	correction_matrix_variation_[index_i] = error_and_parameters.error_.norm();
}
//=================================================================================================//
FirstOrderEvolutionStep::FirstOrderEvolutionStep(BaseInnerRelation& inner_relation, bool level_set_correction)
	: BaseDynamics<void>(inner_relation.getSPHBody()), time_step_size_(0.01),
	get_time_step_(*inner_relation.real_body_)
{
	first_order_consistency_evolution_ =
		makeUnique<InteractionSplit<FirstOrderConsistencyEvolution>>(inner_relation, level_set_correction);
}
//=================================================================================================//
void FirstOrderEvolutionStep::exec(Real dt)
{
	first_order_consistency_evolution_->exec(dt);
}
//=================================================================================================//
UpdateParticleKineticEnergy::
UpdateParticleKineticEnergy(BaseInnerRelation& inner_relation) :
	LocalDynamics(inner_relation.getSPHBody()), GeneralDataDelegateInner(inner_relation),
	mass_(particles_->mass_), acc_(particles_->acc_)
{
	particles_->registerVariable(particle_kinetic_energy, "particle_kinetic_energy");
	particles_->addVariableToWrite<Real>("particle_kinetic_energy");
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
	level_set_shape_ = DynamicCast<LevelSetShape>(this, sph_body_.body_shape_);
	particles_->registerVariable(corrected_zero_order_error_, "corrected_zero_order_error");
	particles_->addVariableToWrite<Real>("corrected_zero_order_error");
	particles_->registerVariable(corrected_zero_order_, "corrected_zero_order");
	particles_->addVariableToWrite<Vecd>("corrected_zero_order");
}
//=================================================================================================//
void CheckCorrectedZeroOrderConsistency::interaction(size_t index_i, Real dt)
{
	Vecd acceleration = Vecd::Zero();
	Neighborhood& inner_neighborhood = inner_configuration_[index_i];
	for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
	{
		size_t index_j = inner_neighborhood.j_[n];
		acceleration -= 0.5 * (B_[index_i] + B_[index_j]) * inner_neighborhood.e_ij_[n] * inner_neighborhood.dW_ijV_j_[n];
	}

	/*if (level_set_correction_)
	{
		acceleration -= 0.5 * (B_[index_i]  + B_[index_i]) * level_set_shape_->computeKernelGradientIntegral(
						pos_[index_i], sph_adaptation_->SmoothingLengthRatio(index_i));
	}*/

	corrected_zero_order_[index_i] = acceleration;
	corrected_zero_order_error_[index_i] = acceleration.norm();
}
//=================================================================================================//
CheckCorrectedFirstOrderConsistency::
CheckCorrectedFirstOrderConsistency(BaseInnerRelation& inner_relation, bool level_set_correction) :
	LocalDynamics(inner_relation.getSPHBody()), GeneralDataDelegateInner(inner_relation),
	level_set_correction_(level_set_correction), pos_(particles_->pos_),
	B_(*particles_->template getVariableByName<Matd>("CorrectionMatrix")),
	sph_adaptation_(sph_body_.sph_adaptation_)
{
	level_set_shape_ = DynamicCast<LevelSetShape>(this, sph_body_.body_shape_);
	particles_->registerVariable(corrected_first_order_error_, "corrected_first_order_error");
	particles_->addVariableToWrite<Real>("corrected_first_order_error");
	particles_->registerVariable(corrected_first_order_, "corrected_first_order");
	particles_->addVariableToWrite<Matd>("corrected_first_order");
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

	corrected_first_order_[index_i] = acceleration;
	corrected_first_order_error_[index_i] = (acceleration - Matd::Identity()).norm();
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
