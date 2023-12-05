#include "relax_dynamics.hpp"
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
UpdateParticleKineticEnergy::
UpdateParticleKineticEnergy(BaseInnerRelation& inner_relation) :
	LocalDynamics(inner_relation.getSPHBody()), GeneralDataDelegateInner(inner_relation),
	mass_(particles_->mass_), acc_(particles_->acc_)
{
	particles_->registerVariable(particle_kinetic_energy, "ParticleKineticEnergy");
	particles_->addVariableToWrite<Real>("ParticleKineticEnergy");
};
//=================================================================================================//
void UpdateParticleKineticEnergy::update(size_t index_i, Real dt)
{
	particle_kinetic_energy[index_i] = acc_[index_i].norm(); /* L2 norm. */
};
//=================================================================================================//
CheckCorrectedZeroOrderConsistency::
CheckCorrectedZeroOrderConsistency(BaseInnerRelation& inner_relation, bool level_set_correction) :
	LocalDynamics(inner_relation.getSPHBody()), GeneralDataDelegateInner(inner_relation),
	level_set_correction_(level_set_correction), pos_(particles_->pos_),
	B_(*particles_->template getVariableByName<Matd>("KernelCorrectionMatrix")),
	sph_adaptation_(sph_body_.sph_adaptation_),
    constrained_distance_(0.5 * sph_body_.sph_adaptation_->MinimumSpacing())
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
		acceleration -= 0.5 * (B_[index_i] + B_[index_j]) * inner_neighborhood.e_ij_[n] * inner_neighborhood.dW_ijV_j_[n];
	}

	if (level_set_correction_)
	{
        Real phi = level_set_shape_->findSignedDistance(pos_[index_i]);
        Real overlap = level_set_shape_->computeKernelIntegral(pos_[index_i], 
                       sph_adaptation_->SmoothingLengthRatio(index_i));

        if (phi > -constrained_distance_)
        {
            acceleration -= 0.5 * (B_[index_i] + B_[index_i]) * level_set_shape_->computeKernelGradientIntegral(
                pos_[index_i], sph_adaptation_->SmoothingLengthRatio(index_i)) * (1 + overlap);
        }
        else
        {
            acceleration -= 0.5 * (B_[index_i]  + B_[index_i]) * level_set_shape_->computeKernelGradientIntegral(
						     pos_[index_i], sph_adaptation_->SmoothingLengthRatio(index_i)) * (1 + overlap);
        };
	}

	corrected_zero_order_error_[index_i] = acceleration;
	corrected_zero_order_error_norm_[index_i] = acceleration.norm();
}
//=================================================================================================//
CheckCorrectedFirstOrderConsistency::
CheckCorrectedFirstOrderConsistency(BaseInnerRelation& inner_relation, bool level_set_correction) :
	LocalDynamics(inner_relation.getSPHBody()), GeneralDataDelegateInner(inner_relation),
	level_set_correction_(level_set_correction), pos_(particles_->pos_),
	B_(*particles_->template getVariableByName<Matd>("KernelCorrectionMatrix")),
	sph_adaptation_(sph_body_.sph_adaptation_),
     constrained_distance_(0.5 * sph_body_.sph_adaptation_->MinimumSpacing())
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

	if (level_set_correction_)
	{
        Real phi = level_set_shape_->findSignedDistance(pos_[index_i]);
        Real overlap = level_set_shape_->computeKernelIntegral(pos_[index_i], 
                       sph_adaptation_->SmoothingLengthRatio(index_i));

        if (phi > -constrained_distance_)
        {
            acceleration -= 0.5 * (B_[index_i] + B_[index_i]) * level_set_shape_->computeDisplacementKernelGradientIntegral(
						pos_[index_i], sph_adaptation_->SmoothingLengthRatio(index_i)) * (1);
        }
        else
        {
            acceleration -= 0.5 * (B_[index_i] + B_[index_i]) * level_set_shape_->computeDisplacementKernelGradientIntegral(
						pos_[index_i], sph_adaptation_->SmoothingLengthRatio(index_i)) * (1);
        };
	}

	corrected_first_order_error_[index_i] = acceleration;
	corrected_first_order_error_norm_[index_i] = (acceleration - Matd::Identity()).norm();
}
//=================================================================================================//
CheckL2NormError::
CheckL2NormError(BaseInnerRelation& inner_relation, bool level_set_correction) :
    LocalDynamics(inner_relation.getSPHBody()), GeneralDataDelegateInner(inner_relation),
    level_set_correction_(level_set_correction), 
    scalar_(*particles_->template getVariableByName<Real>("Scalar")), 
    pos_(particles_->pos_),
    B_(*particles_->registerSharedVariable<Matd>("KernelCorrectionMatrix")),
    sph_adaptation_(sph_body_.sph_adaptation_),
    constrained_distance_(4 * sph_body_.sph_adaptation_->MinimumSpacing())
{
    particles_->addVariableToWrite<Real>("Scalar");
    particles_->registerVariable(average_label_, "AverageLabel");
    particles_->addVariableToWrite<Real>("AverageLabel");
    particles_->registerVariable(analytical_, "Analytical");
    particles_->addVariableToWrite<Real>("Analytical");
    particles_->registerVariable(L2_NKGC_, "NKGCError");
    particles_->addVariableToWrite<Real>("NKGCError");
    particles_->registerVariable(L2_SKGC_, "SKGCError");
    particles_->addVariableToWrite<Real>("SKGCError");
    particles_->registerVariable(L2_CKGC_, "CKGCError");
    particles_->addVariableToWrite<Real>("CKGCError");
    level_set_shape_ = DynamicCast<LevelSetShape>(this, sph_body_.body_shape_);
}
//=================================================================================================//
void CheckL2NormError::interaction(size_t index_i, Real dt)
{
    Vecd NKGC_error_ = Vecd::Zero();
    Vecd SKGC_error_ = Vecd::Zero();
    Vecd CKGC_error_ = Vecd::Zero();
     
    Neighborhood& inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Vecd& e_ij = inner_neighborhood.e_ij_[n];
        Real dW_ijV_j = inner_neighborhood.dW_ijV_j_[n];

        NKGC_error_ += (scalar_[index_i] + scalar_[index_j]) * e_ij * dW_ijV_j;
        SKGC_error_ += (scalar_[index_i] * B_[index_i] + scalar_[index_j] * B_[index_j]) * e_ij * dW_ijV_j;
        CKGC_error_ += (scalar_[index_i] * B_[index_j] + scalar_[index_j] * B_[index_i]) * e_ij * dW_ijV_j;
    }

    analytical_[index_i] = 2 * Pi * cos(2 * Pi * pos_[index_i][0]);
 
    L2_NKGC_[index_i] = (NKGC_error_ - Vec2d(2 * Pi * cos(2 * Pi * pos_[index_i][0]), 0)).norm();
    L2_SKGC_[index_i] = (SKGC_error_ - Vec2d(2 * Pi * cos(2 * Pi * pos_[index_i][0]), 0)).norm();
    L2_CKGC_[index_i] = (CKGC_error_ - Vec2d(2 * Pi * cos(2 * Pi * pos_[index_i][0]), 0)).norm();
};
//=================================================================================================//
CheckConsistencyRealization::
CheckConsistencyRealization(BaseInnerRelation& inner_relation, bool level_set_correction) :
    LocalDynamics(inner_relation.getSPHBody()), GeneralDataDelegateInner(inner_relation),
    level_set_correction_(level_set_correction), pos_(particles_->pos_),
    scalar_(*particles_->template getVariableByName<Real>("Scalar")),
    vector_(*particles_->template getVariableByName<Vecd>("Vector")),
    matrix_(*particles_->template getVariableByName<Matd>("Matrix")),
    B_(*particles_->registerSharedVariable<Matd>("KernelCorrectionMatrix")),
    sph_adaptation_(sph_body_.sph_adaptation_)
{
    particles_->registerVariable(reproduce_scalar_gradient_, "ReproduceScalarGradient");
    particles_->addVariableToWrite<Vecd>("ReproduceScalarGradient");
    particles_->registerVariable(reproduce_vector_gradient_, "ReproduceVectorGradient");
    particles_->addVariableToWrite<Real>("ReproduceVectorGradient");
    particles_->registerVariable(reproduce_matrix_gradient_, "ReproduceMatrixGradient");
    particles_->addVariableToWrite<Vecd>("ReproduceMatrixGradient");

    particles_->registerVariable(ACterm_norm_, "ACTERMNORM");
    particles_->addVariableToWrite<Real>("ACTERMNORM");
    particles_->registerVariable(ACterm_, "ACTERM");
    particles_->addVariableToWrite<Vecd>("ACTERM");
    particles_->registerVariable(ASterm_norm_, "ASTERMNORM");
    particles_->addVariableToWrite<Real>("ASTERMNORM");
    particles_->registerVariable(ASterm_, "ASTERM");
    particles_->addVariableToWrite<Vecd>("ASTERM");
    particles_->registerVariable(Cterm_norm_, "CTERMNORM");
    particles_->addVariableToWrite<Real>("CTERMNORM");
    particles_->registerVariable(Cterm_, "CTERM");
    particles_->addVariableToWrite<Vecd>("CTERM");
    particles_->registerVariable(NKGC_norm_, "NKGCNORM");
    particles_->addVariableToWrite<Real>("NKGCNORM");
    particles_->registerVariable(NKGC_, "NKGC");
    particles_->addVariableToWrite<Vecd>("NKGC");

    particles_->addVariableToWrite<Real>("Scalar");
    particles_->addVariableToWrite<Vecd>("Vector");
    particles_->addVariableToWrite<Matd>("Matrix");

    level_set_shape_ = DynamicCast<LevelSetShape>(this, sph_body_.body_shape_);
}
//=================================================================================================//
void CheckConsistencyRealization::interaction(size_t index_i, Real dt)
{
    Vecd AC = Vecd::Zero();
    Vecd C = Vecd::Zero();
    Vecd AS = Vecd::Zero();
    Vecd NKGC = Vecd::Zero();
    
    Neighborhood& inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Vecd &e_ij = inner_neighborhood.e_ij_[n];
        Real dW_ijV_j = inner_neighborhood.dW_ijV_j_[n];

        AC += (scalar_[index_i] * B_[index_i] + scalar_[index_j] * B_[index_j]) * e_ij * dW_ijV_j;
        AS += scalar_[index_i] * (B_[index_i] + B_[index_j]) * e_ij * dW_ijV_j;
        C -= (scalar_[index_i] - scalar_[index_j]) * B_[index_j] * e_ij * dW_ijV_j;
        NKGC += (scalar_[index_i] + scalar_[index_j]) * e_ij * dW_ijV_j;
    }
    
    ACterm_[index_i] = AC;
    ACterm_norm_[index_i] = (AC - Vec2d(2 * Pi * cos(2 * Pi * pos_[index_i][0]), 0)).norm();
    ASterm_[index_i] = AS;
    ASterm_norm_[index_i] = AS.norm();
    Cterm_[index_i] = C;
    Cterm_norm_[index_i] = (C - Vec2d(2 * Pi * cos(2 * Pi * pos_[index_i][0]), 0)).norm();
    NKGC_[index_i] = NKGC;
    NKGC_norm_[index_i] = NKGC.norm();
}
//=================================================================================================//
CheckReverseConsistencyRealization::
CheckReverseConsistencyRealization(BaseInnerRelation& inner_relation, bool level_set_correction) :
    LocalDynamics(inner_relation.getSPHBody()), GeneralDataDelegateInner(inner_relation),
    level_set_correction_(level_set_correction), pos_(particles_->pos_),
    scalar_(*particles_->template getVariableByName<Real>("Scalar")),
    B_(*particles_->template getVariableByName<Matd>("KernelCorrectionMatrix")),
    sph_adaptation_(sph_body_.sph_adaptation_)
{
    particles_->registerVariable(ABterm_norm_, "ABTERMNORM");
    particles_->addVariableToWrite<Real>("ABTERMNORM");
    particles_->registerVariable(ABterm_, "ABTERM");
    particles_->addVariableToWrite<Vecd>("ABTERM");
    particles_->registerVariable(ARterm_norm_, "ARTERMNORM");
    particles_->addVariableToWrite<Real>("ARTERMNORM");
    particles_->registerVariable(ARterm_, "ARTERM");
    particles_->addVariableToWrite<Vecd>("ARTERM");
    particles_->registerVariable(Bterm_norm_, "BTERMNORM");
    particles_->addVariableToWrite<Real>("BTERMNORM");
    particles_->registerVariable(Bterm_, "BTERM");
    particles_->addVariableToWrite<Vecd>("BTERM");
    level_set_shape_ = DynamicCast<LevelSetShape>(this, sph_body_.body_shape_);
}
//=================================================================================================//
void CheckReverseConsistencyRealization::interaction(size_t index_i, Real dt)
{
    Neighborhood& inner_neighborhood = inner_configuration_[index_i];
    Vecd AB = Vecd::Zero();
    Vecd B = Vecd::Zero();
    Vecd AR = Vecd::Zero();
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Vecd& e_ij = inner_neighborhood.e_ij_[n];
        Real dW_ijV_j = inner_neighborhood.dW_ijV_j_[n];

        AB += (scalar_[index_i] * B_[index_j] + scalar_[index_j] * B_[index_i]) * e_ij * dW_ijV_j;
        AR += scalar_[index_i] * (B_[index_i] + B_[index_j]) * e_ij * dW_ijV_j;
        B -= (scalar_[index_i] - scalar_[index_j]) * B_[index_i] * e_ij * dW_ijV_j;
    }

    ABterm_[index_i] = AB;
    ABterm_norm_[index_i] = (AB - Vec2d(2 * Pi * cos(2 * Pi * pos_[index_i][0]), 0)).norm();
    ARterm_[index_i] = AR;
    ARterm_norm_[index_i] = AR.norm();
    Bterm_[index_i] = B;
    Bterm_norm_[index_i] = (B - Vec2d(2 * Pi * cos(2 * Pi * pos_[index_i][0]), 0)).norm();
}
//=================================================================================================//
ShellMidSurfaceBounding::
    ShellMidSurfaceBounding(NearShapeSurface &body_part, BaseInnerRelation &inner_relation)
    : BaseLocalDynamics<BodyPartByCell>(body_part), RelaxDataDelegateInner(inner_relation),
      pos_(particles_->pos_), constrained_distance_(0.5 * sph_body_.sph_adaptation_->MinimumSpacing()),
      particle_spacing_ref_(sph_body_.sph_adaptation_->MinimumSpacing()),
      level_set_shape_(DynamicCast<LevelSetShape>(this, sph_body_.body_shape_)) {}
//=================================================================================================//
void ShellMidSurfaceBounding::update(size_t index_i, Real dt)
{
    Vecd none_normalized_normal = level_set_shape_->findLevelSetGradient(pos_[index_i]);
    Vecd normal = level_set_shape_->findNormalDirection(pos_[index_i]);
    Real factor = 0.2 * none_normalized_normal.norm();
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
    ShellRelaxationStepInner(BaseInnerRelation &inner_relation, bool level_set_correction)
    : RelaxationStepInner(inner_relation, level_set_correction),
      update_shell_particle_position_(*real_body_),
      mid_surface_bounding_(near_shape_surface_, inner_relation) {}
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
