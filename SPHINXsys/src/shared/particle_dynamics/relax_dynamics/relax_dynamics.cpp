#include "relax_dynamics.h"
#include "base_particles.hpp"
#include "particle_generator_lattice.h"
#include "level_set_shape.h"

//========================================================================================================//
namespace SPH
{
	//=====================================================================================================//
	namespace relax_dynamics
	{
		//=================================================================================================//
		GetTimeStepSizeSquare::GetTimeStepSizeSquare(SPHBody &sph_body)
			: LocalDynamicsReduce<Real, ReduceMax>(sph_body, Real(0))
			, RelaxDataDelegateSimple(sph_body)
			, acc_(particles_->acc_)
			, h_ref_(sph_body.sph_adaptation_->ReferenceSmoothingLength()) 
		{}
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
			: LocalDynamics(inner_relation.sph_body_), RelaxDataDelegateInner(inner_relation),
			  acc_(particles_->acc_), pos_(particles_->pos_) {}
		//=================================================================================================//
		void RelaxationAccelerationInner::interaction(size_t index_i, Real dt)
		{
			Vecd acceleration = Vecd::Zero();
			const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				acceleration -= 2.0 * inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
			}
			acc_[index_i] = acceleration;
		}
		//=================================================================================================//
		RelaxationAccelerationInnerWithLevelSetCorrection::
			RelaxationAccelerationInnerWithLevelSetCorrection(BaseInnerRelation &inner_relation)
			: RelaxationAccelerationInner(inner_relation), sph_adaptation_(sph_body_.sph_adaptation_)
		{
			level_set_shape_ = DynamicCast<LevelSetShape>(this, sph_body_.body_shape_);
		}
		//=================================================================================================//
		void RelaxationAccelerationInnerWithLevelSetCorrection::interaction(size_t index_i, Real dt)
		{
			RelaxationAccelerationInner::interaction(index_i, dt);
			acc_[index_i] -= 2.0 * level_set_shape_->computeKernelGradientIntegral(
									   pos_[index_i], sph_adaptation_->SmoothingLengthRatio(index_i));
		}
		//=================================================================================================//
		UpdateParticlePosition::UpdateParticlePosition(SPHBody &sph_body)
			: LocalDynamics(sph_body)
			, RelaxDataDelegateSimple(sph_body)
			, sph_adaptation_(sph_body.sph_adaptation_)
			, pos_(particles_->pos_)
			, acc_(particles_->acc_) 
		{}
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
			Vol_[index_i] = powerN(local_spacing, Dimensions);
		}
		//=================================================================================================//
		RelaxationAccelerationComplex::
			RelaxationAccelerationComplex(ComplexRelation &complex_relation)
			: LocalDynamics(complex_relation.sph_body_),
			  RelaxDataDelegateComplex(complex_relation),
			  acc_(particles_->acc_), pos_(particles_->pos_) {}
		//=================================================================================================//
		void RelaxationAccelerationComplex::interaction(size_t index_i, Real dt)
		{
			Vecd acceleration = Vecd::Zero();
			Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				acceleration -= 2.0 * inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
			}

			/** Contact interaction. */
			for (size_t k = 0; k < contact_configuration_.size(); ++k)
			{
				Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					size_t index_j = contact_neighborhood.j_[n];

					acceleration -= 2.0 * contact_neighborhood.dW_ijV_j_[n] * contact_neighborhood.e_ij_[n];
				}
			}

			acc_[index_i] = acceleration;
		}
		//=================================================================================================//
		ShapeSurfaceBounding::ShapeSurfaceBounding(NearShapeSurface &near_shape_surface)
			: LocalDynamics(near_shape_surface.getSPHBody())
			, RelaxDataDelegateSimple(sph_body_)
			, pos_(particles_->pos_)
			, constrained_distance_(0.5 * sph_body_.sph_adaptation_->MinimumSpacing())
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
		RelaxationStepInner::
			RelaxationStepInner(BaseInnerRelation &inner_relation, bool level_set_correction)
			: BaseDynamics<void>(), real_body_(inner_relation.real_body_),
			  inner_relation_(inner_relation), near_shape_surface_(*real_body_),
			  get_time_step_square_(*real_body_), update_particle_position_(*real_body_),
			  surface_bounding_(near_shape_surface_)
		{
			if (!level_set_correction)
			{
				relaxation_acceleration_inner_ =
					std::move(makeUnique<InteractionDynamics<RelaxationAccelerationInner>>(inner_relation));
			}
			else
			{
				relaxation_acceleration_inner_ =
					std::move(makeUnique<InteractionDynamics<RelaxationAccelerationInnerWithLevelSetCorrection>>(inner_relation));
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
		void RelaxationStepInner::parallel_exec(Real dt)
		{
			real_body_->updateCellLinkedList();
			inner_relation_.updateConfiguration();
			relaxation_acceleration_inner_->parallel_exec();
			Real dt_square = get_time_step_square_.parallel_exec();
			update_particle_position_.parallel_exec(dt_square);
			surface_bounding_.parallel_exec();
		}
		//=================================================================================================//
		RelaxationAccelerationComplexWithLevelSetCorrection::
			RelaxationAccelerationComplexWithLevelSetCorrection(ComplexRelation &body_complex_relation, const std::string &shape_name)
			: RelaxationAccelerationComplex(body_complex_relation),
			  sph_adaptation_(sph_body_.sph_adaptation_)
		{
			ComplexShape &complex_shape = DynamicCast<ComplexShape>(this, *sph_body_.body_shape_);
			level_set_shape_ = DynamicCast<LevelSetShape>(this, complex_shape.getShapeByName(shape_name));
		}
		//=================================================================================================//
		void RelaxationAccelerationComplexWithLevelSetCorrection::interaction(size_t index_i, Real dt)
		{
			RelaxationAccelerationComplex::interaction(index_i, dt);

			acc_[index_i] -= 2.0 * level_set_shape_->computeKernelGradientIntegral(
									   pos_[index_i], sph_adaptation_->SmoothingLengthRatio(index_i));
		}
		//=================================================================================================//
		RelaxationStepComplex::RelaxationStepComplex(ComplexRelation &body_complex_relation,
													 const std::string &shape_name, bool level_set_correction)
			: BaseDynamics<void>(),
			  real_body_(body_complex_relation.inner_relation_.real_body_),
			  complex_relation_(body_complex_relation),
			  near_shape_surface_(*real_body_, shape_name),
			  get_time_step_square_(*real_body_), update_particle_position_(*real_body_),
			  surface_bounding_(near_shape_surface_)
		{
			if (!level_set_correction)
			{
				relaxation_acceleration_complex_ =
					std::move(makeUnique<InteractionDynamics<RelaxationAccelerationComplex>>(body_complex_relation));
			}
			else
			{
				relaxation_acceleration_complex_ =
					std::move(makeUnique<InteractionDynamics<RelaxationAccelerationComplexWithLevelSetCorrection>>(body_complex_relation, shape_name));
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
		void RelaxationStepComplex::parallel_exec(Real dt)
		{
			real_body_->updateCellLinkedList();
			complex_relation_.updateConfiguration();
			relaxation_acceleration_complex_->parallel_exec();
			Real dt_square = get_time_step_square_.parallel_exec();
			update_particle_position_.parallel_exec(dt_square);
			surface_bounding_.parallel_exec();
		}
		//=================================================================================================//
		ShellMidSurfaceBounding::
			ShellMidSurfaceBounding(NearShapeSurface &body_part, BaseInnerRelation &inner_relation,
									Real thickness, Real level_set_refinement_ratio)
			: LocalDynamics(body_part.getSPHBody()), RelaxDataDelegateInner(inner_relation),
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
			: BaseDynamics<void>(),
			  convergence_criterion_(cos(0.01 * Pi)),
			  consistency_criterion_(consistency_criterion),
			  normal_prediction_(inner_relation.sph_body_, thickness),
			  normal_prediction_convergence_check_(inner_relation.sph_body_, convergence_criterion_),
			  consistency_correction_(inner_relation, consistency_criterion_),
			  consistency_updated_check_(inner_relation.sph_body_),
			  smoothing_normal_(inner_relation) {}
		//=================================================================================================//
		void ShellNormalDirectionPrediction::exec(Real dt)
		{
			predictNormalDirection();
			correctNormalDirection();
			predictNormalDirection();
			smoothing_normal_.parallel_exec();
		}
		//=================================================================================================//
		void ShellNormalDirectionPrediction::predictNormalDirection()
		{
			bool prediction_convergence = false;
			size_t ite_predict = 0;
			while (!prediction_convergence)
			{
				normal_prediction_.parallel_exec();
				prediction_convergence = normal_prediction_convergence_check_.parallel_exec();
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
				consistency_correction_.parallel_exec();
				consistency_updated = consistency_updated_check_.parallel_exec();
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
			: RelaxDataDelegateSimple(sph_body)
			, LocalDynamics(sph_body), thickness_(thickness)
			, level_set_shape_(DynamicCast<LevelSetShape>(this, sph_body.body_shape_))
			, pos_(particles_->pos_)
			, n_(*particles_->getVariableByName<Vecd>("NormalDirection"))
		{
			particles_->registerVariable(n_temp_, "PreviousNormalDirection", [&](size_t i) -> Vecd { return n_[i]; });
		}
		//=================================================================================================//
		void ShellNormalDirectionPrediction::NormalPrediction::update(size_t index_i, Real dt)
		{
			n_temp_[index_i] = n_[index_i];
			n_[index_i] = level_set_shape_->findNormalDirection(pos_[index_i] + 0.3 * thickness_ * n_temp_[index_i]);
		}
		//=================================================================================================//
		ShellNormalDirectionPrediction::PredictionConvergenceCheck::PredictionConvergenceCheck(SPHBody &sph_body, Real convergence_criterion)
			: LocalDynamicsReduce<bool, ReduceAND>(sph_body, true)
			, RelaxDataDelegateSimple(sph_body), convergence_criterion_(convergence_criterion)
			, n_(*particles_->getVariableByName<Vecd>("NormalDirection"))
			, n_temp_(*particles_->getVariableByName<Vecd>("PreviousNormalDirection")) 
		{}
		//=================================================================================================//
		bool ShellNormalDirectionPrediction::PredictionConvergenceCheck::reduce(size_t index_i, Real dt)
		{
			return n_[index_i].dot(n_temp_[index_i]) > convergence_criterion_;
		}
		//=================================================================================================//
		ShellNormalDirectionPrediction::ConsistencyCorrection::
			ConsistencyCorrection(BaseInnerRelation &inner_relation, Real consistency_criterion)
			: LocalDynamics(inner_relation.sph_body_), RelaxDataDelegateInner(inner_relation),
			  consistency_criterion_(consistency_criterion),
			  n_(*particles_->getVariableByName<Vecd>("NormalDirection"))
		{
			particles_->registerVariable(updated_indicator_, "UpdatedIndicator", [&](size_t i) -> int {return 0;});
			updated_indicator_[particles_->total_real_particles_ / 3] = 1;
		}
		//=================================================================================================//
		void ShellNormalDirectionPrediction::ConsistencyCorrection::interaction(size_t index_i, Real dt)
		{
			mutex_modify_neighbor_.lock();
			const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				if (updated_indicator_[index_i] == 1)
				{
					size_t index_j = inner_neighborhood.j_[n];
					if (updated_indicator_[index_j] == 0)
					{
						updated_indicator_[index_j] = 1;
						if (n_[index_i].dot(n_[index_j]) < consistency_criterion_)
						{
							if (n_[index_i].dot(-n_[index_j]) < consistency_criterion_)
							{
								n_[index_j] = n_[index_i];
								updated_indicator_[index_j] = 2;
							}
							else
							{
								n_[index_j] = -n_[index_j];
								updated_indicator_[index_j] = 1;
							}
						}
					}
				}
			}
			mutex_modify_neighbor_.unlock();
		}
		//=================================================================================================//
		ShellNormalDirectionPrediction::ConsistencyUpdatedCheck::ConsistencyUpdatedCheck(SPHBody &sph_body)
			: LocalDynamicsReduce<bool, ReduceAND>(sph_body, true)
			,  RelaxDataDelegateSimple(sph_body)
			,  updated_indicator_(*particles_->getVariableByName<int>("UpdatedIndicator")) 
		{}
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
		void ShellRelaxationStepInner::parallel_exec(Real ite_p)
		{
			real_body_->updateCellLinkedList();
			inner_relation_.updateConfiguration();
			relaxation_acceleration_inner_->parallel_exec();
			Real dt_square = get_time_step_square_.parallel_exec();
			update_shell_particle_position_.parallel_exec(dt_square);
			mid_surface_bounding_.parallel_exec();
		}
		//=================================================================================================//
	}
	//=================================================================================================//
}
//=================================================================================================//
