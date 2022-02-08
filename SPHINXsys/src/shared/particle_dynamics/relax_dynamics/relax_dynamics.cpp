/**
 * @file 	relax_dynamics.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 */

#include "relax_dynamics.h"
#include "particle_generator_lattice.h"
#include "level_set_shape.h"

using namespace SimTK;
//=================================================================================================//
namespace SPH
{
	//=================================================================================================//
	namespace relax_dynamics
	{
		//=================================================================================================//
		GetTimeStepSizeSquare::GetTimeStepSizeSquare(SPHBody &sph_body)
			: ParticleDynamicsReduce<Real, ReduceMax>(sph_body),
			  RelaxDataDelegateSimple(sph_body), dvel_dt_(particles_->dvel_dt_),
			  h_ref_(sph_body.sph_adaptation_->ReferenceSmoothingLength())
		{
			//given by sound criteria by assuming unit speed of sound and zero particle velocity
			initial_reference_ = 1.0 / h_ref_;
		}
		//=================================================================================================//
		Real GetTimeStepSizeSquare::ReduceFunction(size_t index_i, Real dt)
		{
			return dvel_dt_[index_i].norm();
		}
		//=================================================================================================//
		Real GetTimeStepSizeSquare::OutputResult(Real reduced_value)
		{
			return 0.0625 * h_ref_ / (reduced_value + TinyReal);
		}
		//=================================================================================================//
		RelaxationAccelerationInner::RelaxationAccelerationInner(BaseBodyRelationInner &inner_relation)
			: InteractionDynamics(*inner_relation.sph_body_),
			  RelaxDataDelegateInner(inner_relation),
			  Vol_(particles_->Vol_), dvel_dt_(particles_->dvel_dt_), pos_n_(particles_->pos_n_) {}
		//=================================================================================================//
		void RelaxationAccelerationInner::Interaction(size_t index_i, Real dt)
		{
			Vecd acceleration(0);
			const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				acceleration -= 2.0 * inner_neighborhood.dW_ij_[n] * inner_neighborhood.e_ij_[n] * Vol_[index_j];
			}
			dvel_dt_[index_i] = acceleration;
		}
		//=================================================================================================//
		RelaxationAccelerationInnerWithLevelSetCorrection::
			RelaxationAccelerationInnerWithLevelSetCorrection(BaseBodyRelationInner &inner_relation)
			: RelaxationAccelerationInner(inner_relation)
		{
			level_set_shape_ = DynamicCast<LevelSetShape>(this, body_->body_shape_.getShapeByName(body_->getBodyName()));
			if (level_set_shape_ == nullptr)
			{
				std::cout << "\n FAILURE: LevelSetShape is undefined!" << std::endl;
				std::cout << __FILE__ << ':' << __LINE__ << std::endl;
				exit(1);
			}
		}
		//=================================================================================================//
		void RelaxationAccelerationInnerWithLevelSetCorrection::Interaction(size_t index_i, Real dt)
		{
			RelaxationAccelerationInner::Interaction(index_i, dt);
			dvel_dt_[index_i] -= 2.0 * level_set_shape_->computeKernelGradientIntegral(
										   pos_n_[index_i], sph_adaptation_->SmoothingLengthRatio(index_i));
		}
		//=================================================================================================//
		UpdateParticlePosition::UpdateParticlePosition(SPHBody &sph_body)
			: ParticleDynamicsSimple(sph_body), RelaxDataDelegateSimple(sph_body),
			  pos_n_(particles_->pos_n_), dvel_dt_(particles_->dvel_dt_) {}
		//=================================================================================================//
		void UpdateParticlePosition::Update(size_t index_i, Real dt_square)
		{
			pos_n_[index_i] += dvel_dt_[index_i] * dt_square * 0.5 / sph_adaptation_->SmoothingLengthRatio(index_i);
		}
		//=================================================================================================//
		UpdateSolidParticlePosition::UpdateSolidParticlePosition(SPHBody &sph_body)
			: ParticleDynamicsSimple(sph_body), solid_dynamics::SolidDataSimple(sph_body),
			  pos_0_(particles_->pos_0_), pos_n_(particles_->pos_n_), dvel_dt_(particles_->dvel_dt_) {}
		//=================================================================================================//
		void UpdateSolidParticlePosition::Update(size_t index_i, Real dt_square)
		{
			pos_n_[index_i] += dvel_dt_[index_i] * dt_square * 0.5 / sph_adaptation_->SmoothingLengthRatio(index_i);
			pos_0_[index_i] = pos_n_[index_i];
		}
		//=================================================================================================//
		UpdateSmoothingLengthRatioByBodyShape::UpdateSmoothingLengthRatioByBodyShape(SPHBody &sph_body)
			: ParticleDynamicsSimple(sph_body), RelaxDataDelegateSimple(sph_body),
			  h_ratio_(*particles_->getVariableByName<Real>("SmoothingLengthRatio")),
			  Vol_(particles_->Vol_), pos_n_(particles_->pos_n_),
			  body_shape_(sph_body.body_shape_), kernel_(*sph_body.sph_adaptation_->getKernel())
		{
			particle_spacing_by_body_shape_ =
				DynamicCast<ParticleSpacingByBodyShape>(this, sph_body.sph_adaptation_);
		}
		//=================================================================================================//
		void UpdateSmoothingLengthRatioByBodyShape::Update(size_t index_i, Real dt_square)
		{
			Real local_spacing = particle_spacing_by_body_shape_->getLocalSpacing(body_shape_, pos_n_[index_i]);
			h_ratio_[index_i] = sph_adaptation_->ReferenceSpacing() / local_spacing;
			Vol_[index_i] = powerN(local_spacing, Dimensions);
		}
		//=================================================================================================//
		RelaxationAccelerationComplex::
			RelaxationAccelerationComplex(ComplexBodyRelation &complex_relation)
			: InteractionDynamics(*complex_relation.sph_body_),
			  RelaxDataDelegateComplex(complex_relation),
			  Vol_(particles_->Vol_), dvel_dt_(particles_->dvel_dt_), pos_n_(particles_->pos_n_)
		{
			for (size_t k = 0; k != contact_particles_.size(); ++k)
			{
				contact_Vol_.push_back(&(contact_particles_[k]->Vol_));
			}
		}
		//=================================================================================================//
		void RelaxationAccelerationComplex::Interaction(size_t index_i, Real dt)
		{
			Vecd acceleration(0);
			Neighborhood &inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				acceleration -= 2.0 * inner_neighborhood.dW_ij_[n] * inner_neighborhood.e_ij_[n] * Vol_[index_j];
			}

			/** Contact interaction. */
			for (size_t k = 0; k < contact_configuration_.size(); ++k)
			{
				StdLargeVec<Real> &Vol_k = *(contact_Vol_[k]);
				Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					size_t index_j = contact_neighborhood.j_[n];

					acceleration -= 2.0 * Vol_k[index_j] * contact_neighborhood.dW_ij_[n] * contact_neighborhood.e_ij_[n];
				}
			}

			dvel_dt_[index_i] = acceleration;
		}
		//=================================================================================================//
		ShapeSurfaceBounding::
			ShapeSurfaceBounding(SPHBody &sph_body, NearShapeSurface &body_part)
			: PartDynamicsByCell(sph_body, body_part), RelaxDataDelegateSimple(sph_body),
			  pos_n_(particles_->pos_n_),
			  constrained_distance_(0.5 * sph_body.sph_adaptation_->MinimumSpacing())
		{
			level_set_shape_ = DynamicCast<LevelSetShape>(this, body_->body_shape_.getShapeByName(body_->getBodyName()));
		}
		//=================================================================================================//
		void ShapeSurfaceBounding::Update(size_t index_i, Real dt)
		{
			Real phi = level_set_shape_->findSignedDistance(pos_n_[index_i]);

			if (phi > -constrained_distance_)
			{
				Vecd unit_normal = level_set_shape_->findNormalDirection(pos_n_[index_i]);
				unit_normal /= unit_normal.norm() + TinyReal;
				pos_n_[index_i] -= (phi + constrained_distance_) * unit_normal;
			}
		}
		//=================================================================================================//
		ConstraintSurfaceParticles::
			ConstraintSurfaceParticles(SPHBody &sph_body, BodySurface &body_part)
			: PartSimpleDynamicsByParticle(sph_body, body_part), RelaxDataDelegateSimple(sph_body),
			  constrained_distance_(0.5 * sph_body.sph_adaptation_->MinimumSpacing()),
			  pos_n_(particles_->pos_n_)
		{
			level_set_shape_ = DynamicCast<LevelSetShape>(this, body_->body_shape_.getShapeByName(body_->getBodyName()));
		}
		//=================================================================================================//
		void ConstraintSurfaceParticles::Update(size_t index_i, Real dt)
		{
			Real phi = level_set_shape_->findSignedDistance(pos_n_[index_i]);
			Vecd unit_normal = level_set_shape_->findNormalDirection(pos_n_[index_i]);
			unit_normal /= unit_normal.norm() + TinyReal;
			pos_n_[index_i] -= (phi + constrained_distance_) * unit_normal;
		}
		//=================================================================================================//
		RelaxationStepInner::
			RelaxationStepInner(BaseBodyRelationInner &inner_relation, bool level_set_correction)
			: ParticleDynamics<void>(*inner_relation.sph_body_),
			  real_body_(inner_relation.real_body_), 
			inner_relation_(inner_relation),
			  near_shape_surface_(*real_body_),
			  get_time_step_square_(*real_body_), update_particle_position_(*real_body_),
			  surface_bounding_(*real_body_, near_shape_surface_),
			  relaxation_acceleration_inner_(
				  !level_set_correction
					  ? std::move(makeUnique<RelaxationAccelerationInner>(inner_relation))
					  : std::move(makeUnique<RelaxationAccelerationInnerWithLevelSetCorrection>(inner_relation))) {}
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
			RelaxationAccelerationComplexWithLevelSetCorrection(ComplexBodyRelation &body_complex_relation) :
			RelaxationAccelerationComplex(body_complex_relation)
		{
			level_set_complex_shape_ = DynamicCast<LevelSetShape>(this, body_->body_shape_.getShapeByName(body_->getBodyName()));
			if (level_set_complex_shape_ == nullptr)
			{
				std::cout << "\n FAILURE: LevelSetComplexShape is undefined!" << std::endl;
				std::cout << __FILE__ << ':' << __LINE__ << std::endl;
				exit(1);
			}
		}
		//=================================================================================================//
		void RelaxationAccelerationComplexWithLevelSetCorrection::Interaction(size_t index_i, Real dt)
		{
			RelaxationAccelerationComplex::Interaction(index_i, dt);

			for (size_t k = 0; k < contact_configuration_.size(); ++k)
			{
				StdLargeVec<Real>& Vol_k = *(contact_Vol_[k]);
				Neighborhood & contact_neighborhood = (*contact_configuration_[k])[index_i];
				if (contact_neighborhood.current_size_ == 0)
				{
					dvel_dt_[index_i] -= 2.0 * level_set_complex_shape_->
						computeKernelGradientIntegral(pos_n_[index_i], sph_adaptation_->SmoothingLengthRatio(index_i));
				}
			}
		}
		//=================================================================================================//
		RelaxationStepComplex::RelaxationStepComplex(ComplexBodyRelation &body_complex_relation, bool level_set_correction) :
			ParticleDynamics<void>(*body_complex_relation.sph_body_), 
			real_body_(body_complex_relation.inner_relation_.real_body_), 
			complex_relation_(body_complex_relation),
			near_shape_surface_(*real_body_),
			get_time_step_square_(*real_body_), update_particle_position_(*real_body_),
			surface_bounding_(*real_body_, near_shape_surface_),
			relaxation_acceleration_complex_(
				!level_set_correction
				? std::move(makeUnique<RelaxationAccelerationComplex>(body_complex_relation))
				: std::move(makeUnique<RelaxationAccelerationComplexWithLevelSetCorrection>(body_complex_relation))) {}
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
		void SolidRelaxationStepInner::exec(Real dt)
		{
			real_body_->updateCellLinkedList();
			inner_relation_.updateConfiguration();
			relaxation_acceleration_inner_->exec();
			Real dt_square = get_time_step_square_.exec();
			update_solid_particle_position_.exec(dt_square);
			surface_bounding_.exec();
		}
		//=================================================================================================//
		void SolidRelaxationStepInner::parallel_exec(Real dt)
		{
			real_body_->updateCellLinkedList();
			inner_relation_.updateConfiguration();
			relaxation_acceleration_inner_->parallel_exec();
			Real dt_square = get_time_step_square_.parallel_exec();
			update_solid_particle_position_.parallel_exec(dt_square);
			surface_bounding_.parallel_exec();

			// copy the updated position at the end of the relaxation step to avoid bugs
			std::copy(GetParticles()->pos_n_.begin(), GetParticles()->pos_n_.end(), GetParticles()->pos_0_.begin());
		}
		//=================================================================================================//
		ShellMidSurfaceBounding::
			ShellMidSurfaceBounding(SPHBody &sph_body, NearShapeSurface &body_part, BaseBodyRelationInner &inner_relation, 
				                    Real thickness, Real level_set_refinement_ratio) :
			PartDynamicsByCell(sph_body, body_part), RelaxDataDelegateInner(inner_relation),
			solid_particles_(dynamic_cast<SolidParticles*>(particles_)),
			pos_n_(solid_particles_->pos_n_), n_0_(solid_particles_->n_0_),
			constrained_distance_(0.5 * sph_body.sph_adaptation_->MinimumSpacing()),
			particle_spacing_ref_(sph_body.sph_adaptation_->MinimumSpacing()),
			thickness_(thickness), level_set_refinement_ratio_(level_set_refinement_ratio)
		{
			particles_->registerAVariable<Real>(color_, "Color");
			particles_->registerAVariable<Vecd>(temporary_n_0_, "TemporatyNormal");
			particles_->registerAVariable<Vecd>(previous_n_0_, "PreviousNormal");
			level_set_shape_ = DynamicCast<LevelSetShape>(this, body_->body_shape_.getShapeByName(body_->getBodyName()));
		}
		//=================================================================================================//
		void ShellMidSurfaceBounding::Update(size_t index_i, Real dt)
		{
			Vecd none_normalized_normal = level_set_shape_->findNoneNormalizedNormalDirection(pos_n_[index_i]);
			Vecd normal = none_normalized_normal / (none_normalized_normal.norm() + TinyReal);
			Real phi = level_set_shape_->findSignedDistance(pos_n_[index_i]);
			Real relaxation_factor_ratio = particle_spacing_ref_ / level_set_refinement_ratio_ / (0.1 * thickness_);
			Real factor = none_normalized_normal.norm() / relaxation_factor_ratio / thickness_;
			pos_n_[index_i] -= factor * constrained_distance_ * normal;
		}
		//=================================================================================================//
		void ShellMidSurfaceBounding::getNormalDirection()
		{
			parallel_for(blocked_range<size_t>(0, particles_->total_real_particles_),
				[&](const blocked_range<size_t>& r) {
				for (size_t i = r.begin(); i != r.end(); ++i)
				{
					n_0_[i] = level_set_shape_->findNormalDirection(pos_n_[i] + 0.2 * thickness_ * n_0_[i]);
					n_0_[i] /= n_0_[i].norm() + TinyReal;
				}
			}, ap);
		}
		//=================================================================================================//
		void ShellMidSurfaceBounding::setColorFunction()
		{
			parallel_for(blocked_range<size_t>(0, particles_->total_real_particles_),
				[&](const blocked_range<size_t>& r) {
				for (size_t i = r.begin(); i != r.end(); ++i)
				{
					color_[i] = 0;
				}
			}, ap);
			color_[0] = 1;
		}
		//=================================================================================================//
		void ShellMidSurfaceBounding::correctNormalDirection()
		{
			for (size_t index_i = 0; index_i < particles_->total_real_particles_; ++index_i)
			{
				Neighborhood& inner_neighborhood = inner_configuration_[index_i];
				for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
				{
					if (color_[index_i])
					{
						size_t index_j = inner_neighborhood.j_[n];
						if (!color_[index_j])
						{
							color_[index_j] = 1;
							Real included_angle_cosine = dot(n_0_[index_i], n_0_[index_j]);
							if (included_angle_cosine < 0)
							{
								n_0_[index_j] = -n_0_[index_j];
							}
						}
					}
				}
			}
		}
		//=================================================================================================//
		void ShellMidSurfaceBounding::averageNormalDirectionFirstHalf()
		{
			Real W0 = body_->sph_adaptation_->getKernel()->W0(Vecd(0));
			parallel_for(blocked_range<size_t>(0, particles_->total_real_particles_),
				[&](const blocked_range<size_t>& r) {
				for (size_t index_i = r.begin(); index_i != r.end(); ++index_i)
				{
					Vecd average_n = Vecd(0.0);
					Real weight = 0.0;
					Neighborhood& inner_neighborhood = inner_configuration_[index_i];
					for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
					{
						size_t index_j = inner_neighborhood.j_[n];
						average_n += inner_neighborhood.W_ij_[n] * n_0_[index_j];
						weight += inner_neighborhood.W_ij_[n];
					}
					average_n = average_n / (weight + TinyReal);
					temporary_n_0_[index_i] = average_n / (average_n.norm() + TinyReal);
				}
			}, ap);
		}
		//=================================================================================================//
		void ShellMidSurfaceBounding::averageNormalDirectionSecondHalf()
		{
			parallel_for(blocked_range<size_t>(0, particles_->total_real_particles_),
				[&](const blocked_range<size_t>& r) {
				for (size_t index_i = r.begin(); index_i != r.end(); ++index_i)
				{
					n_0_[index_i] = temporary_n_0_[index_i];
				}
			}, ap);
		}
		//=================================================================================================//
		void ShellMidSurfaceBounding::getDirectionCriteria()
		{
			direction_criteria_ = dot(previous_n_0_[0], n_0_[0]);
			parallel_for(blocked_range<size_t>(1, particles_->total_real_particles_),
				[&](const blocked_range<size_t>& r) {
				for (size_t index_i = r.begin(); index_i != r.end(); ++index_i)
				{
					direction_criteria_ = SMIN(direction_criteria_, dot(previous_n_0_[index_i], n_0_[index_i]));
				}
			}, ap);
		}
		//=================================================================================================//
		void ShellMidSurfaceBounding::getIncludedAngleCriteria()
		{
			included_angle_ = 1.0;
			parallel_for(blocked_range<size_t>(0, particles_->total_real_particles_),
				[&](const blocked_range<size_t>& r) {
				for (size_t index_i = r.begin(); index_i != r.end(); ++index_i)
				{
					Neighborhood& inner_neighborhood = inner_configuration_[index_i];
					for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
					{
						size_t index_j = inner_neighborhood.j_[n];
						included_angle_ = SMIN(included_angle_, dot(n_0_[index_i], n_0_[index_j]));
					}
				}
			}, ap);
		}
		//=================================================================================================//
		void ShellMidSurfaceBounding::calculateNormalDirection()
		{
			for (int i = 0; i < 40; ++i)
			{
				getNormalDirection();
				setColorFunction();
				for (int j = 0; j < 80; ++j)
				{
					correctNormalDirection();
					getIncludedAngleCriteria();
					if (included_angle_ >= 0) break;
				}
				for (int k = 0; k < 80; ++k)
				{
					parallel_for(blocked_range<size_t>(0, particles_->total_real_particles_),
						[&](const blocked_range<size_t>& r) {
						for (size_t index_i = r.begin(); index_i != r.end(); ++index_i)
						{
							previous_n_0_[index_i] = n_0_[index_i];
						}
					}, ap);
					getNormalDirection();
					averageNormalDirectionFirstHalf();
					averageNormalDirectionSecondHalf();
					getDirectionCriteria();
					if (direction_criteria_ > cos(0.01 * Pi)) break;
				}
				getIncludedAngleCriteria();
				if (included_angle_ > cos(0.01 * Pi)) break;
			}
		}
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
