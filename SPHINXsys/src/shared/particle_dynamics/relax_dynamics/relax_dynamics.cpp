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
			  h_ratio_(*particles_->getVariableByName<indexScalar, Real>("SmoothingLengthRatio")),
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
			  Vol_(particles_->Vol_), dvel_dt_(particles_->dvel_dt_)
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
			  real_body_(inner_relation.real_body_), inner_relation_(inner_relation),
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
	}
	//=================================================================================================//
}
//=================================================================================================//
