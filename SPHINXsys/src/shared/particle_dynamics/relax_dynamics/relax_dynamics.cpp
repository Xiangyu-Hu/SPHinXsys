/**
 * @file 	relax_dynamics.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 */

#include "relax_dynamics.h"
#include "particle_generator_lattice.h"
#include "geometry_level_set.h"

using namespace SimTK;
//=================================================================================================//
namespace SPH
{
//=================================================================================================//
	namespace relax_dynamics
	{
		//=================================================================================================//
		GetTimeStepSizeSquare::GetTimeStepSizeSquare(SPHBody* body) :
			ParticleDynamicsReduce<Real, ReduceMax>(body),
			RelaxDataDelegateSimple(body), dvel_dt_(particles_->dvel_dt_)
		{
			h_ref_ = body->particle_adaptation_->ReferenceSmoothingLength();
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
		RelaxationAccelerationInner::RelaxationAccelerationInner(BaseBodyRelationInner* body_inner_relation) : 
			InteractionDynamics(body_inner_relation->sph_body_), 
			RelaxDataDelegateInner(body_inner_relation),
			Vol_(particles_->Vol_), dvel_dt_(particles_->dvel_dt_), pos_n_(particles_->pos_n_) {}
		//=================================================================================================//
		void RelaxationAccelerationInner::Interaction(size_t index_i, Real dt)
		{
			Vecd acceleration(0);
			Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				acceleration -= 2.0 * inner_neighborhood.dW_ij_[n] * inner_neighborhood.e_ij_[n] * Vol_[index_j];
			}
			dvel_dt_[index_i] = acceleration;
		}
		//=================================================================================================//
		RelaxationAccelerationInnerWithLevelSetCorrection::
			RelaxationAccelerationInnerWithLevelSetCorrection(BaseBodyRelationInner* body_inner_relation) : 
			RelaxationAccelerationInner(body_inner_relation)
		{
			level_set_complex_shape_ = dynamic_cast<LevelSetComplexShape*>(body_->body_shape_);
			if (level_set_complex_shape_ == nullptr)
			{
				std::cout << "\n FAILURE: LevelSetComplexShape is undefined!" << std::endl;
				std::cout << __FILE__ << ':' << __LINE__ << std::endl;
				exit(1);
			}
		}		
		//=================================================================================================//
		void RelaxationAccelerationInnerWithLevelSetCorrection::Interaction(size_t index_i, Real dt)
		{
			RelaxationAccelerationInner::Interaction(index_i, dt);
			dvel_dt_[index_i] -= 2.0 * level_set_complex_shape_->
				computeKernelGradientIntegral(pos_n_[index_i], particle_adaptation_->SmoothingLengthRatio(index_i));
		}
		//=================================================================================================//
		UpdateParticlePosition::UpdateParticlePosition(SPHBody* body) :
			ParticleDynamicsSimple(body), RelaxDataDelegateSimple(body), 
			pos_n_(particles_->pos_n_), dvel_dt_(particles_->dvel_dt_) {}
		//=================================================================================================//
		void UpdateParticlePosition::Update(size_t index_i, Real dt_square)
		{
			pos_n_[index_i] += dvel_dt_[index_i] * dt_square * 0.5 / particle_adaptation_->SmoothingLengthRatio(index_i);
		}
		//=================================================================================================//
		UpdateSolidParticlePosition::UpdateSolidParticlePosition(SPHBody* body) :
			ParticleDynamicsSimple(body), solid_dynamics::SolidDataSimple(body),
			pos_0_(particles_->pos_0_), pos_n_(particles_->pos_n_), dvel_dt_(particles_->dvel_dt_) {}
		//=================================================================================================//
		void UpdateSolidParticlePosition::Update(size_t index_i, Real dt_square)
		{
			pos_n_[index_i] += dvel_dt_[index_i] * dt_square * 0.5 / particle_adaptation_->SmoothingLengthRatio(index_i);
			pos_0_[index_i] = pos_n_[index_i];
		}
		//=================================================================================================//
		UpdateSmoothingLengthRatioByBodyShape::UpdateSmoothingLengthRatioByBodyShape(SPHBody* body) :
			ParticleDynamicsSimple(body), RelaxDataDelegateSimple(body),
			h_ratio_(*particles_->getVariableByName<indexScalar, Real>("SmoothingLengthRatio")), 
			Vol_(particles_->Vol_), pos_n_(particles_->pos_n_),
			body_shape_(*body->body_shape_), kernel_(*body->particle_adaptation_->getKernel())
		{
			particle_spacing_by_body_shape_ =
				dynamic_cast<ParticleSpacingByBodyShape*>(body->particle_adaptation_);
		}
		//=================================================================================================//
		void UpdateSmoothingLengthRatioByBodyShape::Update(size_t index_i, Real dt_square)
		{
			Real local_spacing = particle_spacing_by_body_shape_->getLocalSpacing(body_shape_, pos_n_[index_i]);
			h_ratio_[index_i] = particle_adaptation_->ReferenceSpacing() / local_spacing;
			Vol_[index_i] = powerN(local_spacing, Dimensions);
		}
		//=================================================================================================//
		RelaxationAccelerationComplex::
			RelaxationAccelerationComplex(ComplexBodyRelation* body_complex_relation) : 
			InteractionDynamics(body_complex_relation->sph_body_),
			RelaxDataDelegateComplex(body_complex_relation), 
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
			Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				acceleration -= 2.0 * inner_neighborhood.dW_ij_[n] * inner_neighborhood.e_ij_[n]* Vol_[index_j];
			}

			/** Contact interaction. */
			for (size_t k = 0; k < contact_configuration_.size(); ++k)
			{
				StdLargeVec<Real>& Vol_k = *(contact_Vol_[k]);
				Neighborhood& contact_neighborhood = (*contact_configuration_[k])[index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					size_t index_j = contact_neighborhood.j_[n];

					acceleration -= 2.0 * Vol_k[index_j] 
								  * contact_neighborhood.dW_ij_[n] * contact_neighborhood.e_ij_[n];
				}
			}
			
			dvel_dt_[index_i] = acceleration;
		}
		//=================================================================================================//
		ShapeSurfaceBounding::
			ShapeSurfaceBounding(SPHBody* body, NearShapeSurface* body_part) :
			PartDynamicsByCell(body, body_part), RelaxDataDelegateSimple(body), 
			pos_n_(particles_->pos_n_), 
			constrained_distance_(0.5 * body->particle_adaptation_->MinimumSpacing()) {}
		//=================================================================================================//
		void ShapeSurfaceBounding::Update(size_t index_i, Real dt)
		{
			Real phi = body_->body_shape_->findSignedDistance(pos_n_[index_i]);

			if (phi > -constrained_distance_)
			{
				Vecd unit_normal = body_->body_shape_->findNormalDirection(pos_n_[index_i]);
				unit_normal /= unit_normal.norm() + TinyReal;
				pos_n_[index_i] -= (phi + constrained_distance_) * unit_normal;
			}
		}
		//=================================================================================================//
		ConstraintSurfaceParticles::
			ConstraintSurfaceParticles(SPHBody* body, ShapeSurface* body_part)
			:PartSimpleDynamicsByParticle(body, body_part), RelaxDataDelegateSimple(body),
			constrained_distance_(0.5 * body->particle_adaptation_->MinimumSpacing()),
			pos_n_(particles_->pos_n_)
		{
		}
		//=================================================================================================//
		void ConstraintSurfaceParticles::Update(size_t index_i, Real dt)
		{
			Real phi = body_->body_shape_->findSignedDistance(pos_n_[index_i]);
			Vecd unit_normal = body_->body_shape_->findNormalDirection(pos_n_[index_i]);
			unit_normal /= unit_normal.norm() + TinyReal;
			pos_n_[index_i] -= (phi + constrained_distance_) * unit_normal;
		}
		//=================================================================================================//
		RelaxationStepInner::
			RelaxationStepInner(BaseBodyRelationInner* body_inner_relation, bool level_set_correction) :
			ParticleDynamics<void>(body_inner_relation->sph_body_),
			real_body_(body_inner_relation->real_body_), inner_relation_(body_inner_relation),
			relaxation_acceleration_inner_(nullptr),
			get_time_step_square_(real_body_), update_particle_position_(real_body_),
			surface_bounding_(real_body_, new NearShapeSurface(real_body_))
		{
			relaxation_acceleration_inner_ = !level_set_correction ?
				new RelaxationAccelerationInner(body_inner_relation) :
				new RelaxationAccelerationInnerWithLevelSetCorrection(body_inner_relation);
		}
		//=================================================================================================//
		void RelaxationStepInner::exec(Real dt)
		{
			real_body_->updateCellLinkedList();
			inner_relation_->updateConfiguration();
			relaxation_acceleration_inner_->exec();
			Real dt_square = get_time_step_square_.exec();
			update_particle_position_.exec(dt_square);
			surface_bounding_.exec();
		}
		//=================================================================================================//
		void RelaxationStepInner::parallel_exec(Real dt)
		{
			real_body_->updateCellLinkedList();
			inner_relation_->updateConfiguration();
			relaxation_acceleration_inner_->parallel_exec();
			Real dt_square = get_time_step_square_.parallel_exec();
			update_particle_position_.parallel_exec(dt_square);
			surface_bounding_.parallel_exec();
		}
		//=================================================================================================//
		void SolidRelaxationStepInner::exec(Real dt)
		{
			real_body_->updateCellLinkedList();
			inner_relation_->updateConfiguration();
			relaxation_acceleration_inner_->exec();
			Real dt_square = get_time_step_square_.exec();
			update_solid_particle_position_.exec(dt_square);
			surface_bounding_.exec();
		}
		//=================================================================================================//
		void SolidRelaxationStepInner::parallel_exec(Real dt)
		{
			real_body_->updateCellLinkedList();
			inner_relation_->updateConfiguration();
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
