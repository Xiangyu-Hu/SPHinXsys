/**
 * @file 	relax_dynamics.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 * @version	0.1
 */

#include "relax_dynamics.h"
#include "level_set.h"

using namespace SimTK;
//=================================================================================================//
namespace SPH
{
//=================================================================================================//
	namespace relax_dynamics
	{
		//=================================================================================================//
		GetTimeStepSizeSquare::GetTimeStepSizeSquare(SPHBody* body) :
			ParticleDynamicsReduce<Real, ReduceMin>(body),
			RelaxDataDelegateSimple(body), dvel_dt_(particles_->dvel_dt_)
		{
			smoothing_length_ = body->kernel_->GetSmoothingLength();
			//given by sound criteria by assuming unit speed of sound and zero particle velocity
			initial_reference_ = 0.0625 * smoothing_length_ * smoothing_length_;
		}
		//=================================================================================================//
		Real GetTimeStepSizeSquare::ReduceFunction(size_t index_i, Real dt)
		{
			return 0.0625 * smoothing_length_ / (dvel_dt_[index_i].norm() + TinyReal);
		}
		//=================================================================================================//
		RelaxationAccelerationInner::RelaxationAccelerationInner(SPHBodyInnerRelation* body_inner_relation) : 
			InteractionDynamics(body_inner_relation->sph_body_), 
			RelaxDataDelegateInner(body_inner_relation),
			Vol_(particles_->Vol_), dvel_dt_(particles_->dvel_dt_), pos_n_(particles_->pos_n_) 
		{
			complex_shape_ = body_->body_shape_;
			kernel_ = body_->kernel_;
		}
		//=================================================================================================//
		void RelaxationAccelerationInner::Interaction(size_t index_i, Real dt)
		{
			Vecd acceleration(0);// = -2.0 * complex_shape_->computeKernelIntegral(pos_n_[index_i], kernel_);
			Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				acceleration -= 2.0 * inner_neighborhood.dW_ij_[n] * inner_neighborhood.e_ij_[n]* Vol_[index_j];
			}
			dvel_dt_[index_i] = acceleration;
		}
		//=================================================================================================//
		UpdateParticlePosition::UpdateParticlePosition(SPHBody* body) :
			ParticleDynamicsSimple(body), RelaxDataDelegateSimple(body), 
			pos_n_(particles_->pos_n_), dvel_dt_(particles_->dvel_dt_)	{}
		//=================================================================================================//
		void UpdateParticlePosition::Update(size_t index_i, Real dt_square)
		{
			pos_n_[index_i] += dvel_dt_[index_i] * dt_square * 0.5;
		}
		//=================================================================================================//
		RelaxationAccelerationComplex::
			RelaxationAccelerationComplex(SPHBodyComplexRelation* body_complex_relation) : 
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
				Neighborhood& contact_neighborhood = contact_configuration_[k][index_i];
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
		BodySurfaceBounding::
			BodySurfaceBounding(SPHBody* body, NearBodySurface* body_part) :
			PartDynamicsByCell(body, body_part),
			RelaxDataDelegateSimple(body), pos_n_(particles_->pos_n_)
		{
		}
		//=================================================================================================//
		void BodySurfaceBounding::Update(size_t index_i, Real dt)
		{
			Real phi = body_->body_shape_->findSignedDistance(pos_n_[index_i]);
			if (phi > -0.5 * body_->particle_spacing_)
			{
				Vecd unit_normal = body_->body_shape_->findNormalDirection(pos_n_[index_i]);
				unit_normal /= unit_normal.norm() + TinyReal;
				pos_n_[index_i] -= (phi + 0.5 * body_->particle_spacing_) * unit_normal;
			}
		}
		//=================================================================================================//
		ConstraintSurfaceParticles::
			ConstraintSurfaceParticles(SPHBody* body, BodySurface* body_part)
			:PartDynamicsByParticle(body, body_part),
			RelaxDataDelegateSimple(body), pos_n_(particles_->pos_n_)
		{
		}
		//=================================================================================================//
		void ConstraintSurfaceParticles::Update(size_t index_i, Real dt)
		{
			Real phi = body_->body_shape_->findSignedDistance(pos_n_[index_i]);
			Vecd unit_normal = body_->body_shape_->findNormalDirection(pos_n_[index_i]);
			unit_normal /= unit_normal.norm() + TinyReal;
			pos_n_[index_i] -= (phi + 0.5 * body_->particle_spacing_) * unit_normal;
		}
		//=================================================================================================//
		RelaxationStepInner::RelaxationStepInner(SPHBodyInnerRelation* body_inner_relation) :
			ParticleDynamics<void>(body_inner_relation->sph_body_),
			sph_body_(body_inner_relation->sph_body_), inner_relation_(body_inner_relation),
			relaxation_acceleration_inner_(inner_relation_),
			get_time_step_square_(sph_body_), update_particle_position_(sph_body_),
			surface_bounding_(sph_body_, new NearBodySurface(sph_body_)) {}
		//=================================================================================================//
		void RelaxationStepInner::exec(Real dt)
		{
			sph_body_->updateCellLinkedList();
			inner_relation_->updateConfiguration();
			relaxation_acceleration_inner_.exec();
			Real dt_square = get_time_step_square_.exec();
			update_particle_position_.exec(dt_square);
			surface_bounding_.exec();
		}
		//=================================================================================================//
		void RelaxationStepInner::parallel_exec(Real dt)
		{
			sph_body_->updateCellLinkedList();
			inner_relation_->updateConfiguration();
			relaxation_acceleration_inner_.parallel_exec();
			Real dt_square = get_time_step_square_.parallel_exec();
			update_particle_position_.parallel_exec(dt_square);
			surface_bounding_.parallel_exec();
		}
		//=================================================================================================//
	}
	//=================================================================================================//
}
//=================================================================================================//
