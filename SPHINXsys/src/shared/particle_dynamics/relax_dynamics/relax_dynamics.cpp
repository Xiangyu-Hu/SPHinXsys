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
			ParticleDynamicsInner(body_inner_relation), RelaxDataDelegateInner(body_inner_relation),
			Vol_(particles_->Vol_), dvel_dt_(particles_->dvel_dt_), pos_n_(particles_->pos_n_) 
		{
			complex_shape_ = body_->body_shape_;
			kernel_ = body_->kernel_;
		}
		//=================================================================================================//
		void RelaxationAccelerationInner::InnerInteraction(size_t index_i, Real dt)
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
			ParticleDynamicsComplex(body_complex_relation), 
			RelaxDataDelegateComplex(body_complex_relation), 
			Vol_(particles_->Vol_), dvel_dt_(particles_->dvel_dt_)
		{
			for (size_t k = 0; k != contact_particles_.size(); ++k)
			{
				contact_Vol_.push_back(&(contact_particles_[k]->Vol_));
			}
		}
		//=================================================================================================//
		void RelaxationAccelerationComplex::ComplexInteraction(size_t index_i, Real dt)
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
		computeNumberDensityBySummation::
			computeNumberDensityBySummation(SPHBodyComplexRelation* body_complex_relation)
			: ParticleDynamicsComplex(body_complex_relation), 
			RelaxDataDelegateComplex(body_complex_relation), Vol_0_(particles_->Vol_0_),
			sigma_0_(particles_->sigma_0_)
		{
			for (size_t k = 0; k != contact_particles_.size(); ++k)
			{
				contact_Vol_0_.push_back(&(contact_particles_[k]->Vol_0_));
			}
			W0_ = body_->kernel_->W(Vecd(0));
		}
		//=================================================================================================//
		void computeNumberDensityBySummation::ComplexInteraction(size_t index_i, Real dt)
		{
			Real Vol_0_i = Vol_0_[index_i];

			/** Inner interaction. */
			Real sigma = W0_;
			Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
				sigma += inner_neighborhood.W_ij_[n];

			/** Contact interaction. */
			for (size_t k = 0; k < contact_configuration_.size(); ++k)
			{
				StdLargeVec<Real>& Vol_0_k = *(contact_Vol_0_[k]);
				Neighborhood& contact_neighborhood = contact_configuration_[k][index_i];
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					size_t index_j = contact_neighborhood.j_[n];

					sigma += contact_neighborhood.W_ij_[n] * Vol_0_k[index_j] / Vol_0_i;
				}
			}

			/** Particle summation. */
			sigma_0_[index_i] = sigma;
		}
		//=================================================================================================//
		getAveragedParticleNumberDensity::getAveragedParticleNumberDensity(SPHBody* body)
			: ParticleDynamicsReduce <Real, ReduceSum<Real>>(body),
			RelaxDataDelegateSimple(body), sigma_0_(particles_->sigma_0_)
		{
			initial_reference_ = 0.0;
			average_farctor_ = 1.0 / Real(body_->number_of_particles_);
		}
		//=================================================================================================//
		Real getAveragedParticleNumberDensity::ReduceFunction(size_t index_i, Real dt)
		{
			return average_farctor_ * sigma_0_[index_i];
		}
		//=================================================================================================//
		FinalizingParticleRelaxation::
			FinalizingParticleRelaxation(SPHBody* body) : 
			ParticleDynamicsSimple(body), RelaxDataDelegateSimple(body), 
			sigma_0_(particles_->sigma_0_), pos_n_(particles_->pos_n_), sigma_(0.0)
		{
			get_average_number_density_ = new getAveragedParticleNumberDensity(body);
		}
		//=================================================================================================//
		void FinalizingParticleRelaxation::Update(size_t index_i, Real dt)
		{
			sigma_0_[index_i] = sigma_;
		}
		//=================================================================================================//
	}
	//=================================================================================================//
}
//=================================================================================================//
