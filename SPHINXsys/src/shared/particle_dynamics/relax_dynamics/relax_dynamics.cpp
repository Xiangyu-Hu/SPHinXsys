/**
 * @file 	relax_dynamics.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 * @version	0.1
 */

#include "relax_dynamics.h"

using namespace SimTK;
//=================================================================================================//
namespace SPH
{
//=================================================================================================//
	namespace relax_dynamics
	{
		//=================================================================================================//
		GetTimeStepSize::GetTimeStepSize(SPHBody* body) : RelaxDynamicsMin<BaseParticles>(body)
		{
			smoothing_length_ = body->kernel_->GetSmoothingLength();
			Real eta = 0.1 * smoothing_length_;
			initial_reference_ = 0.125 * smoothing_length_ * smoothing_length_ * 1.0 / (eta + TinyReal);
		}
		//=================================================================================================//
		Real GetTimeStepSize::ReduceFunction(size_t index_particle_i, Real dt)
		{
			BaseParticleData& base_particle_data_i
				= particles_->base_particle_data_[index_particle_i];

			return 0.25 * SMIN(sqrt(smoothing_length_ / (base_particle_data_i.dvel_dt_.norm() + TinyReal)),
				smoothing_length_ / (40.0 * base_particle_data_i.vel_n_.norm() + TinyReal));
		}
		//=================================================================================================//
		PhysicsRelaxationInner::PhysicsRelaxationInner(SPHBodyInnerRelation* body_inner_relation)
			: ParticleDynamicsInner<SPHBody>(body_inner_relation),
			p0_(1.0), sound_speed_(1.0), particle_spacing_(body_->particle_spacing_)
		{
			smoothing_length_ = body_->kernel_->GetSmoothingLength();
			p_star_ = p0_ * sound_speed_ * sound_speed_;
			mass_ = 1.0 * powern(body_->particle_spacing_, Vecd(0).size());
			std::cout << "The background pressure for relaxation:" << " " << p_star_ << std::endl;
		}
		//=================================================================================================//
		void PhysicsRelaxationInner::InnerInteraction(size_t index_particle_i, Real dt)
		{
			BaseParticleData& base_particle_data_i = particles_->base_particle_data_[index_particle_i];

			Vecd acceleration(0);
			Neighborhood& inner_neighborhood = inner_configuration_[index_particle_i];
			CommonRelationList& inner_common_relations = inner_neighborhood.common_relation_list_;
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				CommonRelation& common_relation = inner_common_relations[n];
				size_t index_particle_j = common_relation.j_;
				BaseParticleData& base_particle_data_j = particles_->base_particle_data_[index_particle_j];

				acceleration -= 2.0 * p_star_ * common_relation.dW_ij_ * common_relation.e_ij_
					* base_particle_data_j.Vol_ * base_particle_data_i.Vol_;
			}
			base_particle_data_i.dvel_dt_ = acceleration / mass_;
			base_particle_data_i.pos_n_ += acceleration / mass_* dt * dt * 0.5;
		}
		//=================================================================================================//
		PhysicsRelaxationComplex::PhysicsRelaxationComplex(SPHBodyComplexRelation* body_complex_relation)
			: RelaxDynamicsComplex1Level<BaseParticles>(body_complex_relation)
		{
			p0_ = 1.0;
			sound_speed_ = 1.0;
			eta_ = 0.1 * body_->kernel_->GetSmoothingLength();
			p_star_ = p0_ * sound_speed_ * sound_speed_;
			mass_ = 1.0 * powern(body_->particle_spacing_, Vecd(0).size());
			std::cout << "The background pressure for relaxation:" << " " << p_star_ << std::endl;
			std::cout << "The viscosity for relaxation:" << " " << eta_ << std::endl;
		}
		//=================================================================================================//
		void PhysicsRelaxationComplex::Initialization(size_t index_particle_i, Real dt)
		{
			BaseParticleData& base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			base_particle_data_i.pos_n_ += base_particle_data_i.vel_n_ * dt * 0.5;
		}
		//=================================================================================================//
		void PhysicsRelaxationComplex::ComplexInteraction(size_t index_particle_i, Real dt)
		{
			BaseParticleData& base_particle_data_i = particles_->base_particle_data_[index_particle_i];

			/** inner interaction*/
			Vecd acceleration = base_particle_data_i.dvel_dt_others_;
			Neighborhood& inner_neighborhood = inner_configuration_[index_particle_i];
			CommonRelationList& inner_common_relations = inner_neighborhood.common_relation_list_;
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				CommonRelation& common_relation = inner_common_relations[n];
				size_t index_particle_j = common_relation.j_;
				BaseParticleData& base_particle_data_j = particles_->base_particle_data_[index_particle_j];

				acceleration -= 2.0 * p_star_ * base_particle_data_j.Vol_ * base_particle_data_i.Vol_
					* common_relation.dW_ij_ * common_relation.e_ij_;
				//viscous force
				Vecd vel_detivative = (base_particle_data_i.vel_n_ - base_particle_data_j.vel_n_)
					/ common_relation.r_ij_;
				acceleration += 0.5 * eta_ * (base_particle_data_i.vel_n_.norm() + base_particle_data_j.vel_n_.norm())
					* vel_detivative * base_particle_data_j.Vol_ * base_particle_data_i.Vol_ * common_relation.dW_ij_;
			}

			/** Contact interaction. */
			for (size_t k = 0; k < contact_configuration_.size(); ++k)
			{
				Neighborhood& contact_neighborhood = contact_configuration_[k][index_particle_i];
				CommonRelationList& contact_common_relations = contact_neighborhood.common_relation_list_;
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					CommonRelation& common_relation = contact_common_relations[n];
					size_t index_particle_j = common_relation.j_;
					BaseParticleData& base_particle_data_j
						= contact_particles_[k]->base_particle_data_[index_particle_j];

					//pressure force
					acceleration -= 2.0 * p_star_ * base_particle_data_j.Vol_ * base_particle_data_i.Vol_
						* common_relation.dW_ij_ * common_relation.e_ij_;

					//viscous force
					Vecd vel_detivative = 2.0 * base_particle_data_i.vel_n_ / common_relation.r_ij_;
					acceleration += 0.5 * eta_ * (base_particle_data_i.vel_n_.norm() + base_particle_data_j.vel_n_.norm())
						* vel_detivative * base_particle_data_j.Vol_ * base_particle_data_i.Vol_ * common_relation.dW_ij_;
				}
			}
			base_particle_data_i.dvel_dt_ = acceleration / mass_;
		}
		//=================================================================================================//
		void PhysicsRelaxationComplex::Update(size_t index_particle_i, Real dt)
		{
			BaseParticleData& base_particle_data_i = particles_->base_particle_data_[index_particle_i];

			base_particle_data_i.vel_n_ = base_particle_data_i.dvel_dt_ * dt;
			base_particle_data_i.pos_n_ += base_particle_data_i.vel_n_ * dt * 0.5;
		}
		//=================================================================================================//
		void BodySurfaceBounding::Update(size_t index_particle_i, Real dt)
		{
			BaseParticleData& base_particle_data_i = particles_->base_particle_data_[index_particle_i];

			Real phi = body_->levelset_mesh_->probeLevelSet(base_particle_data_i.pos_n_);
			if (phi > -0.5 * body_->particle_spacing_)
			{
				Vecd unit_normal = body_->levelset_mesh_->probeNormalDirection(base_particle_data_i.pos_n_);
				unit_normal /= unit_normal.norm() + TinyReal;
				base_particle_data_i.pos_n_ -= (phi + 0.5 * body_->particle_spacing_) * unit_normal;
			}
		}
		//=================================================================================================//
		void ConstraintSurfaceParticles::Update(size_t index_particle_i, Real dt)
		{
			BaseParticleData& base_particle_data_i = particles_->base_particle_data_[index_particle_i];

			Real phi = body_->levelset_mesh_->probeLevelSet(base_particle_data_i.pos_n_);
			Vecd norm = body_->levelset_mesh_->probeNormalDirection(base_particle_data_i.pos_n_);

			if (phi >= 0.0)
			{
				base_particle_data_i.pos_n_ += (phi - 0.5 * body_->particle_spacing_) * norm;
			}
			else
			{
				base_particle_data_i.pos_n_ += (ABS(phi) + 0.5 * body_->particle_spacing_) * norm;
			}
		}
		//=================================================================================================//
		void computeNumberDensityBySummation::ComplexInteraction(size_t index_particle_i, Real dt)
		{
			BaseParticleData& base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			Real Vol_0_i = base_particle_data_i.Vol_0_;

			/** Inner interaction. */
			Real sigma = W0_;
			Neighborhood& inner_neighborhood = inner_configuration_[index_particle_i];
			KernelValueList& kernel_value_list = inner_neighborhood.kernel_value_list_;
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
				sigma += kernel_value_list[n];

			/** Contact interaction. */
			for (size_t k = 0; k < contact_configuration_.size(); ++k)
			{
				Neighborhood& contact_neighborhood = contact_configuration_[k][index_particle_i];
				KernelValueList& contact_kernel_values = contact_neighborhood.kernel_value_list_;
				CommonRelationList& contact_common_relations = contact_neighborhood.common_relation_list_;
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					CommonRelation& common_relation = contact_common_relations[n];
					size_t index_particle_j = common_relation.j_;
					BaseParticleData& base_particle_data_j
						= contact_particles_[k]->base_particle_data_[index_particle_j];

					sigma += contact_kernel_values[n] * base_particle_data_j.Vol_0_ / Vol_0_i;
				}
			}

			/** Particle summation. */
			base_particle_data_i.sigma_0_ = sigma;
		}
		//=================================================================================================//
		getAveragedParticleNumberDensity::getAveragedParticleNumberDensity(SPHBody* body)
			: RelaxDynamicsSum<Real>(body)
		{
			initial_reference_ = 0.0;
			average_farctor_ = 1.0 / Real(body_->number_of_particles_);
		}
		//=================================================================================================//
		Real getAveragedParticleNumberDensity::ReduceFunction(size_t index_particle_i, Real dt)
		{
			BaseParticleData& base_particle_data_i = particles_->base_particle_data_[index_particle_i];

			return average_farctor_ * (base_particle_data_i.sigma_0_);
		}
		//=================================================================================================//
		void FinalizingParticleRelaxation::Update(size_t index_particle_i, Real dt)
		{
			BaseParticleData& base_particle_data_i = particles_->base_particle_data_[index_particle_i];

			base_particle_data_i.sigma_0_ = sigma_;
			base_particle_data_i.pos_0_ = base_particle_data_i.pos_n_;
		}
		//=================================================================================================//
	}
	//=================================================================================================//
}
//=================================================================================================//
