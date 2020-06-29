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
			initial_reference_ = 0.125 * smoothing_length_ * smoothing_length_ * 1.0 / (eta + 1.0e-15);
		}
		//=================================================================================================//
		Real GetTimeStepSize::ReduceFunction(size_t index_particle_i, Real dt)
		{
			BaseParticleData& base_particle_data_i
				= particles_->base_particle_data_[index_particle_i];

			return 0.25 * SMIN(sqrt(smoothing_length_ / (base_particle_data_i.dvel_dt_.norm() + 1.0e-15)),
				smoothing_length_ / (40.0 * base_particle_data_i.vel_n_.norm() + 1.0e-15));
		}
		//=================================================================================================//
		PhysicsRelaxationInner::PhysicsRelaxationInner(SPHBody* body) 
			: ParticleDynamicsInner<SPHBody>(body),
			p0_(1.0), sound_speed_(1.0), mesh_background_(body->mesh_background_),
			particle_spacing_(body->particle_spacing_)
		{
			smoothing_length_ = body_->kernel_->GetSmoothingLength();
			p_star_ = p0_ * sound_speed_ * sound_speed_;
			back_mesh_spacing_ = mesh_background_->getGridSpacing();
			mass_ = 1.0 * powern(body_->particle_spacing_, Vecd(0).size());
			std::cout << "The background pressure for relaxation:" << " " << p_star_ << std::endl;
		};
			//=================================================================================================//
		void PhysicsRelaxationInner::InnerInteraction(size_t index_particle_i, Real dt)
		{
			BaseParticleData& base_particle_data_i = particles_->base_particle_data_[index_particle_i];

			Vecd acceleration(0);
			Neighborhood& inner_neighborhood = (*inner_configuration_)[index_particle_i];
			NeighborList& inner_neighors = std::get<0>(inner_neighborhood);
			for (size_t n = 0; n != std::get<2>(inner_neighborhood); ++n)
			{
				BaseNeighborRelation* neighboring_particle = inner_neighors[n];
				size_t index_particle_j = neighboring_particle->j_;
				BaseParticleData& base_particle_data_j = particles_->base_particle_data_[index_particle_j];

				acceleration -= 2.0 * p_star_ * neighboring_particle->dW_ij_ * neighboring_particle->e_ij_
					* base_particle_data_j.Vol_ * base_particle_data_i.Vol_;
			}
			base_particle_data_i.dvel_dt_ = acceleration / mass_;
			base_particle_data_i.pos_n_ += acceleration / mass_* dt * dt * 0.5;
		}
		//=================================================================================================//
		PhysicsRelaxationComplex::PhysicsRelaxationComplex(SPHBody *body, StdVec<SPHBody*> interacting_bodies)
			: RelaxDynamicsComplex1Level<BaseParticles>(body, interacting_bodies)
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
			Neighborhood& inner_neighborhood = (*inner_configuration_)[index_particle_i];
			NeighborList& inner_neighors = std::get<0>(inner_neighborhood);
			for (size_t n = 0; n != std::get<2>(inner_neighborhood); ++n)
			{
				BaseNeighborRelation* neighboring_particle = inner_neighors[n];
				size_t index_particle_j = neighboring_particle->j_;
				BaseParticleData& base_particle_data_j = particles_->base_particle_data_[index_particle_j];

				acceleration -= 2.0 * p_star_ * base_particle_data_j.Vol_ * base_particle_data_i.Vol_
					* neighboring_particle->dW_ij_ * neighboring_particle->e_ij_;
				//viscous force
				Vecd vel_detivative = (base_particle_data_i.vel_n_ - base_particle_data_j.vel_n_)
					/ neighboring_particle->r_ij_;
				acceleration += 0.5 * eta_ * (base_particle_data_i.vel_n_.norm() + base_particle_data_j.vel_n_.norm())
					* vel_detivative * base_particle_data_j.Vol_ * base_particle_data_i.Vol_ * neighboring_particle->dW_ij_;
			}

			/** Contact interaction. */
			for (size_t k = 0; k < current_interacting_configuration_.size(); ++k)
			{
				Neighborhood& contact_neighborhood = (*current_interacting_configuration_[k])[index_particle_i];
				NeighborList& contact_neighors = std::get<0>(contact_neighborhood);
				for (size_t n = 0; n != std::get<2>(contact_neighborhood); ++n)
				{
					BaseNeighborRelation* neighboring_particle = contact_neighors[n];
					size_t index_particle_j = neighboring_particle->j_;
					BaseParticleData& base_particle_data_j
						= (*interacting_particles_[k]).base_particle_data_[index_particle_j];

					//pressure force
					acceleration -= 2.0 * p_star_ * base_particle_data_j.Vol_ * base_particle_data_i.Vol_
						* neighboring_particle->dW_ij_ * neighboring_particle->e_ij_;

					//viscous force
					Vecd vel_detivative = 2.0 * base_particle_data_i.vel_n_ / neighboring_particle->r_ij_;
					acceleration += 0.5 * eta_ * (base_particle_data_i.vel_n_.norm() + base_particle_data_j.vel_n_.norm())
						* vel_detivative * base_particle_data_j.Vol_ * base_particle_data_i.Vol_ * neighboring_particle->dW_ij_;
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
		void BodySurfaceBounding::ConstraintAParticle(size_t index_particle_i, Real dt)
		{
			BaseParticleData& base_particle_data_i = particles_->base_particle_data_[index_particle_i];

			Real phi = body_->mesh_background_->ProbeLevelSet(base_particle_data_i.pos_n_);
			if (phi > -0.5 * body_->particle_spacing_)
			{
				Vecd unit_normal = body_->mesh_background_->ProbeNormalDirection(base_particle_data_i.pos_n_);
				unit_normal /= unit_normal.norm() + 1.0e-16;
				base_particle_data_i.pos_n_ -= (phi + 0.5 * body_->particle_spacing_) * unit_normal;
			}
		}
		//=================================================================================================//
		void ConstriantSurfaceParticles::ConstraintAParticle(size_t index_particle_i, Real dt)
		{
			BaseParticleData& base_particle_data_i = particles_->base_particle_data_[index_particle_i];

			Real phi = body_->mesh_background_->ProbeLevelSet(base_particle_data_i.pos_n_);
			Vecd norm = body_->mesh_background_->ProbeNormalDirection(base_particle_data_i.pos_n_);

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
		void computeNumberDensitybySummation::ComplexInteraction(size_t index_particle_i, Real dt)
		{
			BaseParticleData& base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			Real Vol_0_i = base_particle_data_i.Vol_0_;

			/** Inner interaction. */
			Real sigma = W0_;
			Neighborhood& inner_neighborhood = (*inner_configuration_)[index_particle_i];
			NeighborList& inner_neighors = std::get<0>(inner_neighborhood);
			for (size_t n = 0; n != std::get<2>(inner_neighborhood); ++n)
			{
				BaseNeighborRelation* neighboring_particle = inner_neighors[n];

				sigma += neighboring_particle->W_ij_;
			}

			/** Contact interaction. */
			for (size_t k = 0; k < current_interacting_configuration_.size(); ++k)
			{
				Neighborhood& contact_neighborhood = (*current_interacting_configuration_[k])[index_particle_i];
				NeighborList& contact_neighors = std::get<0>(contact_neighborhood);
				for (size_t n = 0; n != std::get<2>(contact_neighborhood); ++n)
				{
					BaseNeighborRelation* neighboring_particle = contact_neighors[n];
					size_t index_particle_j = neighboring_particle->j_;
					BaseParticleData& base_particle_data_j
						= (*interacting_particles_[k]).base_particle_data_[index_particle_j];

					sigma += neighboring_particle->W_ij_ * base_particle_data_j.Vol_0_ / Vol_0_i;
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
