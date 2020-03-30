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
		GetTimeStepSize::GetTimeStepSize(RelaxBody* body)
			: RelaxDynamicsMin(body)
		{
			smoothing_length_ = body->kernel_->GetSmoothingLength();
			Real eta = 0.1 * 1.3 * body->particle_spacing_;
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
		PhysicsRelaxationInner::PhysicsRelaxationInner(RelaxBody* body)
			: RelaxDynamicsInner1Level(body)
		{
			p0_ = 1.0;
			sound_speed_ = 1.0;
			smoothing_length_ = body_->kernel_->GetSmoothingLength();
			eta_ = 0.1 * smoothing_length_;
			p_star_ = p0_ * sound_speed_ * sound_speed_;
			mass_ = 1.0 * powern(body_->particle_spacing_, Vecd(0).size());
			std::cout << "The background pressure for relaxation:" << " " << p_star_ << std::endl;
			std::cout << "The viscosity for relaxation:" << " " << eta_ << std::endl;
		}
		//=================================================================================================//
		void PhysicsRelaxationInner::Initialization(size_t index_particle_i, Real dt)
		{
			BaseParticleData& base_particle_data_i
				= particles_->base_particle_data_[index_particle_i];

			base_particle_data_i.pos_n_ += base_particle_data_i.vel_n_ * dt * 0.5;
		}
		//=================================================================================================//
		void PhysicsRelaxationInner::InnerInteraction(size_t index_particle_i, Real dt)
		{
			BaseParticleData& base_particle_data_i
				= particles_->base_particle_data_[index_particle_i];
			RelaxBodyParticleData& relax_data_i
				= particles_->relax_body_data_[index_particle_i];

			//including gravity and force from fluid
			Vecd acceleration(0);

			NeighborList& inner_neighors
				= getNeighborList(current_configuration_, index_particle_i);
			for (size_t n = 0; n < inner_neighors.size(); ++n)
			{
				BaseNeighborRelation* neighboring_particle = inner_neighors[n];
				size_t index_particle_j = neighboring_particle->j_;

				BaseParticleData& base_particle_data_j
					= particles_->base_particle_data_[index_particle_j];
				RelaxBodyParticleData& relax_data_j
					= particles_->relax_body_data_[index_particle_j];

				acceleration -= 2.0 * p_star_ * neighboring_particle->dW_ij_ * neighboring_particle->e_ij_
					* base_particle_data_j.Vol_ * base_particle_data_i.Vol_;
				//viscous force
				Vecd vel_detivative = (base_particle_data_i.vel_n_ - base_particle_data_j.vel_n_)
					/ (neighboring_particle->r_ij_ + 0.01 * smoothing_length_);
				acceleration += 0.5 * eta_ * (base_particle_data_i.vel_n_.norm() + base_particle_data_j.vel_n_.norm())
					* vel_detivative * base_particle_data_j.Vol_ * base_particle_data_i.Vol_ * neighboring_particle->dW_ij_;
				if (isnan(acceleration.norm()))
				{
					std::cout << "\n Error: the viscous acceleration is not valid" << std::endl;
					std::cout << "\n Particle ID :" << index_particle_i << "\n Particle ID :" << index_particle_j << std::endl;
					std::cout << "\n Particle distance :" << neighboring_particle->r_ij_ << std::endl;
					std::cout << __FILE__ << ':' << __LINE__ << std::endl;
					exit(1);
				}
			}

			Real phi = body_->mesh_background_->ProbeLevelSet(base_particle_data_i.pos_n_);
			//Vecd norm  = SGN(phi) *  elastic_body_.mesh_background_->ProbeNormalDirection(base_particle_data_i.pos_n_);
			Vecd dist_2_face(0);
			Vecd norm(0);
			dist_2_face = body_->mesh_background_->ProbeNormalDirection(base_particle_data_i.pos_n_);
			norm = -dist_2_face / (dist_2_face.norm() + 1.0e-15);

			if ((phi - 0.5 * body_->particle_spacing_) <= smoothing_length_)
			{
				Real q = (phi - 0.5 * body_->particle_spacing_) / smoothing_length_;
				if(Vecd(0).size() == 2)
				{
					Real gamma_ = 1.0 + (0.0625 - 0.0531 * q) * (q - 2.0) * (q - 2.0) * (q - 2.0);
					Real grad_gamma_ = (0.2937 - 0.2124 * q) * (q - 2.0) * (q - 2.0);
					acceleration += p_star_ * base_particle_data_i.Vol_ / smoothing_length_ / gamma_ * grad_gamma_ * norm;
				}else if(Vecd(0).size() == 3)
				{
					Real kernel_factor = 21.0 / (16.0 * pi) / smoothing_length_ / smoothing_length_ / smoothing_length_;
					Real grad_gamma = kernel_factor * 2.0 * 3.1415926 * (2.0 + 5.0 * q + 4.0 * q * q) * pow((1.0 - 0.5 * q), 5.0) / 7.0; 
					acceleration += p_star_ * base_particle_data_i.Vol_ * grad_gamma * norm;
				}

				if ((phi - 0.5 * body_->particle_spacing_) <= back_mesh_spacing_)
				{
					acceleration -= dot(acceleration, norm) * norm;
				}
			}
			
			base_particle_data_i.dvel_dt_ = acceleration / mass_;
		}
		//=================================================================================================//
		void PhysicsRelaxationInner::Update(size_t index_particle_i, Real dt)
		{
			BaseParticleData& base_particle_data_i
				= particles_->base_particle_data_[index_particle_i];
			RelaxBodyParticleData& relax_data_i
				= particles_->relax_body_data_[index_particle_i];

			base_particle_data_i.vel_n_ = base_particle_data_i.dvel_dt_ * dt;
			base_particle_data_i.pos_n_ += base_particle_data_i.vel_n_ * dt * 0.5;
		}
//=================================================================================================//
		PhysicsRelaxationComplex::PhysicsRelaxationComplex(RelaxBody *body, StdVec<RelaxBody*> interacting_bodies)
			: RelaxDynamicsComplex1Level(body, interacting_bodies)
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
			RelaxBodyParticleData& relax_data_i = particles_->relax_body_data_[index_particle_i];

			Vecd acceleration = base_particle_data_i.dvel_dt_others_;

			NeighborList& inner_neighors
				= getNeighborList(current_configuration_, index_particle_i);
			for (size_t n = 0; n < inner_neighors.size(); ++n)
			{
				BaseNeighborRelation* neighboring_particle = inner_neighors[n];
				size_t index_particle_j = neighboring_particle->j_;
				BaseParticleData& base_particle_data_j = particles_->base_particle_data_[index_particle_j];
				RelaxBodyParticleData& relax_data_j = particles_->relax_body_data_[index_particle_j];

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
				NeighborList& contact_neighors
					= getNeighborList(current_interacting_configuration_[k], index_particle_i);
				for (size_t n = 0; n < contact_neighors.size(); ++n)
				{
					BaseNeighborRelation* neighboring_particle = contact_neighors[n];
					size_t index_particle_j = neighboring_particle->j_;
					BaseParticleData& base_particle_data_j
						= (*interacting_particles_[k]).base_particle_data_[index_particle_j];
					RelaxBodyParticleData& relax_data_j
						= (*interacting_particles_[k]).relax_body_data_[index_particle_j];

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
		void ConstriantSurfaceParticles::ConstraintAParticle(size_t index_particle_i, Real dt)
		{
			BaseParticleData& base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			RelaxBodyParticleData& relax_particle_data_i = particles_->relax_body_data_[index_particle_i];

			Real phi = body_->mesh_background_->ProbeLevelSet(base_particle_data_i.pos_n_);
			Vecd dist_2_face = body_->mesh_background_->ProbeNormalDirection(base_particle_data_i.pos_n_);
			Vecd norm = dist_2_face / (dist_2_face.norm() + 1.0e-15);

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
		void ConstriantSingularityParticles::ConstraintAParticle(size_t index_particle_i, Real dt)
		{
			BaseParticleData& base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			RelaxBodyParticleData& relax_particle_data_i = particles_->relax_body_data_[index_particle_i];

			base_particle_data_i.pos_n_ = relax_particle_data_i.pos_0_;
		}
		//=================================================================================================//
		void computeNumberDensitybySummation::ComplexInteraction(size_t index_particle_i, Real dt)
		{
			BaseParticleData& base_particle_data_i = particles_->base_particle_data_[index_particle_i];
			Real Vol_0_i = base_particle_data_i.Vol_0_;

			/** Inner interaction. */
			Real sigma = W0_;
			NeighborList& inner_neighors
				= getNeighborList(current_configuration_, index_particle_i);
			for (size_t n = 0; n < inner_neighors.size(); ++n)
			{
				BaseNeighborRelation* neighboring_particle = inner_neighors[n];

				sigma += neighboring_particle->W_ij_;
			}

			/** Contact interaction. */
			for (size_t k = 0; k < current_interacting_configuration_.size(); ++k)
			{
				NeighborList& contact_neighors
					= getNeighborList(current_interacting_configuration_[k], index_particle_i);
				for (size_t n = 0; n < contact_neighors.size(); ++n)
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
		getAveragedParticleNumberDensity::getAveragedParticleNumberDensity(RelaxBody* body)
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
		void updateNumberDensity::Update(size_t index_particle_i, Real dt)
		{
			BaseParticleData& base_particle_data_i = particles_->base_particle_data_[index_particle_i];

			base_particle_data_i.sigma_0_ = sigma_;
		}
		//=================================================================================================//
	}
	//=================================================================================================//
}
//=================================================================================================//
